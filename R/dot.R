#' An S4 class to represent a DOT object.
#'
#' @slot srt A list containing spatial data as processed by setup.srt
#' @slot ref A list containing reference data as processed by setup.ref
#' @slot weights A matrix containing absolute abundance of categories
#' @slot solution A matrix containing raw solution
#' @slot history A data.frame containing solution history
#' @export
Dot <- setClass("Dot", slots = list(srt = "list", ref = "list",
                                    weights = "matrix", solution = "matrix", history = "data.frame"))

#' Processing the reference single-cell data
#'
#' @param ref_data A gene x cell matrix of gene expressions. Can be a matrix-like object or a Seurat/AnnData object.
#' @param ref_annotations A character vector (one for each cell) or a single vector pointing to the slot in the Seurat/AnnData object
#' @param ref_subcluster_size An integer. Maximum number of sub-clusters per sub-population.
#' @param max_genes An integer. Maximum number of genes to pick.
#' @param remove_mt Boolean. Whether mitochondrial genes must be removed.
#' @param verbose Boolean. Whether progress should be displayed.
#' @return A list containing the processed ref data.
#' @export
#'
#' @examples
#' data(dot.sample)
#' dot.ref <- setup.ref(dot.sample$ref$counts[, 1:1000], dot.sample$ref$labels[1:1000], 2)
setup.ref <- function(ref_data, ref_annotations = NULL, ref_subcluster_size = 10,
                      max_genes = 5000, remove_mt = TRUE, verbose = FALSE)
{
  if(methods::is(ref_data, "Seurat"))
  {
    if(is.null(ref_annotations))
    {
      ref_annotations <- Seurat::Idents(ref_data)
    }else if(is.character(ref_annotations) && length(ref_annotations) == 1)
    {
      ref_annotations <- ref_data[[ref_annotations]]
    }

    # cell x gene
    ref_data <- t(Seurat::GetAssayData(ref_data, slot = "counts"))
  }else if(methods::is(ref_data, "AnnDataR6"))
  {
    if(is.null(ref_annotations))
    {
      for (o in c("cluster", "cell_type", "cell_subclass", "cell_class",
                  "cluster", "cell_types", "cell_subclasses", "cell_classes")) {
        if(!is.null(ref_data$obs[[o]]))
        {
          message(sprintf("Picking %s as the reference annotations.", o))
          ref_annotations <- ref_data$obs[[o]]
          break
        }
      }
    }else if(is.character(ref_annotations) && length(ref_annotations) == 1)
    {
      ref_annotations <- ref_data$obs[[ref_annotations]]
    }

    # X is already cell X gene in AnnData6
    ref_data <- ref_data$X
  }else if(methods::is(ref_data, "data.frame"))
  {
    # cell x gene
    ref_data <- t(as.matrix(ref_data))
  }else if(methods::is(ref_data, "matrix") | methods::is(ref_data, "Matrix"))
  {
    # cell x gene
    ref_data <- t(ref_data)
  }else
  {
    stop("Unsupported data type for ref_data")
  }

  if(is.null(ref_annotations))
    stop("Annotations must be supplied.")

  if(methods::is(ref_annotations, "data.frame") || methods::is(ref_annotations, "matrix"))
  {
    if(ncol(ref_annotations) > 1)
    {
      warning("ref_annotations contains more than one column. Picking the first column.")
      ref_annotations <- ref_annotations[, 1]
    }
  }

  ref_annotations <- as.character(ref_annotations)

  if(length(ref_annotations) != nrow(ref_data))
  {
    stop("Inconsistent number of annotations and cells in the reference data.")
  }

  if(verbose)
    message("Pre-processing")

  if(remove_mt)
  {
    mt <- is_mt(colnames(ref_data))
    ref_data <- ref_data[, which(!mt)]
  }

  vg_genes <- max(5000, max_genes)
  if(ncol(ref_data) > vg_genes)
  {
    vg <- Seurat::FindVariableFeatures(t(ref_data))
    if("vst.variance.standardized" %in% colnames(vg))
    {
      vg <- rownames(vg)[order(vg$vst.variance.standardized, decreasing = TRUE)[1:vg_genes]]
    }else{
      vg <- rownames(vg)[order(vg$variance.standardized, decreasing = TRUE)[1:vg_genes]]
    }

    ref_data <- ref_data[, vg]

    gc()
  }

  gc()

  nonempty_cells <- which(rowSums(ref_data) > 0)
  if(length(nonempty_cells) < nrow(ref_data))
  {
    ref_data <- ref_data[nonempty_cells, ]
    ref_annotations <- ref_annotations[nonempty_cells]
  }

  if(verbose)
    message("Aggregating")

  ref_agg <-  aggregate_ref(ref_data, ref_annotations, ref_subcluster_size, verbose = verbose)

  rm(ref_data)

  de_info <- get_de_genes(ref_agg$sub_centroids, ref_agg$sub_ratios, max_genes, verbose = verbose)

  r <- list(X = as.matrix(ref_agg$sub_centroids[, de_info$de_genes]), C = ref_agg$clusters, R = ref_agg$major_ratios)

  gc()

  return(r)
}


#' Processing the target spatial data
#'
#' @param srt_data A gene x spot/cell matrix of gene expressions. Can be a matrix-like object or a Seurat/AnnData object
#' @param srt_coords A matrix-like object with two columns (x and y coordinates of spots). If set to NULL and srt_data is Seurat/AnnData object, coords will be extracted from srt_data
#' @param th.spatial A value between 0 and 1. Threshold on similarity of adjacent spots
#' @param th.nonspatial A value between 0 and 1. Threshold on similarity of non-adjacent spots
#' @param th.gene.low Minimum percentage of spots that a valid gene must be expressed in.
#' @param th.gene.high Maximum percentage of spots that a valid gene must be expressed in.
#' @param remove_mt Boolean. Whether mitochondrial genes must be removed.
#' @param radius Adjacency radius. If set to 'auto' it is computed using the coordinates of spots
#' @param verbose Boolean. Whether progress should be displayed.
#' @return A list containing the processed srt data.
#' @export
#'
#' @examples
#' data(dot.sample)
#' dot.srt <- setup.srt(dot.sample$srt$counts, dot.sample$srt$coordinates)
setup.srt <- function(srt_data, srt_coords = NULL, th.spatial = 0.84, th.nonspatial = 0,
                      th.gene.low = 0.01, th.gene.high = 0.99, remove_mt = TRUE, radius = 'auto', verbose = FALSE)
{
  if(methods::is(srt_data, "Seurat"))
  {
    if(is.null(srt_coords))
    {
      srt_coords <- tryCatch(expr = Seurat::GetTissueCoordinates(srt_data),
                             error = function(e) {
                               message("Spatial coordinates could not be extracted from the Seurat object.")
                             })
    }

    # spot x gene
    srt_data <- t(Seurat::GetAssayData(srt_data, slot = "counts"))
  }else if(methods::is(srt_data, "AnnDataR6"))
  {
    if(is.null(srt_coords))
    {
      srt_coords <- srt_data$obsm$spatial
    }
    # X is already spot X gene in AnnData6
    srt_data <- srt_data$X
  }else if(methods::is(srt_data, "data.frame"))
  {
    # spot x gene
    srt_data <- t(as.matrix(srt_data))
  }else if(methods::is(srt_data, "matrix") | methods::is(srt_data, "Matrix"))
  {
    # spot x gene
    srt_data <- t(srt_data)
  }else
  {
    stop("Unsupported data type for srt_data")
  }

  if(!is.null(srt_coords))
  {
    if(!methods::is(srt_coords, "data.frame") && !methods::is(srt_coords, "matrix"))
      stop("Invalid coordinates supplied.")

    if(ncol(srt_coords) == 1)
      stop("Invalid coordinates supplied.")

    if(nrow(srt_coords) != nrow(srt_data))
      stop("Invalid coordinates supplied.")

    if(all(c("x", "y") %in% colnames(srt_coords)))
      srt_coords <- srt_coords[, c("x", "y")]
    else if(all(c("row", "col") %in% colnames(srt_coords)))
      srt_coords <- srt_coords[, c("row", "col")]
    else if(all(c("image_row", "image_col") %in% colnames(srt_coords)))
      srt_coords <- srt_coords[, c("image_row", "image_col")]
    else
    {
      if(ncol(srt_coords) > 2)
        warning("Picking the first two columns of coords as the spatial coordinates.")
      srt_coords <- srt_coords[, c(1,2)]
    }

    colnames(srt_coords) <- c("x", "y")
  }

  if(remove_mt)
  {
    mt <- is_mt(colnames(srt_data))
    srt_data <- as.matrix(srt_data[, which(!mt)])
  }

  if(th.gene.high < 1 | th.gene.low > 0)
  {
    cm <- colMeans(srt_data > 0)
    srt_data <- srt_data[, which(cm > th.gene.low & cm < th.gene.high)]
  }

  s <- list(X = srt_data, C = srt_coords)

  rm(srt_data)

  if(th.spatial > 0)
  {
    s$P <- get_pairs(s, th.spatial, th.nonspatial, nrow(s$X), verbose = verbose)
  }

  gc()

  return(s)
}


#' Creating a DOT object based on the processed ref and srt data
#'
#' @param srt A list containing the processed srt data produced by setup.srt
#' @param ref A list containing the processed ref data produced by setup.ref
#' @param ls_solution Boolean. Whether an initial solution based on LS should be produced
#' @return A DOT object ready to be fed to the algorithm
#' @export
#'
#' @examples
#' data(dot.sample)
#' dot.ref <- setup.ref(dot.sample$ref$counts[, 1:1000], dot.sample$ref$labels[1:1000], 2)
#' dot.srt <- setup.srt(dot.sample$srt$counts, dot.sample$srt$coordinates)
#' dot <- create.DOT(dot.srt, dot.ref)
create.DOT <- function(srt, ref, ls_solution = TRUE)
{
  # Now both ref_data and srt_data must have genes for columns
  cg <- intersect(colnames(ref$X[, which(colSums(ref$X) > 0)]), colnames(srt$X))

  if(length(cg) == 0)
  {
    stop("No common genes found between ref_data and srt_data.")
  }

  srt$X <- as.matrix(srt$X[, cg])
  ref$X <- ref$X[, cg]

  initial_Y <- matrix(0, 0, 0)
  if(ls_solution)
  {
    initial_Y <- ls_sol(ref$X, srt$X, lambda = 100)
  }

  return(Dot(ref = ref, srt = srt, solution = initial_Y, weights = matrix(0, 0, 0)))
}


#' A wrapper for running the DOT algorithm for high-resolution spatial data with suggested parameters
#'
#' @param object A DOT object created using create.DOT().
#' @param ratios_weight A value between 0 and 1 for matching ratio of cell types
#' @param iterations Integer. Maximum number of iterations of FW
#' @param verbose Boolean. Whether progress should be displayed.
#' @return A DOT object with the produced results contained in the weights slot
#' @export
#'
#' @examples
#' data(dot.sample)
#' dot.ref <- setup.ref(dot.sample$ref$counts[, 1:1000], dot.sample$ref$labels[1:1000], 2)
#' dot.srt <- setup.srt(dot.sample$srt$counts, dot.sample$srt$coordinates)
#' dot <- create.DOT(dot.srt, dot.ref)
#' # No. iterations is reduced to 10 for this example (default is 100)
#' dot <- run.DOT.highresolution(dot, iterations = 10)
run.DOT.highresolution <- function(object, ratios_weight = 0, iterations = 100, verbose = FALSE)
{
  return(.run.DOT(object, ratios_weight = ratios_weight,
                  sparsity_coef = 0.6, max_size = 1, verbose = verbose, iterations = iterations))
}


#' A wrapper for running the DOT algorithm for low-resolution spatial data with suggested parameters
#'
#' @param object A DOT object created using create.DOT().
#' @param ratios_weight A value between 0 and 1 for matching ratio of cell types
#' @param max_spot_size An upper bound on the size of spots. Default is 20, can be set to a higher value for lower resolution (eg, 200 for ST).
#' @param iterations Integer. Maximum number of iterations of FW
#' @param verbose Boolean. Whether progress should be displayed.
#' @return A DOT object with the produced results contained in the weights slot
#' @export
#'
#' @examples
#' data(dot.sample)
#' dot.ref <- setup.ref(dot.sample$ref$counts[, 1:1000], dot.sample$ref$labels[1:1000], 2)
#' dot.srt <- setup.srt(dot.sample$srt$counts, dot.sample$srt$coordinates)
#' dot <- create.DOT(dot.srt, dot.ref)
#' # No. iterations is reduced to 10 for this example (default is 100)
#' dot <- run.DOT.lowresolution(dot, iterations = 10)
run.DOT.lowresolution <- function(object, ratios_weight = 0,
                                  max_spot_size = 20, iterations = 100, verbose = FALSE)
{
  return(.run.DOT(object, ratios_weight = ratios_weight,
                  max_size = max_spot_size, sparsity_coef = 0.4,
                  verbose = verbose, iterations = iterations))
}

#' The internal DOT algorithm
#'
#' @param object A DOT object created using create.DOT().
#' @param ratios_weight Penalty weight for matching abundances
#' @param sparsity_coef A value between 0 (mixed) and 1 (sparse)
#' @param max_size An upper bound on the size of spots.
#' @param min_size A lower bound on the size of spots.
#' @param gap_threshold Upper bound on relative optimality gap
#' @param iterations Integer. Maximum number of iterations of FW
#' @param verbose Boolean. Whether progress should be displayed.
#' @return A DOT object with the produced results contained in the weights slot
#'
#' @keywords internal
#' @noRd
#'
.run.DOT <- function(object, ratios_weight = 1, sparsity_coef = 1, max_size = 20, min_size = 1,
                     iterations = 100, gap_threshold = 0.01, verbose = TRUE, ...)
{
  ST_X <- object@srt$X # S * G

  if(!is.null(object@srt$b))
  {
    ST_X <- fast_sweep(ST_X, MARGIN = 2, object@srt$b, "*")
  }

  ST_sizes <- object@srt$S # S

  SC_X <- object@ref$X # C * G
  SC_clusters <- object@ref$C # K
  if(is.null(SC_clusters))
    SC_clusters <- stats::setNames(as.list(1:nrow(object@ref$X)), rownames(object@ref$X))
  SC_ratios <- object@ref$R # K

  G <- ncol(ST_X)
  if(G != ncol(SC_X))
    stop("Invalid arguments")

  S <- nrow(ST_X)
  C <- nrow(SC_X)
  K <- length(SC_clusters)
  Km <- rep(0, C)
  for (k in 1:K) {
    Km[SC_clusters[[k]]] <- k
  }

  min_iterations <- min(10, iterations)
  marginal_improvement <- gap_threshold
  neighborhood_obj <- "JS"
  abundance_obj <- "JS"
  lambda_e <- 0
  fit.genes <- FALSE
  sqrtd <- TRUE
  diminishing_entropy <- FALSE

  epsilon <- 0
  epsilon_threshold <- 1e-2
  epsilon_step <- 0.1
  if(diminishing_entropy)
    epsilon <- 1e-1

  r_ST <- rep(0.9*min_size + 0.1*max_size, S)
  if(!is.null(ST_sizes))
  {
    r_ST <- ST_sizes
  }

  n_ST <- sum(r_ST)

  if(is.null(SC_ratios))
  {
    SC_ratios <- rep(1/K, K)
    names(SC_ratios) <- names(SC_clusters)
    ratios_weight <- 0 #cannot penalize abundance if ratios are not available
  }

  SC_ratios <- SC_ratios/sum(SC_ratios)

  r_SC <- SC_ratios * n_ST #expected cell types
  r_SC_ex <- r_SC[Km]
  w_SC <- rep(0, C)
  for(k in 1:K)
  {
    w_SC[SC_clusters[[k]]] <- r_SC[k]/length(SC_clusters[[k]])
  }
  w_SC <- w_SC / sum(w_SC) * C

  arguments <- list(...)
  inner_params <- arguments[["inner_params"]]

  if(is.null(inner_params))
  {
    inner_params <- c(1, 0.25, 0, 0.01)

    if(max_size == 1)
      inner_params[2] <- 1
  }else if(length(inner_params) < 4)
  {
    inner_params <- c(inner_params, rep(0, 4-length(inner_params)))
  }

  l_A <- ratios_weight / max_size
  l_E <- lambda_e * n_ST
  l_G <- inner_params[1] * S / G
  l_I <- inner_params[2]
  l_sp <- l_I * sparsity_coef / max_size
  l_C <- inner_params[3] * S / C

  if(is.null(object@srt$P))
  {
    l_S <- 0
  }else
  {
    l_S <- inner_params[4] * S / (max_size* nrow(object@srt$P))

    pairs_i <- object@srt$P$i
    pairs_j <- object@srt$P$j
    pairs_w <- object@srt$P$w
  }

  if(l_G > 0)
  {
    ST_Xg <- t(normalize(t(ST_X)))
  }

  SC_Xn <- normalize(SC_X)
  ST_Xn <- normalize(ST_X)

  env_min <- 1e-2
  env_slope <- 0.25/env_min
  env_threshold <- 4*env_min^2

  sqrt_env <- Vectorize(function(v)
  {
    if(v < env_threshold)
    {
      return(env_slope*v+env_min)
    }else
    {
      return(sqrt(v))
    }
  })

  sqrt_env_grad <- Vectorize(function(v)
  {
    if(v < env_threshold)
    {
      return(env_slope)
    }else
    {
      return(0.5/sqrt(v))
    }
  })

  bt <- object@srt$b
  if(is.null(bt))
    bt <- rep(1, G)
  gene_scale_sum <- sum(bt)
  gene_scale_eps <- 1
  gene_scale_min <- 0.1*(gene_scale_sum / G)

  if(fit.genes)
  {
    l_C <- 0
    sparsity_coef <- 0
    l_sp <- 0

    if(l_I == 0)
      l_I <- 1
  }

  Yt <- NULL
  if(nrow(object@solution) > 0)
  {
    Yt <- t(object@solution)
    if(nrow(Yt) != C | ncol(Yt) != S)
    {
      Yt <- NULL
    }else
    {
      Yt[which(Yt < 0)] <- 0
      cs <- colSums(Yt)
      Yt[, which(cs < 1e-3)] <- 1/C

      cs_factors <- rep(1, S)
      if(sparsity_coef > 0.5)
      {
        cs_high <- which(cs >= 1e-3)
        cs_factors[cs_high] <- 1 / cs[cs_high]
      }else
      {
        cs_high <- which(cs > max_size)
        cs_factors[cs_high] <- max_size / cs[cs_high]
        cs_low <- which(cs < min_size & cs >= 1e-3)
        cs_factors[cs_low] <- min_size / cs[cs_low]
      }

      cs_weight <- 0.99
      Yt <- Yt * matrix(cs_factors, nrow = C, ncol = S, byrow = TRUE)
      if(any(Yt == 0))
        Yt <- Yt*cs_weight + matrix((1-cs_weight)/C, C, S)
    }
  }

  linear_dcosine <- NULL
  if(l_sp > 0 | is.null(Yt))
  {
    linear_dcosine <- 1-(SC_Xn %*% t(ST_Xn))

    if(sqrtd)
    {
      linear_dcosine[which(linear_dcosine < 0)] <- 0
      linear_dcosine <- sqrt(linear_dcosine)
    }
  }

  iteration_start <- Sys.time()

  if(is.null(Yt))
  {
    initial_ratios <- rep(0, C)
    for(k in 1:K)
      initial_ratios[SC_clusters[[k]]] <- SC_ratios[k]/length(SC_clusters[[k]])
    Yt <- outer(initial_ratios, r_ST)

    mix_weight <- 0.1
    Yt <- Yt*mix_weight
    for (i in 1:S)
    {
      c_min <- which.min(linear_dcosine[, i])
      Yt[c_min, i] <- Yt[c_min, i] + (1-mix_weight)*r_ST[i]
    }
  }

  f <- Inf #best value (upper bound)
  f_en <- Inf
  LB <- -Inf #best lower bound
  Y <- NULL #best solution
  b <- NULL

  iteration_cols <- c("ft", "UB", "LB", "Gap", "Time", "d_ST", "d_SC", "d_Lin", "d_G", "d_S", "err_R", "err_E", "step_size")
  iteration_info <- matrix(NA, nrow = iterations, ncol = length(iteration_cols))
  colnames(iteration_info) <- iteration_cols

  lg2 <- log(2)

  iteration <- 1
  converged <- FALSE
  while(!converged)
  {
    Ytk <- rowsum(Yt, Km)
    rho_tk <- rowSums(Ytk)
    rho_tk[which(rho_tk < 1e-10)] <- 1e-10

    rho_t_ex <- rho_tk[Km]

    ratio_error <- 0
    if(l_A > 0)
    {
      if(abundance_obj == "JS")
      {
        rho_avg <- (rho_tk + r_SC)/2
        log_rho <- 0.5*safelog2(rho_tk/rho_avg)
        log_rsc <- 0.5*safelog2(r_SC/rho_avg)

        ratio_error <- sum(rho_tk*log_rho) + sum(r_SC*log_rsc)

        log_rho_ex <- log_rho[Km]
        d_ratio <- l_A*log_rho_ex - l_E/rho_t_ex
      }else
      {
        ratio_error <- sum((rho_tk-r_SC)^2)
        d_ratio <- l_A*2*sqrt_env_grad(ratio_error)*(rho_t_ex-r_SC_ex) - l_E/rho_t_ex
        ratio_error <- sqrt_env(ratio_error)
      }

      Dt <- matrix(d_ratio, nrow = C, ncol = S, byrow = F)
    }else
    {
      Dt <- matrix(0, nrow = C, ncol = S)
    }

    Dg <- rep(0, G)

    dcosine_ST <- 0
    dcosine_G <- 0
    if((sparsity_coef < 1 & l_I > 0) | l_G > 0 | fit.genes)
    {
      ST_Xt <- t(Yt) %*% SC_X

      ST_De <- NULL

      if((sparsity_coef < 1 & l_I > 0 ) | fit.genes)
      {
        ST_Xt_inorms <- matrix_norm(ST_Xt, 1) #apply(ST_Xt, 1, l2)

        # ST_Xn ==> ST_X normalized by i

        if(fit.genes)
        {
          ST_Xn <- normalize(fast_sweep(ST_X, MARGIN=2, bt, "*"))
        }

        ST_Xt_n <- fast_sweep(ST_Xt, 1, ST_Xt_inorms, "/")
        ST_XX <- ST_Xn * ST_Xt_n
        csi <- rowSums(ST_XX)
        di <- 1 - csi

        d_i_grad <- rep(1, S)
        if(sqrtd)
        {
          d_i_grad <- sqrt_env_grad(di)
          di <- sqrt_env(di)
        }

        dcosine_ST <- sum(di)

        if(sparsity_coef < 1 & l_I > 0)
        {
          ST_De <- l_I*(1-sparsity_coef)*fast_sweep(ST_Xn - fast_sweep(ST_Xt_n, 1, csi, "*"), 1, d_i_grad/ST_Xt_inorms, "*")
        }

        if(fit.genes)
        {
          Dg <- -as.vector((d_i_grad %*% ST_XX - (csi * d_i_grad) %*% ST_Xn^2)) / bt
        }
      }

      if(l_G > 0)
      {
        # start <- Sys.time()
        # ST_Xg ==> ST_X normalized by g
        ST_Xt_gnorms <- matrix_norm(ST_Xt, 2) #apply(ST_Xt, 2, l2)
        ST_Xt_gn <- fast_sweep(ST_Xt, 2, ST_Xt_gnorms, "/")
        csg <- colSums(ST_Xt_gn * ST_Xg)
        dg <- 1 - csg

        if(sqrtd)
        {
          dg_coefs <- sqrt_env_grad(dg)/ST_Xt_gnorms
          dg <- sqrt_env(dg)
        }else
        {
          dg_coefs <- 1/ST_Xt_gnorms
        }

        if(is.null(ST_De))
        {
          ST_De <- l_G*fast_sweep(ST_Xg - fast_sweep(ST_Xt_gn, 2, csg, "*"), 2, dg_coefs, "*")
        }else
        {
          ST_De <- ST_De + l_G*fast_sweep(ST_Xg - fast_sweep(ST_Xt_gn, 2, csg, "*"), 2, dg_coefs, "*")
        }

        dcosine_G <- sum(dg)
      }

      # negative comes from 1 - cosine
      Dt <- Dt - t(ST_De %*% t(SC_X))
    }

    dcosine_SC <- 0
    if(l_C > 0)
    {
      SC_Xt <- Yt %*% ST_X

      SC_Xt_inorms <- matrix_norm(SC_Xt, 1) #apply(SC_Xt, 1, l2)
      SC_Xt_n <- fast_sweep(SC_Xt, 1, SC_Xt_inorms, "/")
      csc <- rowSums(SC_Xt_n * SC_Xn)
      dc <- 1 - csc
      dc_coefs <- w_SC/SC_Xt_inorms

      if(sqrtd)
      {
        dc_coefs <- dc_coefs*sqrt_env_grad(dc)
        dc <- sqrt_env(dc)
      }

      SC_De <- t(fast_sweep(SC_Xn - fast_sweep(SC_Xt_n, 1, csc, "*"), 1, dc_coefs, "*"))
      dcosine_SC <- sum(dc*w_SC)

      Dt <- Dt - l_C*t(ST_X %*% SC_De)
    }

    dcosine_Lin <- 0
    if(l_sp > 0)
    {
      Dt <- Dt + l_sp*linear_dcosine
      dcosine_Lin <- sum(Yt * linear_dcosine)
    }

    d_S <- 0
    if(l_S > 0)
    {
      Dtk <- matrix(0, nrow = K, ncol = S)

      if(neighborhood_obj == "JS")
      {
        for (p in 1:nrow(object@srt$P)) {
          i <- pairs_i[p]
          j <- pairs_j[p]
          w <- pairs_w[p]*0.5/lg2 # 0.5 is because of definition of JS; lg2=log(2) is because JS is in base 2
          Ym <- 0.5*(Ytk[, i] + Ytk[, j])

          for(ii in c(i,j))
          {
            L_ii <- safelog2(Ytk[,ii]/Ym)
            d_S <- d_S + w*sum(Ytk[,ii]*L_ii)
            Dtk[, ii] <- Dtk[, ii] + l_S*w*L_ii
          }
        }
      }else
      {
        for (p in 1:nrow(object@srt$P)) {
          i <- pairs_i[p]
          j <- pairs_j[p]
          w <- pairs_w[p]
          Yd <- (Ytk[, i] - Ytk[, j])
          d_S <- d_S + w*sum(Yd^2)

          Dtk[, i] <- Dtk[, i] + (2*l_S*w)*Yd
          Dtk[, j] <- Dtk[, j] - (2*l_S*w)*Yd
        }
      }

      Dt <- Dt + Dtk[Km, ]
    }

    ratio_equity_error <- 0
    if(l_E > 0)
      ratio_equity_error <- -sum(log(rho_tk))

    Yt_h <- matrix(0, nrow = C, ncol = S)
    bt_h <- rep(0, G)

    if(iteration == 1 & epsilon > 0)
    {
      epsilon_scale <- mean(abs(Dt))
      epsilon <- epsilon_scale * epsilon
      epsilon_threshold <- epsilon_scale * epsilon_threshold
    }

    eps_failed <- TRUE
    if(epsilon > epsilon_threshold)
    {
      eps_failed <- FALSE
      for (i in 1:S)
      {
        sm <- exp(-Dt[, i]/epsilon)
        sms <- sum(sm)

        if(!is.finite(sms))
        {
          if(i > 1)
            Yt_h[, 1:i] <- 0
          eps_failed <- TRUE

          break()
        }else
        {
          if(sms > max_size)
          {
            sm <- sm * (max_size/sms)
          }else if(sms < min_size)
          {
            sm <- sm * (min_size/sms)
          }
        }
        sm[which(sm < 1e-10)] <- 0
        Yt_h[, i] <- sm
      }
    }

    if(eps_failed)
    {
      epsilon <- 0

      for (i in 1:S)
      {
        kk <- which.min(Dt[, i])
        Yt_h[kk, i] <- ifelse(Dt[kk,i] >= 0, min_size, max_size)
      }
    }

    Y_diff <- Yt - Yt_h

    ft <- l_I*(1-sparsity_coef)*dcosine_ST + l_sp*dcosine_Lin + l_C*dcosine_SC + l_G*dcosine_G +
      l_S*d_S + l_A*ratio_error + l_E*ratio_equity_error

    #gap:
    gap <- sum(Dt * Y_diff)

    entropy_t <- 0
    if(epsilon > 0)
    {
      entropy_t <- epsilon*sum(Entropy(Yt))
      entropy_h <- epsilon*sum(Entropy(Yt_h))
      gap <- gap + (entropy_t - entropy_h)
    }

    if(fit.genes)
    {
      failed <- TRUE
      # gene_scale_eps <- 0
      if(gene_scale_eps > 0)
      {
        # if(iteration == 1)
        # gs_eps <- max(1, mean(abs(Dg)))*gene_scale_eps
        gs_eps <- max(abs(Dg))/10
        bt_h <- exp(-Dg/gs_eps)
        bt_h <- bt_h/sum(bt_h) * gene_scale_sum

        if(all(is.finite(bt_h)))
        {
          entropy_bt <- gs_eps*sum(Entropy(bt))
          entropy_bh <- gs_eps*sum(Entropy(bt_h))

          gap <- gap + (entropy_bt - entropy_bh)

          failed <- FALSE
        }else
        {
          failed <- TRUE
        }
      }

      if(failed)
      {
        bt_h <- rep(gene_scale_min, G)
        bt_h[which.min(Dg)] <- gene_scale_sum - (G-1)*gene_scale_min
      }

      b_diff <- bt - bt_h

      gap <- gap + sum(Dg * b_diff)
    }

    if(ft < f)
    {
      f <- ft
      Y <- Yt
      b <- bt
    }

    if(ft + entropy_t < f_en)
      f_en <- ft + entropy_t

    # if(LB < ft + entropy_t - gap)
    LB <- ft + entropy_t - gap

    # gap <- f_en - LB

    if(abs(f_en) > 1e-10)
      gap <- gap/abs(f_en)

    # step_i <- ifelse(epsilon > 0, 6*(iteration + 1)/((iteration + 2)*(2*iteration + 3)), 2/(iteration + 1))
    step_i <- 2/(iteration + 1)
    step <- min(0.99, step_i)

    time <- as.numeric(Sys.time() - iteration_start, units = "secs")
    if(verbose)
    {
      if(iteration %% 10 == 1)
        cat("#: obj, UB, LB, d-ST, d-SC, d-Lin, d-G, d-N, d-R, d-E, d-H, gap, alpha, time\n")
      cat(sprintf("%d: %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n",
                  iteration, ft, f, LB, dcosine_ST, dcosine_SC, dcosine_Lin, dcosine_G, d_S, ratio_error, ratio_equity_error,
                  entropy_t, gap, step, time))
    }

    # c("ft", "UB", "LB", "Gap", "Time", "d_ST", "d_SC", "d_Lin", "d_S", "err_R", "err_E", "step_size")
    iteration_info[iteration, ] <- c(ft, f, LB, gap, time, dcosine_ST, dcosine_SC, dcosine_Lin, dcosine_G, d_S, ratio_error, ratio_equity_error, step)

    iteration_start <- Sys.time()

    if(gap <= gap_threshold)
    {
      converged <- TRUE
    }

    if(converged & diminishing_entropy & epsilon > 0)
    {
      converged <- FALSE
      LB <- -Inf

      epsilon <- epsilon * epsilon_step

      if(epsilon < epsilon_threshold)
      {
        epsilon <- 0
        diminishing_entropy <- FALSE
      }

      message(sprintf("Resetting epsilon to %g ...", epsilon))
    }

    if(converged & min_iterations > 0 & iteration < min_iterations)
    {
      converged <- FALSE
      message("Resetting lower-bound to avoid local optima...")
      LB <- -Inf
    }

    if(converged & marginal_improvement > 0)
    {
      if((iteration_info[max(2, iteration - 5), 2]-f)/f >= marginal_improvement)
      {
        converged <- FALSE
        message("Resetting lower-bound as still making improvement...")
        LB <- -Inf
      }
    }

    if(step <= 1e-5)
      converged <- TRUE

    if(iterations > 0 & iteration >= iterations)
      converged <- TRUE

    if(converged)
      break

    Yt <- Yt - step * Y_diff
    if(fit.genes)
      bt <- bt - step * b_diff

    iteration <- iteration + 1
  }

  rownames(Y) <- rownames(SC_X)
  colnames(Y) <- rownames(ST_X)

  if(K > 0)
  {
    weights <- c()
    for(ct in names(SC_clusters))
      weights <- cbind(weights, colSums(Y[SC_clusters[[ct]], ,drop = FALSE]))
    colnames(weights) <- names(SC_clusters)
    rownames(weights) <- rownames(ST_X)
  }else
  {
    weights <- t(Y)
  }

  object@weights <- weights
  object@solution <- t(Y)
  object@srt$b <- stats::setNames(b, colnames(object@srt$X))
  object@history <- as.data.frame(iteration_info[1:iteration, ])

  return(object)
}
