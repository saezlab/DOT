#' An internal method for aggregating ref data, sub-clustering and computing the centroids
#' @param ref_data A cell x gene matrix of gene expressions.
#' @param ref_annotations A character vector (one for each cell)
#' @param cluster_size An integer. Maximum number of sub-clusters per sub-population
#' @param th.inner_logfold log-fold threshold for selecting genes for sub-clustering
#' @param verbose Boolean. Whether progress should be displayed.
#' @return A list containing the processed ref data.
#' @keywords internal
#' @noRd
#'
aggregate_ref <- function(ref_data, ref_annotations, cluster_size, th.inner_logfold = 0.75, verbose = FALSE)
{
  major_types <- sort(as.character(unique(ref_annotations)))

  candidate_ind <- c()
  candidate_type <- c()
  for (u in major_types) {
    u_ind <- which(ref_annotations == u)
    if(length(u_ind) > 1000)
      u_ind <- sample(u_ind, 1000, replace = FALSE)
    candidate_ind <- c(candidate_ind, u_ind)
    candidate_type <- c(candidate_type, rep(u, length(u_ind)))
  }

  major_centroids <- fast_centroid(ref_data[candidate_ind, ], candidate_type)
  major_ratios <- int.table(as.character(ref_annotations), TRUE)

  if(cluster_size > 1)
  {
    if(verbose)
      cat("Clustering cell types\n")
    sub_annotations <- rep(NA, nrow(ref_data))

    for(u in major_types)
    {
      u_ind <- which(ref_annotations == u)
      if(length(u_ind) <= 1)
        next

      if(length(u_ind) > 10000)
        u_ind <- sample(u_ind, 10000, replace = FALSE)

      x <- ref_data[u_ind, ]

      K <- round(min(cluster_size, 2*log(nrow(x))-7))

      if(K <= 1)
      {
        sub_annotations[u_ind] <- 1
      }else
      {
        kmeans_genes <- colnames(x)

        if(length(kmeans_genes) > 500)
        {
          this_ct <- major_centroids[u, ]
          others <- setdiff(rownames(major_centroids), u)
          other_ct <- t(major_ratios[others] / sum(major_ratios[others])) %*% major_centroids[others, ]

          logfc_ct <- log((this_ct + 1e-9) / (other_ct[1, ] + 1e-9))

          kmeans_genes <- names(which(logfc_ct > th.inner_logfold)) # Without abs
          if(length(kmeans_genes) > 500)
          {
            kmeans_genes <- names(sort(logfc_ct, decreasing = T)[1:500])
          }
        }

        x <- as.matrix(x[, kmeans_genes])

        if(verbose)
          cat(sprintf("Clustering %d %s cells into %d clusters ...\n", nrow(x), u, K))

        km <- stats::kmeans(x, K, nstart = 10)

        noise <- int.table(km$cluster, TRUE)

        km$cluster[which(km$cluster %in% as.integer(which(noise < 0.025)))] <- NA

        sub_annotations[u_ind] <- km$cluster
      }
    }

    nonna <- which(!is.na(sub_annotations))
    ref_data <- ref_data[nonna, ]
    sub_annotations <-  sub_annotations[nonna]
    ref_annotations <-  ref_annotations[nonna]

    sub_annotations <- sprintf("%s__SC%d", ref_annotations, sub_annotations)

    sub_centroids <- fast_centroid(ref_data, sub_annotations)
    sub_types <- rownames(sub_centroids)

    major_types <- unique(ref_annotations)
    major_ratios <- major_ratios[major_types]
    major_ratios <- major_ratios / sum(major_ratios)

    clusters <- lapply(major_types, function(ct) which(startsWith(sub_types, sprintf("%s__SC", ct))))
    names(clusters) <- major_types

    sub_ratios <- int.table(as.character(sub_annotations), TRUE)
  }else
  {
    clusters <- stats::setNames(as.list(1:length(major_types)),  major_types)
    sub_centroids <- major_centroids
    sub_ratios <- major_ratios
  }

  return(list(major_centroids =  major_centroids, major_ratios = major_ratios, clusters = clusters,
              sub_centroids = sub_centroids, sub_ratios = sub_ratios))
}

#' An internal method for fast de gene analysis
#' @param ref_centroids A cluster x gene matrix of gene expression centroids.
#' @param ref_ratios Relative abundance of sub-populations
#' @param max_genes Maximum number of genes to select
#' @param verbose Boolean. Whether progress should be displayed.
#' @return A list
#' @keywords internal
#' @noRd
#'
get_de_genes <- function(ref_centroids, ref_ratios, max_genes, verbose = FALSE)
{
  Y_greedy <- NULL
  de_genes <- colnames(ref_centroids)

  if(length(de_genes) > max_genes)
  {
    if(verbose)
      cat("Gene filtering...\n")

    C <- nrow(ref_centroids)
    G <- ncol(ref_centroids)

    ref_centroids <- ref_centroids + 1e-9

    gene_scors <- matrix(0, nrow = C, ncol = G, dimnames = dimnames(ref_centroids))

    for(i in 1:C)
    {
      this_ct <- matrix(ref_centroids[i, ], nrow = C-1, ncol = G, byrow = T)
      other_ct <- ref_centroids[-i, ]
      logfc_ct <- log(this_ct / other_ct)
      lo <- apply(logfc_ct, 1, rank)
      lo <- t(max(lo) - lo + 1)
      # logfc_ct[which(lo > genes_per_pair)] <- 0
      # llg <- logfc_ct[, which(apply(logfc_ct, 2, max) > 0)]
      gene_scors[i, ] <- apply(lo, 2, stats::median)
    }

    de_genes <- names(sort(apply(gene_scors, 2, min))[1:max_genes])
    if(verbose)
      cat("Extractd", length(de_genes), "de genes.\n")
  }

  return(list(de_genes = de_genes, Y = Y_greedy))
}

#' An internal method for computing the spatial neighborhood radius
#' @param coordinates Coordinates of spots
#' @param neighbors An internal parameter
#' @param sample_size An internal parameter
#' @param th.quantile An internal parameter
#' @return A numeric value
#' @keywords internal
#' @noRd
#'
find_radius <- function(coordinates, neighbors = 8, sample_size = 2000, th.quantile = 0.9)
{
  N <- nrow(coordinates)
  if(N > sample_size)
  {
    sample_spots <- sample(1:N, sample_size)
    dists <- fields::rdist(coordinates[sample_spots, ], coordinates)
  }else
  {
    dists <- as.matrix(stats::dist(coordinates))
  }
  min_dists <- matrix(NA, nrow = nrow(dists), ncol = neighbors)
  for(r in 1:nrow(dists))
  {
    min_dists[r, ] <- sort(dists[r, which(dists[r,]>0)])[1:neighbors]
  }

  mid_dists <- as.vector(min_dists)

  radius <- as.numeric(stats::quantile(min_dists, th.quantile)) # 99%

  return(radius)
}

#' An internal method for finding the set of spatial pairs
#' @param nrm_X Normalized spot x gene counts
#' @param coordinates Coordinates of spots
#' @param radius Spatial radius. If 'auto' it's computed automatically
#' @param neighbors An internal parameter
#' @param th.quantile An internal parameter
#' @return A list
#' @keywords internal
#' @noRd
#'
find_spatial_pairs <- function(nrm_X, coordinates, radius = "auto", neighbors = 8, th.quantile = 0.9)
{
  # nrm_X is assumed to be normalized so that l2-norm of each row is 1

  N <- nrow(nrm_X)
  ST_DS <- NULL
  message("Computing spatial radius")
  if(radius == "auto")
  {
    radius <- find_radius(coordinates, neighbors = neighbors,
                          sample_size = 1000, th.quantile = th.quantile)
    radius <- radius * 1.05
    # print(radius)
    gc()
  }

  srtp <- c()
  batch_size <- 500
  for(b in 1:ceiling(N/batch_size))
  {
    spot_range <- ((b-1)*batch_size+1):(min(b*batch_size,N))

    ST_DS <- fields::rdist(coordinates[spot_range, ], coordinates)
    adj_pairs <- which(ST_DS < radius, arr.ind = TRUE)
    if(nrow(adj_pairs) > 0)
    {
      adj_pairs[, 1] <- spot_range[adj_pairs[, 1]]
      adj_w <- rowSums(nrm_X[adj_pairs[, 1],] * nrm_X[adj_pairs[, 2],])
      srtp <- rbind(srtp, cbind(adj_pairs, adj_w))
    }
  }

  colnames(srtp) <- c("i", "j", "w")
  srtp <- as.data.frame(srtp)
  srtp <- srtp[which(srtp$i < srtp$j), ]

  return(list(pairs = srtp, radius = radius))
}

#' An internal method for finding the set of spatial/non-spatial pairs
#' @param srt List containing data related to srt
#' @param th.spatial Threshold on similarity of adjacent pairs
#' @param th.nonspatial Threshold on similarity of non-adjacent pairs
#' @param max_size Maximum number of pairs
#' @param radius Spatial radius. If 'auto' it's computed automatically
#' @return A list
#' @keywords internal
#' @noRd
#'
get_pairs <- function(srt, th.spatial, th.nonspatial, max_size, radius = 'auto')
{
  nrm_X <- normalize(srt$X)
  srtp_info <- find_spatial_pairs(nrm_X, srt$C, radius = radius, th.quantile = 0.9)

  srtp <- srtp_info$pairs
  srtp <- srtp[which(srtp$w >= th.spatial), ]

  if(nrow(srtp) > max_size)
  {
    srtp <- srtp[order(srtp$w, decreasing = TRUE)[1:max_size], ]
  }
  else if(th.nonspatial > 0)
  {
    N <- nrow(nrm_X)
    srta <- as.data.frame(cbind(rep(1:N, N), rep(1:N, each = N)))
    colnames(srta) <- c("i", "j")
    srta$w <- as.vector(nrm_X %*% t(nrm_X))
    srta <- srta[which(srta$w >= th.nonspatial), ]
    srta <- srta[which(srta$i < srta$j), ]

    rownames(srtp) <- sprintf("%d_%d", srtp$i, srtp$j)
    rownames(srta) <- sprintf("%d_%d", srta$i, srta$j)

    extra_pairs <- srta[setdiff(rownames(srta), intersect(rownames(srtp), rownames(srta))), ]

    if(nrow(extra_pairs) > max_size - nrow(srtp))
    {
      extra_pairs <- extra_pairs[order(extra_pairs$w, decreasing = TRUE)[1:max_size], ]
    }

    srtp <- rbind(srtp, extra_pairs)
  }

  if(nrow(srtp) == 0)
  {
    message("No pairs satisfied the thresholds.")
    srtp <- NULL
  }

  return(srtp)
}

#' An internal method filtering MT-like genes
#' @param genes Vector of gene symbols
#' @return A Boolean vector
#' @keywords internal
#' @noRd
#'
is_mt <- function(genes)
{
  return(startsWith(genes, "MT-") | startsWith(genes, "HLA-") | startsWith(genes, "RPL"))
}

#' A plotting wrapper for drawing gene/abundance maps on tissue
#' @param spatial Coordinates
#' @param weights Weights to draw (spot x feature)
#' @param normalize Whether weights should sum up to 1 for each spot
#' @param ncol Number of columns in the grid
#' @param trans Color transformation scale
#' @param size Size of points
#' @param shape Shape of points
#' @param flip_y Whether y axis should be negated
#' @param viridis_option Viridis color options
#' @param legend_title Title of color legend
#' @param background Background color of plot
#' @return A ggplot object
#' @export
#'
#' @examples
#' data(dot.sample)
#' draw_maps(dot.sample$srt$coordinates, t(as.matrix(dot.sample$srt$counts[1:8, ])), normalize = FALSE)
draw_maps <- function(spatial, weights, normalize = TRUE, ncol = 4, trans = 'identity',
                      size = 0.5, shape = 19, flip_y = TRUE, viridis_option = "magma",
                      legend_title = "abundance", background = 'gray90')
{
  if(normalize)
  {
    weights <- sweep(weights, 1, rowSums(weights), "/")
  }

  plot_dt <- c()
  for (s in colnames(weights)) {
    sdt <- as.data.frame(spatial)
    sdt$abundance <- weights[, s]
    sdt$celltype <- s
    sdt <- sdt[order(sdt$abundance), ]
    plot_dt <- rbind(plot_dt, sdt)
  }

  plot_dt <- stats::setNames(as.data.frame(plot_dt), c("x", "y", "abundance", "celltype"))
  if(flip_y)
    plot_dt$y <- -plot_dt$y
  p <- ggplot2::ggplot(plot_dt, ggplot2::aes_string(x = "x", y = "y", color = "abundance"))+
    ggplot2::geom_point(size =  size, shape = shape)+
    ggplot2::scale_color_viridis_c(legend_title, option = viridis_option, trans = trans)+
    ggplot2::xlab("") + ggplot2::ylab("")+
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
          strip.background = ggplot2::element_rect(fill = 'transparent'),
          panel.background = ggplot2::element_rect(fill = background)
    )+
    ggplot2::facet_wrap(~celltype, ncol = ncol)

  return(p)
}
