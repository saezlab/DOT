#' @keywords internal
l2 <- function(x) norm(x, type="2")

#' @keywords internal
fast_centroid <- function(X, ann, MARGIN = 1)
{
  if(methods::is(X, "Matrix"))
  {
    group <- as.factor(ann)
    group_mat <- Matrix::sparse.model.matrix(~ 0 + group, transpose = TRUE)
    rownames(group_mat) <- stringr::str_extract(rownames(group_mat), "(?<=^group).+")
    sizes <- rowSums(group_mat)

    if(MARGIN == 1)
    {
      return(sweep(group_mat %*% X, MARGIN, sizes, "/"))
    }else
    {
      return(sweep(X %*% t(group_mat), MARGIN, sizes, "/"))
    }
  }else
  {
    sizes <- int.table(ann, FALSE)

    if(MARGIN == 1)
    {
      centroids <- rowsum(as.matrix(X), ann)
      return(fast_sweep(centroids, 1, sizes[rownames(centroids)], "/"))
    }else
    {
      centroids <- t(rowsum(t(as.matrix(X)), ann))
      return(fast_sweep(centroids, 2, sizes[colnames(centroids)], "/"))
    }
  }
}

#' @keywords internal
int.table <- function(vals, normalize = FALSE)
{
  counts <- table(vals)
  counts <- stats::setNames(as.integer(counts), names(counts))
  if(normalize)
    counts <- counts / sum(counts)

  return(counts)
}

#' @keywords internal
matrix_norm <- function(x, MARGIN)
{
  return(apply(x, MARGIN, l2))
}

#' @keywords internal
fast_sweep <- function(x, MARGIN, STATS, FUN)
{
  return(sweep(x, MARGIN, STATS, FUN))
}

#' @keywords internal
normalize <- function(x, MARGIN = 1){
  x <- fast_sweep(x, MARGIN = MARGIN, matrix_norm(x, MARGIN), "/")
  x[which(!is.finite(x))] <- 0

  return(x)
}

#' @keywords internal
E10 <- exp(-10)

#' @keywords internal
Entropy <- Vectorize(function(y){
  if(y < E10)
  {
    return(-E10*10)
  }else
  {
    return(y*log(y))
  }
})

#' @keywords internal
safelog2 <- function(x)
{
  lg <- log2(x)
  lg[which(is.nan(lg))] <- 0
  lg[which(is.infinite(lg))] <- -20
  return(lg)
}

#' @keywords internal
ls_sol <- function(ref_X, srt_X, lambda = 10)
{
  C <- nrow(ref_X)
  Y <- solve((ref_X %*% t(ref_X)) + diag(lambda, C, C), ref_X %*% t(srt_X))
  Y[which(Y < 0)] <- 0
  rownames(Y) <- rownames(ref_X)
  colnames(Y) <- rownames(srt_X)

  return(t(Y))
}
