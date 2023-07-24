#' A utility function for computing l2 of a vector
#' @param x A vector
#' @return A numeric value
#' @keywords internal
#' @noRd
#'
l2 <- function(x) norm(x, type="2")

#' A utility function for computing row/column-wise means of a matrix
#' @param X A numeric matrix
#' @param ann A character vector
#' @param MARGIN An integer (1 or 2)
#' @return A matrix of centroids
#' @keywords internal
#' @noRd
#'
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

#' A utility function for computing absolute/relateive frequency of unique items in a vector
#' @param vals A vector
#' @param normalize Whether frequencies are relative
#' @return A vector
#' @keywords internal
#' @noRd
#'
int.table <- function(vals, normalize = FALSE)
{
  counts <- table(vals)
  counts <- stats::setNames(as.integer(counts), names(counts))
  if(normalize)
    counts <- counts / sum(counts)

  return(counts)
}

#' A utility function for computing row/column-wise l2 norms of a matrix
#' @param X A numeric matrix
#' @param MARGIN An integer (1 or 2)
#' @return A vector
#' @keywords internal
#' @noRd
#'
matrix_norm <- function(x, MARGIN)
{
  return(apply(x, MARGIN, l2))
}

#' A wrapper for sweep
#' @param x A numeric matrix
#' @param MARGIN An integer (1 or 2)
#' @param STATS A vector
#' @param FUN A character representation of operator
#' @return A matrix
#' @keywords internal
#' @noRd
#'
fast_sweep <- function(x, MARGIN, STATS, FUN)
{
  return(sweep(x, MARGIN, STATS, FUN))
}

#' A utility function for row/column-wise normalization of a matrix
#' @param x A numeric matrix
#' @param MARGIN An integer (1 or 2)
#' @return A matrix
#' @keywords internal
#' @noRd
#'
normalize <- function(x, MARGIN = 1){
  x <- fast_sweep(x, MARGIN = MARGIN, matrix_norm(x, MARGIN), "/")
  x[which(!is.finite(x))] <- 0

  return(x)
}

#' @keywords internal
E10 <- exp(-10)

#' A utility function for computing entropy of a vector
#' @param y A numeric vector
#' @return A vector
#' @keywords internal
#' @noRd
#'
Entropy <- Vectorize(function(y){
  if(y < E10)
  {
    return(-E10*10)
  }else
  {
    return(y*log(y))
  }
})

#' A utility function for safe log
#' @param x A numeric vector
#' @return A vector
#' @keywords internal
#' @noRd
#'
safelog2 <- function(x)
{
  lg <- log2(x)
  lg[which(is.nan(lg))] <- 0
  lg[which(is.infinite(lg))] <- -20
  return(lg)
}
