#' @keywords internal
l2 <- function(x) norm(x, type="2")

#' @keywords internal
fast_centroid <- function(X, ann, MARGIN = 1)
{
  if(is(X, "Matrix"))
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

  # if(!(FUN %in% c("*", "+", "/", "-")))
  #   stop("Operator not supported")
  #
  # xs <- NULL
  # if(FUN == "-")
  # {
  #   xs <- rcpp_sweep(x, MARGIN, -STATS, "+")
  # }else if(FUN == "/")
  # {
  #   xs <- rcpp_sweep(x, MARGIN, 1/STATS, "*")
  # }else
  # {
  #   xs <- rcpp_sweep(x, MARGIN, STATS, FUN)
  # }
  #
  # if(!is.null(xs))
  # {
  #   dimnames(xs) <- dimnames(x)
  # }
  #
  # return(xs)
}

#' @keywords internal
normalize <- function(x, MARGIN = 1){
  x <- fast_sweep(x, MARGIN = MARGIN, matrix_norm(x, MARGIN), "/")

  x[which(!is.finite(x))] <- 0

  return(x)
}

#' @keywords internal
e10 <- exp(-10)

#' @keywords internal
Entropy <- Vectorize(function(y){
  if(y < e10)
  {
    return(-e10*10)
  }else
  {
    return(y*log(y))
  }
})

#' @keywords internal
softmin <- function(x, epsilon = 1e-2){
  ex <- exp((max(x)-x)/epsilon)
  return(ex/sum(ex))
}

#' @keywords internal
softmin_safe <- function(x, epsilon){
  s <- softmin(x, epsilon)
  if(any(is.nan(s)))
  {
    ss <- which.max(x)
    s[ss] <- 0
    s[-ss] <- softmin_safe(x[-ss], epsilon)
  }
  return(s)
}

#' @keywords internal
smallL = -20

#' @keywords internal
safelog2 <- function(x)
{
  lg <- log2(x)
  lg[which(is.nan(lg))] <- 0
  lg[which(is.infinite(lg))] <- smallL
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

# # @keywords internal
# rcpp_sweep <- Rcpp::cppFunction('
# NumericMatrix rcpp_sweep(NumericMatrix x, int MARGIN, NumericVector STATS, String FUN) {
#   int n = x.nrow();
#   int m = x.ncol();
#   NumericMatrix out(n, m);
#
#   if(MARGIN == 1) {
#     if(FUN == "+") {
#       for(int i = 0; i < n; ++i)
#         for (int j = 0; j < m; ++j)
#           out(i,j) = x(i,j) + STATS[i];
#     }
#     else if(FUN == "*") {
#       for(int i = 0; i < n; ++i)
#         for (int j = 0; j < m; ++j)
#           out(i,j) = x(i,j) * STATS[i];
#     }
#   }
#   else if(MARGIN == 2) {
#     if(FUN == "+") {
#       for(int i = 0; i < n; ++i)
#         for (int j = 0; j < m; ++j)
#           out(i,j) = x(i,j) + STATS[j];
#     }
#     else if(FUN == "*") {
#       for(int i = 0; i < n; ++i)
#         for (int j = 0; j < m; ++j)
#           out(i,j) = x(i,j) * STATS[j];
#     }
#   }
#
#   return out;
# }
# ')

# matrix_norm <- Rcpp::cppFunction('
# NumericVector matrix_norm(NumericMatrix x, int MARGIN) {
#
#   int n = x.nrow();
#   int m = x.ncol();
#
#   if(MARGIN == 1){
#     NumericVector nrm(n);
#     for(int i = 0; i < n; ++i){
#       double nr = 0;
#       for (int j = 0; j < m; ++j)
#         nr += x(i,j)*x(i,j);
#       nrm[i] = sqrt(nr);
#     }
#
#     return(nrm);
#   }
#   else{
#     NumericVector nrm(m);
#
#     for(int j = 0; j < m; ++j){
#       double nr = 0;
#       for (int i = 0; i < n; ++i)
#         nr += x(i,j)*x(i,j);
#       nrm[j] = sqrt(nr);
#     }
#     return(nrm);
#   }
# }
# ')
