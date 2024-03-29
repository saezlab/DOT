% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dot.R
\name{setup.srt}
\alias{setup.srt}
\title{Processing the target spatial data}
\usage{
setup.srt(
  srt_data,
  srt_coords = NULL,
  th.spatial = 0.84,
  th.nonspatial = 0,
  th.gene.low = 0.01,
  th.gene.high = 0.99,
  remove_mt = TRUE,
  radius = "auto",
  verbose = FALSE
)
}
\arguments{
\item{srt_data}{A gene x spot/cell matrix of gene expressions. Can be a matrix-like object or a Seurat/AnnData object}

\item{srt_coords}{A matrix-like object with two columns (x and y coordinates of spots). If set to NULL and srt_data is Seurat/AnnData object, coords will be extracted from srt_data}

\item{th.spatial}{A value between 0 and 1. Threshold on similarity of adjacent spots}

\item{th.nonspatial}{A value between 0 and 1. Threshold on similarity of non-adjacent spots}

\item{th.gene.low}{Minimum percentage of spots that a valid gene must be expressed in.}

\item{th.gene.high}{Maximum percentage of spots that a valid gene must be expressed in.}

\item{remove_mt}{Boolean. Whether mitochondrial genes must be removed.}

\item{radius}{Adjacency radius. If set to 'auto' it is computed using the coordinates of spots}

\item{verbose}{Boolean. Whether progress should be displayed.}
}
\value{
A list containing the processed srt data.
}
\description{
Processing the target spatial data
}
\examples{
data(dot.sample)
dot.srt <- setup.srt(dot.sample$srt$counts, dot.sample$srt$coordinates)
}
