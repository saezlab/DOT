# Processing the target spatial data

Processing the target spatial data

## Usage

``` r
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
```

## Arguments

- srt_data:

  A gene x spot/cell matrix of gene expressions. Can be a matrix-like
  object or a Seurat/AnnData object

- srt_coords:

  A matrix-like object with two columns (x and y coordinates of spots).
  If set to NULL and srt_data is Seurat/AnnData object, coords will be
  extracted from srt_data

- th.spatial:

  A value between 0 and 1. Threshold on similarity of adjacent spots

- th.nonspatial:

  A value between 0 and 1. Threshold on similarity of non-adjacent spots

- th.gene.low:

  Minimum percentage of spots that a valid gene must be expressed in.

- th.gene.high:

  Maximum percentage of spots that a valid gene must be expressed in.

- remove_mt:

  Boolean. Whether mitochondrial genes must be removed.

- radius:

  Adjacency radius. If set to 'auto' it is computed using the
  coordinates of spots

- verbose:

  Boolean. Whether progress should be displayed.

## Value

A list containing the processed srt data.

## Examples

``` r
data(dot.sample)
dot.srt <- setup.srt(dot.sample$srt$counts, dot.sample$srt$coordinates)
#> Computing spatial radius
```
