# Processing the reference single-cell data

Processing the reference single-cell data

## Usage

``` r
setup.ref(
  ref_data,
  ref_annotations = NULL,
  ref_subcluster_size = 10,
  max_genes = 5000,
  remove_mt = TRUE,
  verbose = FALSE
)
```

## Arguments

- ref_data:

  A gene x cell matrix of gene expressions. Can be a matrix-like object
  or a Seurat/AnnData object.

- ref_annotations:

  A character vector (one for each cell) or a single vector pointing to
  the slot in the Seurat/AnnData object

- ref_subcluster_size:

  An integer. Maximum number of sub-clusters per sub-population.

- max_genes:

  An integer. Maximum number of genes to pick.

- remove_mt:

  Boolean. Whether mitochondrial genes must be removed.

- verbose:

  Boolean. Whether progress should be displayed.

## Value

A list containing the processed ref data.

## Examples

``` r
data(dot.sample)
dot.ref <- setup.ref(dot.sample$ref$counts[, 1:1000], dot.sample$ref$labels[1:1000], 2)
```
