# A wrapper for running the DOT algorithm for high-resolution spatial data with suggested parameters

A wrapper for running the DOT algorithm for high-resolution spatial data
with suggested parameters

## Usage

``` r
run.DOT.highresolution(
  object,
  ratios_weight = 0,
  iterations = 100,
  verbose = FALSE
)
```

## Arguments

- object:

  A DOT object created using create.DOT().

- ratios_weight:

  A value between 0 and 1 for matching ratio of cell types

- iterations:

  Integer. Maximum number of iterations of FW

- verbose:

  Boolean. Whether progress should be displayed.

## Value

A DOT object with the produced results contained in the weights slot

## Examples

``` r
data(dot.sample)
dot.ref <- setup.ref(dot.sample$ref$counts[, 1:1000], dot.sample$ref$labels[1:1000], 2)
dot.srt <- setup.srt(dot.sample$srt$counts, dot.sample$srt$coordinates)
#> Computing spatial radius
dot <- create.DOT(dot.srt, dot.ref)
# No. iterations is reduced to 10 for this example (default is 100)
dot <- run.DOT.highresolution(dot, iterations = 10)
```
