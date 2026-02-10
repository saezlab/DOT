# Creating a DOT object based on the processed ref and srt data

Creating a DOT object based on the processed ref and srt data

## Usage

``` r
create.DOT(srt, ref, ls_solution = TRUE)
```

## Arguments

- srt:

  A list containing the processed srt data produced by setup.srt

- ref:

  A list containing the processed ref data produced by setup.ref

- ls_solution:

  Boolean. Whether an initial solution based on LS should be produced

## Value

A DOT object ready to be fed to the algorithm

## Examples

``` r
data(dot.sample)
dot.ref <- setup.ref(dot.sample$ref$counts[, 1:1000], dot.sample$ref$labels[1:1000], 2)
dot.srt <- setup.srt(dot.sample$srt$counts, dot.sample$srt$coordinates)
#> Computing spatial radius
dot <- create.DOT(dot.srt, dot.ref)
```
