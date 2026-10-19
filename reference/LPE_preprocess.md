# Preprocess RNA-seq Count Data for LPE-ANOVA

This function checks raw count data, filters low-count genes, normalizes
sample-level library size differences, and optionally applies log2
transformation.

## Usage

``` r
LPE_preprocess(
  counts,
  colData,
  design = ~group,
  normalize.method = c("library_size", "TMM", "DESeq2", "none"),
  log.transform = TRUE,
  min.count = 5,
  prior.count = 1,
  verbose = TRUE
)
```

## Arguments

- counts:

  A numeric matrix of raw counts with genes as rows and samples as
  columns.

- colData:

  A data.frame containing sample-level metadata. Row names must match
  column names of counts.

- design:

  A design formula with no left-hand side, e.g. `~ group` or
  `~ genotype + treatment`. One or more main-effect terms are supported
  (multi-way ANOVA without interaction); interaction terms (`A:B` or
  `A*B`) are not supported and will raise an error. Every term must name
  a column of `colData` that can be coerced to a factor with at least
  two levels.

- normalize.method:

  Normalization method. One of `"TMM"`, `"library_size"`, `"DESeq2"`, or
  `"none"`.

- log.transform:

  Logical. Whether to apply log2 transformation.

- min.count:

  Minimum count threshold for low-count filtering.

- prior.count:

  Prior count added before log transformation.

- verbose:

  Logical. Whether to print progress messages.

## Value

A list containing normalized expression matrix, a factor identifying
each unique combination of the design factors (`group`, kept for
backward compatibility with one-way designs), a named list of the
individual design factors (`factors`), design, colData, and
preprocessing options.
