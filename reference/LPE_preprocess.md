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

  A one-way design formula, such as `~ group`.

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

A list containing normalized expression matrix, group factor, design,
colData, and preprocessing options.
