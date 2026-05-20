# Run LPE-ANOVA

Performs multi-group differential expression analysis using an
intensity-dependent local pooled error variance trend.

## Usage

``` r
LPE_ANOVA(
  object,
  n.bin = 100,
  df = 10,
  trim.method = c("fixed", "local_fixed", "none"),
  d = 1.2,
  local.k = 3,
  min.local.bin.size = 10,
  use_weighted_between = FALSE,
  verbose = TRUE,
  p.method = c("chisq", "F_inf")
)
```

## Arguments

- object:

  A list returned by
  [`LPE_preprocess()`](https://ldy0012.github.io/LPEseq2/reference/LPE_preprocess.md).

- n.bin:

  Number of quantile bins used for variance trend estimation.

- df:

  Degrees of freedom for smoothing spline.

- trim.method:

  Outlier trimming method applied only to between-group differences. One
  of `"fixed"`, `"local_fixed"`, or `"none"`.

- d:

  Fixed trimming threshold on the raw log2-scale between-group
  difference.

- local.k:

  Multiplier for the local MAD-based threshold used when
  `trim.method = "local_fixed"`.

- min.local.bin.size:

  Minimum number of within-group differences required in an A-bin to
  estimate a local threshold.

- use_weighted_between:

  Logical. Whether to include weighted between-group differences in
  variance trend estimation.

- verbose:

  Logical. Whether to print progress messages.

- p.method:

  P-value calculation method. One of `"chisq"` or `"F_inf"`.

## Value

A data.frame containing gene-level LPE-ANOVA statistics. Trimming
information is stored in `attr(result, "trim.info")`.
