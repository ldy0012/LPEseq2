# Run LPE-ANOVA

Performs multi-group differential expression analysis using an
intensity-dependent local pooled error variance trend.

## Usage

``` r
LPE_ANOVA(
  object,
  n.bin = 100,
  df = 10,
  trim.method = c("mad", "quantile", "fixed"),
  d = 1.2,
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

  Outlier trimming method. One of `"mad"`, `"quantile"`, or `"fixed"`.

- d:

  Fixed trimming threshold used when `trim.method = "fixed"`.

- use_weighted_between:

  Logical. Whether to include weighted between-group differences in
  variance trend estimation.

- verbose:

  Logical. Whether to print progress messages.

- p.method:

  P-value calculation method. One of `"chisq"` or `"F_inf"`.

## Value

A data.frame containing gene-level LPE-ANOVA statistics.
