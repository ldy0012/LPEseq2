# Estimate Intensity-Dependent Variance Trend for LPE-ANOVA

This function estimates an intensity-dependent variance trend using
within-group pairwise differences. Optionally, weighted between-group
differences can be included when replicate information is limited.

## Usage

``` r
LPE_ANOVA_var(
  expr,
  group,
  n.bin = 100,
  df = 10,
  trim.method = c("mad", "quantile", "fixed"),
  d = 1.2,
  use_weighted_between = FALSE
)
```

## Arguments

- expr:

  A numeric expression matrix with genes as rows and samples as columns.

- group:

  A factor or vector indicating sample group labels.

- n.bin:

  Number of quantile bins used for variance trend estimation.

- df:

  Degrees of freedom for smoothing spline.

- trim.method:

  Outlier trimming method. One of `"mad"`, `"quantile"`, or `"fixed"`.

- d:

  Fixed trimming threshold used when `trim.method = "fixed"`.

- use_weighted_between:

  Logical. Whether to include weighted between-group differences.

## Value

A `smooth.spline` object representing the estimated variance trend.
