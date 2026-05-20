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
  trim.method = c("fixed", "local_fixed", "none"),
  d = 1.2,
  local.k = 3,
  min.local.bin.size = 10,
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

  Outlier trimming method applied only to between-group differences. One
  of `"fixed"`, `"local_fixed"`, or `"none"`. `"fixed"` excludes
  between-group differences with `|D_between| >= d`. `"local_fixed"`
  uses an A-bin-specific threshold estimated from within-group
  differences.

- d:

  Fixed trimming threshold on the raw log2-scale difference. This
  threshold is applied to `D_between`, not to the variance-scaled
  `M_between`.

- local.k:

  Multiplier for the local MAD-based threshold used when
  `trim.method = "local_fixed"`. The local threshold is computed as
  `max(d, local.k * MAD(D_within_bin))`.

- min.local.bin.size:

  Minimum number of within-group differences required in an A-bin to
  estimate a local threshold. If insufficient, the method falls back to
  the fixed threshold `d`.

- use_weighted_between:

  Logical. Whether to include weighted between-group differences in
  variance trend estimation.

## Value

A `smooth.spline` object representing the estimated variance trend.
Trimming information is stored in `attr(object, "trim.info")`.
