# Estimate Intensity-Dependent Variance Trend for LPE-ANOVA

This function estimates an intensity-dependent variance trend using
within-group pairwise differences. Optionally, weighted between-group
differences can be included when replicate information is limited. When
trimming is enabled, the conventional 1.5\*IQR rule is applied within
expression-intensity bins to pairwise values used for variance trend
estimation.

## Usage

``` r
LPE_ANOVA_var(
  expr,
  group,
  n.bin = 100,
  df = 10,
  trim.method = c("iqr", "none"),
  trend.method = c("mean_spline", "quantile_regression"),
  tau = 0.75,
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

  Outlier trimming method applied to pairwise values used for variance
  trend estimation. One of `"iqr"` or `"none"`. `"iqr"` applies the
  conventional 1.5\*IQR boxplot rule within expression-intensity A-bins
  after pooling within-group and between-group values. `"none"` performs
  no outlier trimming.

- trend.method:

  Method used to fit the variance trend. One of `"mean_spline"` or
  `"quantile_regression"`.

- tau:

  Quantile level used when `trend.method = "quantile_regression"`.
  Default is `0.75`.

- use_weighted_between:

  Logical. Whether to include weighted between-group differences in
  variance trend estimation.

## Value

A list containing `type`, `object`, `x_min`, and `x_max`. When
`trend.method = "mean_spline"` or fallback, `object` is a
`smooth.spline`. When `trend.method = "quantile_regression"`, `object`
is an `rq` object. Trimming information is stored in
`attr(result, "trim.info")`, trend information is stored in
`attr(result, "trend.info")`, and bin-level variance points are stored
in `attr(result, "base.var")`.
