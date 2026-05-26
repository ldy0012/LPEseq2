# Run LPE-ANOVA

Performs multi-group differential expression analysis using an
intensity-dependent local pooled error variance trend.

## Usage

``` r
LPE_ANOVA(
  object,
  n.bin = 100,
  df = 10,
  trim.method = c("iqr", "none"),
  trend.method = c("mean_spline", "quantile_regression"),
  tau = 0.75,
  use_weighted_between = FALSE,
  analysis.method = c("LPE", "standard_anova", "auto"),
  standard.min.group.n = 5,
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

  Outlier trimming method applied to pairwise values used for variance
  trend estimation. One of `"iqr"` or `"none"`. `"iqr"` applies the
  conventional 1.5\*IQR boxplot rule within expression-intensity A-bins
  after pooling within-group and between-group values. `"none"` performs
  no outlier trimming.

- trend.method:

  Method used to fit the variance trend. One of `"mean_spline"` or
  `"quantile_regression"`. `"mean_spline"` fits a smoothing spline to
  bin-level variance estimates. `"quantile_regression"` fits an
  upper-quantile regression trend to reduce potential variance
  underestimation.

- tau:

  Quantile level used when `trend.method = "quantile_regression"`.
  Default is `0.75`.

- use_weighted_between:

  Logical. Whether to include weighted between-group differences in
  variance trend estimation.

- analysis.method:

  Analysis method. One of `"LPE"`, `"standard_anova"`, or `"auto"`.
  `"LPE"` uses the local pooled error-based ANOVA. `"standard_anova"`
  uses conventional gene-wise one-way ANOVA. `"auto"` selects standard
  ANOVA when all groups have at least `standard.min.group.n` samples;
  otherwise, LPE-ANOVA is used.

- standard.min.group.n:

  Minimum per-group sample size required to use standard one-way ANOVA
  when `analysis.method = "auto"`.

- verbose:

  Logical. Whether to print progress messages.

- p.method:

  P-value calculation method. One of `"chisq"` or `"F_inf"`.

## Value

A data.frame containing gene-level test statistics. The selected
analysis method is stored in `attr(result, "analysis.method")`. For
LPE-ANOVA, trimming information, variance trend information, and
bin-level variance points are stored in `attr(result, "trim.info")`,
`attr(result, "trend.info")`, and `attr(result, "base.var")`.
