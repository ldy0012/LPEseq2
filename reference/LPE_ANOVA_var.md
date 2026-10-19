# Estimate Intensity-Dependent Variance Trend for LPE-ANOVA

This function estimates an intensity-dependent variance trend using
within-group pairwise differences. Following LPEseq1 (Gim et al. 2016),
each raw pairwise difference D is symmetrized to \\\pm D\\ before
quantile binning, so that the pooled difference distribution used for
variance estimation is centered at zero by construction rather than by
subtracting an empirically estimated mean. Optionally, weighted
between-group differences can be included when replicate information is
limited. When trimming is enabled, the conventional 1.5\*IQR rule is
applied within expression-intensity bins to pairwise values used for
variance trend estimation.

## Usage

``` r
LPE_ANOVA_var(
  expr,
  group,
  n.bin = 100,
  df = 10,
  trim.method = c("dvalue", "iqr", "none"),
  use_weighted_between = TRUE,
  d.threshold = 1.2
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
  trend estimation. One of `"iqr"`, `"dvalue"`, or `"none"`. `"iqr"`
  applies the conventional 1.5\*IQR boxplot rule within
  expression-intensity A-bins after pooling within-group and
  between-group values. `"dvalue"` applies a fixed threshold to the M
  value (the rescaled pairwise difference actually used for variance
  estimation), converted from `d.threshold` as `d.threshold/sqrt(2)`, as
  in LPEseq1's non-replicate outlier procedure (Gim et al. 2016): any
  pairwise value with \|M\| greater than this converted threshold is
  removed, independent of bin and of group size. `"none"` performs no
  outlier trimming.

- use_weighted_between:

  Logical. Whether to include weighted between-group differences in
  variance trend estimation.

- d.threshold:

  Numeric. Threshold on the LPEseq1 raw-D scale, used to derive the
  fixed M-scale cutoff (`d.threshold/sqrt(2)`) applied when
  `trim.method = "dvalue"`. Default 1.2, matching the default reported
  in LPEseq1 (Gim et al. 2016). Ignored for other `trim.method` values.
  See `LPE_ANOVA_var` for details.

## Note

Pairwise D/M values are symmetrized (\\\pm D\\) before trimming and
quantile binning, following LPEseq1's approach of fixing the difference
distribution's center at zero. Consequently, the pairwise counts
reported in the returned `trim.info` (e.g. `n_total_before`,
`n_within_before`, `n_between_before`, and their `_after`/`_removed`
counterparts) are twice the number of underlying replicate/group pairs.
