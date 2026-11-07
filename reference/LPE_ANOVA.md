# Run LPE-ANOVA

Performs multi-group (and multi-way, interaction-free) differential
expression analysis using an intensity-dependent local pooled error
variance trend.

## Usage

``` r
LPE_ANOVA(
  object,
  n.bin = 100,
  df = 10,
  trim.method = c("auto", "dvalue", "iqr", "none"),
  use_weighted_between = "auto",
  d.threshold = 1.2,
  analysis.method = c("LPE", "standard_anova", "auto"),
  standard.min.group.n = 5,
  verbose = TRUE,
  p.method = c("chisq", "F_inf"),
  variance.eval = c("grand_mean", "per_group"),
  variance.trend = c("pooled", "per_group")
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
  trend estimation. One of `"auto"` (default), `"dvalue"`, `"iqr"`, or
  `"none"`. `"auto"` chooses by replication (see *Automatic
  variance-training settings* in Details): `"none"` when every design
  cell has at least two samples, otherwise `"dvalue"`. See
  [`LPE_ANOVA_var`](https://ldy0012.github.io/LPEseq2/reference/LPE_ANOVA_var.md)
  for the trimming rules.

- use_weighted_between:

  `"auto"` (default), `TRUE` or `FALSE`. Whether to include weighted
  between-group differences in variance trend estimation. `"auto"` uses
  `FALSE` when every design cell has at least two samples (replicate
  pairs only), otherwise `TRUE`.

- d.threshold:

  Numeric. Threshold on the LPEseq1 raw-D scale, used to derive the
  fixed M-scale cutoff (`d.threshold/sqrt(2)`) applied when
  `trim.method = "dvalue"`. Default 1.2, matching the default reported
  in LPEseq1 (Gim et al. 2016). Ignored for other `trim.method` values.
  See
  [`LPE_ANOVA_var`](https://ldy0012.github.io/LPEseq2/reference/LPE_ANOVA_var.md)
  for details.

- analysis.method:

  Analysis method. One of `"LPE"`, `"standard_anova"`, or `"auto"`.
  `"LPE"` uses the local pooled error-based ANOVA. `"standard_anova"`
  uses conventional gene-wise ANOVA (Type II sums of squares for
  multi-way designs). `"auto"` selects standard ANOVA when every
  factor-level combination (design "cell") has at least
  `standard.min.group.n` samples; otherwise, LPE-ANOVA is used.

- standard.min.group.n:

  Minimum per-cell sample size required to use standard ANOVA when
  `analysis.method = "auto"`. A "cell" is a unique combination of the
  design factor levels (for a one-way design, this is simply a group).

- verbose:

  Logical. Whether to print progress messages.

- p.method:

  P-value calculation method. One of `"chisq"` or `"F_inf"`.

- variance.eval:

  Character string specifying how the fitted variance trend is evaluated
  for the test statistic. One of `"grand_mean"` or `"per_group"`.
  `"grand_mean"` (default) evaluates the variance trend once at each
  gene's overall (grand) mean expression and uses it as a single pooled
  variance for the test statistic; this mode supports both one-way and
  multi-way (main-effects-only) designs. `"per_group"` evaluates the
  variance trend separately at each group's mean expression and combines
  the group-specific variances via inverse-variance weighting into a
  Welch-type unequal-variance test statistic (returned with
  `method = "LPE_Welch"` and a different `var`/`MS_between` definition;
  see Details). `"per_group"` is only supported for one-way designs (a
  single design term); it raises an error for multi-way designs.

- variance.trend:

  Character string specifying which variance trend is read for each
  group when `variance.eval = "per_group"`. One of `"pooled"` (default)
  or `"per_group"`. `"pooled"` reads every group's variance from the
  single trend estimated from all design cells (behaviour of previous
  versions). `"per_group"` fits a separate variance trend for every
  group with at least two samples, using only that group's replicate
  pairs, as in LPEseq1 and the original LPE method; groups with a single
  sample read their variance from the pooled trend. Only used with
  `variance.eval = "per_group"`; requesting `"per_group"` together with
  `variance.eval = "grand_mean"` raises an error. See Details.

## Value

A data.frame containing gene-level test statistics. For a one-way
design, columns are unsuffixed (`MS_between`, `F`, `p.value`, `q.value`,
...), matching previous versions. For a multi-way design (e.g.
`design = ~ A + B`), one main-effect test is reported per design term,
with columns suffixed by the term name (e.g. `F_A`, `p.value_A`,
`q.value_A`, `F_B`, `p.value_B`, `q.value_B`, ...); `gene`, `mean`,
`var`, and `method` are shared across terms. The selected analysis
method is stored in `attr(result, "analysis.method")`, and the design
term names are stored in `attr(result, "design.terms")`. For LPE-ANOVA,
trimming information, variance trend information, bin-level variance
points, and the fitted variance trend object are stored in
`attr(result, "trim.info")`, `attr(result, "trend.info")`,
`attr(result, "base.var")`, and `attr(result, "var.spline")`. The
variance-training settings actually used (after resolving `"auto"`) are
stored in `attr(result, "variance.settings")`. With
`variance.trend = "per_group"`, the group-specific trends and their
bin-level variance points are stored in
`attr(result, "var.spline.group")` and `attr(result, "base.var.group")`
(named by group; `NULL` for groups that used the pooled trend), and the
trend used by each group is recorded in
`attr(result, "variance.settings")$trend.source`.

## Details

**Automatic variance-training settings.** With the defaults
`trim.method = "auto"` and `use_weighted_between = "auto"`, the pairwise
values used to estimate the variance trend are chosen from the
replication of the design cells (unique combinations of the design
factor levels):

- every cell has at least two samples: within-cell replicate pairs only,
  without trimming (`use_weighted_between = FALSE`,
  `trim.method = "none"`). Replicate pairs contain no treatment effect,
  so trimming them would only remove genuine noise and under-estimate
  the variance.

- at least one cell has a single sample: within-cell pairs (if any) plus
  weighted between-group pairs with fixed-threshold trimming
  (`use_weighted_between = TRUE`, `trim.method = "dvalue"`), the
  non-replicate procedure of LPEseq1 (Gim et al. 2016).

This mirrors LPEseq1, which uses within-group replicate pairs without
outlier removal whenever a group is replicated and the outlier-trimmed
between-condition pairs only for non-replicated groups. A value given
explicitly for either argument is always used as given.

**Group-specific variance trends.** With `variance.eval = "per_group"`
and `variance.trend = "per_group"` (one-way designs only), the variance
of group \\i\\ is read from a trend \\\sigma^2_i(A)\\ fitted to the
replicate pairs of that group alone, evaluated at the group mean:
\$\$w\_{gi} = n_i / \hat\sigma^2_i(\bar y\_{gi}), \qquad T_g = \sum_i
w\_{gi} (\bar y\_{gi} - \bar y^w_g)^2 \sim \chi^2\_{k-1}.\$\$ This
follows LPEseq1 (Gim et al. 2016) and the original LPE method (Jain et
al. 2003), which estimate a baseline error curve for each replicated
condition, and allows the noise level to differ between conditions at
the same expression intensity. The rules are:

- a group with at least two samples uses its own trend, fitted with
  within-group pairs only. Trimming follows `trim.method` when it is
  given explicitly; with `trim.method = "auto"` no trimming is applied,
  because replicate pairs contain no treatment effect;

- a group with a single sample uses the pooled trend (estimated as for
  `variance.trend = "pooled"`, i.e. with weighted between-group pairs
  and D-value trimming under the automatic settings), as LPEseq1 does
  for non-replicated conditions;

- if a group-specific trend cannot be fitted (too few valid pairwise
  values), that group falls back to the pooled trend with a warning.

When every group has a single sample, all groups use the pooled trend
and the result equals `variance.trend = "pooled"`. The test statistic
itself is unchanged; only the trend each group's variance is read from
differs.

**Multi-way designs.** For a multi-way design `~ A + B + ...` (no
interaction terms allowed), each design term's main effect is tested by
comparing, per gene, the full additive model against the model with that
term dropped: `SS_term = RSS(reduced) - RSS(full)`, with
`df_term = rank(X_full) - rank(X_reduced)`. This is the standard Type II
sum-of-squares decomposition, which for a main-effects-only model does
not depend on term order. The intensity-dependent pooled variance (from
the same LPE variance-trend spline used for one-way designs) is used as
the common denominator for every term's test statistic. This reduces
exactly to the original one-way LPE-ANOVA when the design has a single
term.

## Examples

``` r
testdat <- system.file("testdat", package = "LPEseq2")
read_counts <- function(f) as.matrix(read.csv(file.path(testdat, f),
                                              row.names = 1, check.names = FALSE))
read_meta <- function(f) read.csv(file.path(testdat, f), row.names = 1)

## Two groups (6 vs 6), one-way design
prep2 <- LPE_preprocess(read_counts("LPEseq2_2group_test_counts.csv"),
                        read_meta("LPEseq2_2group_test_metadata.csv"),
                        design = ~ group, verbose = FALSE)
res2 <- LPE_ANOVA(prep2, verbose = FALSE)
head(res2[order(res2$p.value), ])
#>         gene     mean       var MS_between        F      p.value      q.value
#> 12 gene_0012 7.073082 0.2899427   18.05682 62.27719 2.983604e-15 2.983604e-12
#> 19 gene_0019 7.802183 0.2598390   13.35244 51.38738 7.582448e-13 3.791224e-10
#> 38 gene_0038 7.117452 0.2878370   14.11076 49.02344 2.529213e-12 8.430709e-10
#> 22 gene_0022 7.010981 0.2934340   13.66770 46.57846 8.802260e-12 2.200565e-09
#> 11 gene_0011 7.820358 0.2589403   11.78999 45.53169 1.501890e-11 3.003780e-09
#> 24 gene_0024 7.448984 0.2759103   10.95124 39.69129 2.974497e-10 4.957495e-08
#>    method
#> 12    LPE
#> 19    LPE
#> 38    LPE
#> 22    LPE
#> 11    LPE
#> 24    LPE
sum(res2$q.value < 0.05)
#> [1] 82
attr(res2, "variance.settings")$rule
#> [1] "within-group (replicate) pairs only, no trimming (smallest cell n = 6)"

## Welch-type statistic with group-specific variance trends
res2_pg <- LPE_ANOVA(prep2, variance.eval = "per_group",
                     variance.trend = "per_group", verbose = FALSE)
attr(res2_pg, "variance.settings")$trend.source
#>     Control Treatment_A 
#>       "own"       "own" 

## Three groups with a batch factor, additive design ~ group + batch
prep3 <- LPE_preprocess(read_counts("LPEseq2_test_counts.csv"),
                        read_meta("LPEseq2_test_metadata.csv"),
                        design = ~ group + batch, verbose = FALSE)
res3 <- LPE_ANOVA(prep3, verbose = FALSE)
head(res3[order(res3$p.value_group), ])
#>         gene     mean       var MS_between_group  F_group p.value_group
#> 11 gene_0011 7.383532 0.2868123         9.580432 33.40314  3.113159e-15
#> 12 gene_0012 6.888715 0.3068180         9.391096 30.60804  5.094461e-14
#> 38 gene_0038 6.802491 0.3119396         8.997715 28.84442  2.971858e-13
#> 19 gene_0019 7.543674 0.2804833         7.886401 28.11719  6.149783e-13
#> 22 gene_0022 6.711890 0.3181820         8.593651 27.00860  1.863438e-12
#> 25 gene_0025 8.482922 0.2395040         6.363551 26.56971  2.890156e-12
#>    q.value_group MS_between_batch    F_batch p.value_batch q.value_batch method
#> 11  3.113159e-12       0.21966581 0.76588690    0.38149249     0.9983678    LPE
#> 12  2.547231e-11       0.02524831 0.08229086    0.77421659     0.9983678    LPE
#> 38  9.906192e-11       0.07251272 0.23245758    0.62970753     0.9983678    LPE
#> 19  1.537446e-10       0.10840430 0.38649114    0.53414960     0.9983678    LPE
#> 22  3.726877e-10       0.07196450 0.22617398    0.63437539     0.9983678    LPE
#> 25  4.816926e-10       0.75357457 3.14639711    0.07609481     0.7939179    LPE

## Expected results shipped with the package
expected <- read.csv(file.path(testdat, "LPEseq2_2group_expected_results.csv"))
all.equal(res2$p.value, expected$p.value)
#> [1] TRUE
```
