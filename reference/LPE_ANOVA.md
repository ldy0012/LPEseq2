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
  variance.eval = c("grand_mean", "per_group")
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
stored in `attr(result, "variance.settings")`.

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
