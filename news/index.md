# Changelog

## LPEseq2 0.1.0

- [`LPE_ANOVA()`](https://ldy0012.github.io/LPEseq2/reference/LPE_ANOVA.md)
  /
  [`LPE_pseudobulk()`](https://ldy0012.github.io/LPEseq2/reference/LPE_pseudobulk.md):
  new default `trim.method = "auto"` and
  `use_weighted_between = "auto"`. The variance-training pairs are
  chosen from the replication of the design cells: when every cell has
  at least two samples, only within-cell replicate pairs are used,
  without trimming (`use_weighted_between = FALSE`,
  `trim.method = "none"`); when any cell has a single sample, weighted
  between-group pairs are added and D-value trimming is applied (`TRUE`,
  `"dvalue"`, the previous default). This follows LPEseq1, which uses
  replicate pairs without outlier removal for replicated groups. In
  simulations of cell-line-like data, trimming replicate pairs
  under-estimated the variance (mean FDR about 0.12 at a nominal 0.05
  with two or three replicates), while replicate pairs without trimming
  kept the FDR closer to nominal (about 0.075) with a similar F1.
  **Results with replicates in every group differ from version 0.0.1
  when the defaults are used.** To reproduce 0.0.1, set
  `trim.method = "dvalue", use_weighted_between = TRUE`. Explicitly
  supplied values are used unchanged.

- The settings actually used are returned in
  `attr(res, "variance.settings")` (`trim.method`,
  `use_weighted_between`, `auto`, `replicated`, `min_cell_n`, `rule`);
  [`LPE_pseudobulk()`](https://ldy0012.github.io/LPEseq2/reference/LPE_pseudobulk.md)
  adds a `variance_training` column per cell type.

- Shiny web tool: “Pairwise outlier trimming method” and “Use weighted
  between-group differences” offer “Auto by replication (recommended)”
  as the default (the latter is now a drop-down instead of a check box).
  The resolved setting is shown under Method info and written to the
  downloaded results (`variance_training` column).

## LPEseq2 0.0.1

- Shiny web tool (`app.R`): the single “Group variable” dropdown is
  replaced by a dynamic set of dropdowns controlled by “+” / “-”
  buttons, so several testing-group factors can be added to the design
  (`~ A + B + ...`, no interaction). Selections are preserved when a
  factor is added, a new dropdown defaults to the first unused metadata
  column, duplicate / NA / single-level selections give a clear message,
  and “Variance evaluation” is restricted to `grand_mean` when more than
  one factor is chosen. The design formula is built with backquoted
  names, so metadata columns containing spaces now work (also in one-way
  analyses).

- [`LPE_preprocess()`](https://ldy0012.github.io/LPEseq2/reference/LPE_preprocess.md)
  /
  [`LPE_ANOVA()`](https://ldy0012.github.io/LPEseq2/reference/LPE_ANOVA.md)
  / `standard_ANOVA_expr()`: design variables with non-syntactic names
  (e.g. `time point`) are supported in multi-way formulas via
  backquotes; output columns use the plain name (`F_time point`).

- [`LPE_preprocess()`](https://ldy0012.github.io/LPEseq2/reference/LPE_preprocess.md)
  /
  [`LPE_ANOVA()`](https://ldy0012.github.io/LPEseq2/reference/LPE_ANOVA.md)
  / `standard_ANOVA_expr()`: added support for multi-way
  (main-effects-only) designs. `design` can now be a formula with more
  than one term, e.g. `~ genotype + treatment`, in addition to a one-way
  `~ group`. Interaction terms (`A:B` or `A*B`) are not supported and
  raise an error, since the LPE variance-trend / test-statistic
  machinery assumes an additive model. Each design term’s main effect is
  tested via a Type II sum-of-squares decomposition (RSS of the model
  with that term dropped, minus RSS of the full model), which is
  order-independent for an additive design and reduces exactly to the
  original one-way computation when there is a single term (verified
  numerically identical). For a one-way design,
  [`LPE_ANOVA()`](https://ldy0012.github.io/LPEseq2/reference/LPE_ANOVA.md)
  output columns are unchanged (`MS_between`, `F`, `p.value`, `q.value`,
  …); for a multi-way design, one suffixed set of these columns is
  reported per design term (e.g. `F_genotype`, `p.value_genotype`,
  `F_treatment`, `p.value_treatment`, …).
  [`LPE_preprocess()`](https://ldy0012.github.io/LPEseq2/reference/LPE_preprocess.md)’s
  returned list gains a `factors` field (named list of the individual
  design factors); `group` is kept for backward compatibility and now
  represents the crossing of all design factors (identical to the single
  factor for one-way designs). `variance.eval = "per_group"`
  (Welch-style) remains one-way only and raises an informative error for
  multi-way designs; `analysis.method = "auto"` now bases its group-size
  check on the minimum factor-level combination (“cell”) size.
  [`LPE_ANOVA_var()`](https://ldy0012.github.io/LPEseq2/reference/LPE_ANOVA_var.md)
  is unchanged.

- [`LPE_ANOVA_var()`](https://ldy0012.github.io/LPEseq2/reference/LPE_ANOVA_var.md):
  pairwise D/M values used for intensity-dependent variance trend
  estimation are now symmetrized to ±D before trimming and quantile
  binning, matching LPEseq1’s original approach (Gim et al. 2016) of
  fixing the difference distribution’s center at zero. This changes the
  variance estimation formula from a data-estimated-mean variance to a
  known-mean (zero) variance, which will shift `var`, `F`, `p.value`,
  and `q.value` in
  [`LPE_ANOVA()`](https://ldy0012.github.io/LPEseq2/reference/LPE_ANOVA.md)
  results (LPE-ANOVA mode) compared to previous versions. `trim.info`
  pairwise counts (`n_total_before`, `n_within_before`,
  `n_between_before`, and their `_after`/`_removed` counterparts) are
  now twice the number of underlying pairs.

- [`LPE_pseudobulk()`](https://ldy0012.github.io/LPEseq2/reference/LPE_pseudobulk.md)
  now supports `trim.method = "dvalue"` and a `d.threshold` argument,
  matching
  [`LPE_ANOVA()`](https://ldy0012.github.io/LPEseq2/reference/LPE_ANOVA.md)
  /
  [`LPE_ANOVA_var()`](https://ldy0012.github.io/LPEseq2/reference/LPE_ANOVA_var.md).
  Previously only `"iqr"` and `"none"` were available at the pseudobulk
  level even though `"dvalue"` was already implemented in the underlying
  [`LPE_ANOVA_var()`](https://ldy0012.github.io/LPEseq2/reference/LPE_ANOVA_var.md).

- [`LPE_ANOVA()`](https://ldy0012.github.io/LPEseq2/reference/LPE_ANOVA.md)
  /
  [`LPE_ANOVA_var()`](https://ldy0012.github.io/LPEseq2/reference/LPE_ANOVA_var.md)
  /
  [`LPE_pseudobulk()`](https://ldy0012.github.io/LPEseq2/reference/LPE_pseudobulk.md):
  `trim.method` default changed from `"iqr"` to `"dvalue"`.

- [`LPE_ANOVA()`](https://ldy0012.github.io/LPEseq2/reference/LPE_ANOVA.md)
  /
  [`LPE_ANOVA_var()`](https://ldy0012.github.io/LPEseq2/reference/LPE_ANOVA_var.md)
  /
  [`LPE_pseudobulk()`](https://ldy0012.github.io/LPEseq2/reference/LPE_pseudobulk.md):
  `use_weighted_between` default changed from `FALSE` to `TRUE`.

- [`LPE_ANOVA()`](https://ldy0012.github.io/LPEseq2/reference/LPE_ANOVA.md)
  gained a `variance.eval` argument (`"grand_mean"` (default) or
  `"per_group"`), controlling whether the intensity-dependent variance
  is evaluated once at the grand mean or separately per group
  (Welch-style, equivalent to LPEseq1’s Z-test at k = 2).

- [`LPE_ANOVA_var()`](https://ldy0012.github.io/LPEseq2/reference/LPE_ANOVA_var.md):
  `trim.method = "dvalue"` now applies its threshold on the M scale
  (`d.threshold / sqrt(2)`) instead of the raw D scale. Under the
  previous raw-D thresholding, the effective strictness of outlier
  removal for between-group pairs became more lenient as group sample
  size increased; thresholding on M keeps the cutoff equally strict
  regardless of group size, and is numerically identical to the previous
  behavior when all groups have exactly one sample (LPEseq1’s original
  design).
