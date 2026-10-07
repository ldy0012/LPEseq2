# LPEseq2 0.0.1

* Shiny web tool (`app.R`): the single "Group variable" dropdown is replaced by a
  dynamic set of dropdowns controlled by "+" / "-" buttons, so several
  testing-group factors can be added to the design (`~ A + B + ...`, no
  interaction). Selections are preserved when a factor is added, a new
  dropdown defaults to the first unused metadata column, duplicate / NA /
  single-level selections give a clear message, and "Variance evaluation" is
  restricted to `grand_mean` when more than one factor is chosen. The design
  formula is built with backquoted names, so metadata columns containing
  spaces now work (also in one-way analyses).

* `LPE_preprocess()` / `LPE_ANOVA()` / `standard_ANOVA_expr()`: design variables
  with non-syntactic names (e.g. `time point`) are supported in multi-way
  formulas via backquotes; output columns use the plain name (`F_time point`).

* `LPE_preprocess()` / `LPE_ANOVA()` / `standard_ANOVA_expr()`: added support
  for multi-way (main-effects-only) designs. `design` can now be a formula
  with more than one term, e.g. `~ genotype + treatment`, in addition to a
  one-way `~ group`. Interaction terms (`A:B` or `A*B`) are not supported
  and raise an error, since the LPE variance-trend / test-statistic
  machinery assumes an additive model. Each design term's main effect is
  tested via a Type II sum-of-squares decomposition (RSS of the model with
  that term dropped, minus RSS of the full model), which is
  order-independent for an additive design and reduces exactly to the
  original one-way computation when there is a single term (verified
  numerically identical). For a one-way design, `LPE_ANOVA()` output columns
  are unchanged (`MS_between`, `F`, `p.value`, `q.value`, ...); for a
  multi-way design, one suffixed set of these columns is reported per
  design term (e.g. `F_genotype`, `p.value_genotype`, `F_treatment`,
  `p.value_treatment`, ...). `LPE_preprocess()`'s returned list gains a
  `factors` field (named list of the individual design factors); `group`
  is kept for backward compatibility and now represents the crossing of
  all design factors (identical to the single factor for one-way designs).
  `variance.eval = "per_group"` (Welch-style) remains one-way only and
  raises an informative error for multi-way designs;
  `analysis.method = "auto"` now bases its group-size check on the minimum
  factor-level combination ("cell") size. `LPE_ANOVA_var()` is unchanged.

* `LPE_ANOVA_var()`: pairwise D/M values used for intensity-dependent
  variance trend estimation are now symmetrized to ±D before trimming and
  quantile binning, matching LPEseq1's original approach (Gim et al. 2016)
  of fixing the difference distribution's center at zero. This changes the
  variance estimation formula from a data-estimated-mean variance to a
  known-mean (zero) variance, which will shift `var`, `F`, `p.value`, and
  `q.value` in `LPE_ANOVA()` results (LPE-ANOVA mode) compared to previous
  versions. `trim.info` pairwise counts (`n_total_before`,
  `n_within_before`, `n_between_before`, and their `_after`/`_removed`
  counterparts) are now twice the number of underlying pairs.

* `LPE_pseudobulk()` now supports `trim.method = "dvalue"` and a
  `d.threshold` argument, matching `LPE_ANOVA()` / `LPE_ANOVA_var()`.
  Previously only `"iqr"` and `"none"` were available at the pseudobulk
  level even though `"dvalue"` was already implemented in the underlying
  `LPE_ANOVA_var()`.

* `LPE_ANOVA()` / `LPE_ANOVA_var()` / `LPE_pseudobulk()`: `trim.method`
  default changed from `"iqr"` to `"dvalue"`.

* `LPE_ANOVA()` / `LPE_ANOVA_var()` / `LPE_pseudobulk()`: `use_weighted_between`
  default changed from `FALSE` to `TRUE`.

* `LPE_ANOVA()` gained a `variance.eval` argument (`"grand_mean"` (default)
  or `"per_group"`), controlling whether the intensity-dependent variance is
  evaluated once at the grand mean or separately per group (Welch-style,
  equivalent to LPEseq1's Z-test at k = 2).

* `LPE_ANOVA_var()`: `trim.method = "dvalue"` now applies its threshold on
  the M scale (`d.threshold / sqrt(2)`) instead of the raw D scale. Under
  the previous raw-D thresholding, the effective strictness of outlier
  removal for between-group pairs became more lenient as group sample size
  increased; thresholding on M keeps the cutoff equally strict regardless
  of group size, and is numerically identical to the previous behavior when
  all groups have exactly one sample (LPEseq1's original design).
