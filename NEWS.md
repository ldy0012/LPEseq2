# LPEseq2 0.0.1

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
