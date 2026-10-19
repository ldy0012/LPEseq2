# Run LPEseq2 on pseudo-bulk RNA-seq data generated from a Seurat object

This function creates pseudo-bulk RNA-seq count matrices from an
annotated Seurat object and runs LPEseq2 differential expression
analysis for each cell type. Pseudo-bulk aggregation can be performed
either at the sample level or at the condition level.

## Usage

``` r
LPE_pseudobulk(
  seurat_obj,
  method = c("sample", "condition"),
  assay = "RNA",
  counts_layer = "counts",
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "celltype",
  celltypes = NULL,
  min_cells = 10,
  normalize_by_ncells = FALSE,
  ncells_scale = 1000,
  normalize.method = c("TMM", "library_size", "DESeq2", "none"),
  log.transform = TRUE,
  min.count = 5,
  prior.count = 1,
  analysis.method = c("auto", "LPE", "standard_anova"),
  standard.min.group.n = 5,
  n.bin = 100,
  df = 10,
  trim.method = c("auto", "dvalue", "iqr", "none"),
  d.threshold = 1.2,
  use_weighted_between = "auto",
  p.method = c("chisq", "F_inf"),
  verbose = TRUE
)
```

## Arguments

- seurat_obj:

  A Seurat object containing raw counts and cell-level metadata.

- method:

  Character string specifying the pseudo-bulk aggregation method. One of
  `"sample"` or `"condition"`.

- assay:

  Character string specifying the assay to use. Default is `"RNA"`.

- counts_layer:

  Character string specifying the count layer or slot to use. Default is
  `"counts"`.

- sample_col:

  Character string specifying the metadata column containing sample IDs.
  Required when `method = "sample"`.

- condition_col:

  Character string specifying the metadata column containing
  experimental condition or group information.

- celltype_col:

  Character string specifying the metadata column containing cell type
  annotations.

- celltypes:

  Optional character vector specifying which cell types to analyze. If
  `NULL`, all available cell types are analyzed.

- min_cells:

  Integer specifying the minimum number of cells required for a
  pseudo-bulk sample to be retained. Default is 10.

- normalize_by_ncells:

  Logical. If `TRUE`, each pseudo-bulk sample's aggregated count vector
  is divided by its cell count and multiplied by `ncells_scale`, scaling
  to a per-`ncells_scale`-cells basis. This reduces the effect of cell
  number imbalance across samples. When `TRUE`, `normalize.method` is
  automatically forced to `"none"` and a warning is issued if the user
  specifies `"TMM"` or `"DESeq2"`, because those methods require raw
  integer counts. Only meaningful when `method = "sample"`. Default is
  `FALSE`.

- ncells_scale:

  Positive numeric. Scale factor used when `normalize_by_ncells = TRUE`.
  Default is `1000`.

- normalize.method:

  Character string specifying the normalization method passed to
  [`LPE_preprocess`](https://ldy0012.github.io/LPEseq2/reference/LPE_preprocess.md).
  One of `"TMM"`, `"library_size"`, `"DESeq2"`, or `"none"`. When
  `normalize_by_ncells = TRUE`, this is automatically set to `"none"`.

- log.transform:

  Logical. Whether to apply log transformation during preprocessing.
  Passed to
  [`LPE_preprocess`](https://ldy0012.github.io/LPEseq2/reference/LPE_preprocess.md).

- min.count:

  Numeric value specifying the minimum count threshold used during
  preprocessing. Passed to
  [`LPE_preprocess`](https://ldy0012.github.io/LPEseq2/reference/LPE_preprocess.md).
  When `normalize_by_ncells = TRUE`, the count scale changes and this
  value may need to be adjusted accordingly.

- prior.count:

  Numeric prior count added before log transformation. Passed to
  [`LPE_preprocess`](https://ldy0012.github.io/LPEseq2/reference/LPE_preprocess.md).

- analysis.method:

  Character string specifying the LPEseq2 analysis method. One of
  `"auto"`, `"LPE"`, or `"standard_anova"`.

- standard.min.group.n:

  Integer specifying the minimum group sample size required for standard
  ANOVA when `analysis.method = "auto"`.

- n.bin:

  Integer specifying the number of bins used for intensity-dependent
  variance trend estimation.

- df:

  Numeric value specifying the degrees of freedom for spline smoothing.

- trim.method:

  Character string specifying the trimming method used during variance
  trend estimation. One of `"auto"` (default), `"dvalue"`, `"iqr"`, or
  `"none"`. `"auto"` is resolved separately for each cell type by
  [`LPE_ANOVA`](https://ldy0012.github.io/LPEseq2/reference/LPE_ANOVA.md):
  `"none"` when every design cell has at least two pseudo-bulk samples,
  otherwise `"dvalue"`. `"dvalue"` applies a fixed threshold to the M
  value (the rescaled pairwise difference actually used for variance
  estimation), derived from `d.threshold` as `d.threshold/sqrt(2)`,
  following LPEseq1's non-replicate outlier procedure (Gim et al. 2016).
  See
  [`LPE_ANOVA_var`](https://ldy0012.github.io/LPEseq2/reference/LPE_ANOVA_var.md)
  for details.

- d.threshold:

  Numeric. Threshold on the LPEseq1 raw-D scale, used to derive the
  fixed M-scale cutoff applied when `trim.method = "dvalue"` (also when
  `"auto"` resolves to it). Default 1.2, matching
  [`LPE_ANOVA`](https://ldy0012.github.io/LPEseq2/reference/LPE_ANOVA.md).
  Ignored otherwise.

- use_weighted_between:

  `"auto"` (default), `TRUE` or `FALSE`. Whether to use weighted
  between-group variance information in LPE-ANOVA. `"auto"` uses `FALSE`
  when every design cell has at least two pseudo-bulk samples, otherwise
  `TRUE` (see
  [`LPE_ANOVA`](https://ldy0012.github.io/LPEseq2/reference/LPE_ANOVA.md)).
  The setting used for each cell type is reported in the
  `variance_training` column.

- p.method:

  Character string specifying the p-value calculation method. One of
  `"chisq"` or `"F_inf"`.

- verbose:

  Logical. Whether to print progress messages.

## Value

A list containing:

- `pseudobulk`:

  The full pseudo-bulk object returned by
  [`make_pseudobulk_from_seurat`](https://ldy0012.github.io/LPEseq2/reference/make_pseudobulk_from_seurat.md).

- `counts_by_celltype`:

  A list of pseudo-bulk count matrices used for LPEseq2 analysis, split
  by cell type. When `normalize_by_ncells = TRUE`, these are
  cell-number-normalized non-integer matrices.

- `metadata_by_celltype`:

  A list of pseudo-bulk metadata tables used for LPEseq2 analysis, split
  by cell type.

- `prep_by_celltype`:

  A list of preprocessed LPEseq2 objects returned by
  [`LPE_preprocess`](https://ldy0012.github.io/LPEseq2/reference/LPE_preprocess.md).

- `results_by_celltype`:

  A list of cell-type-specific differential expression result tables
  returned by
  [`LPE_ANOVA`](https://ldy0012.github.io/LPEseq2/reference/LPE_ANOVA.md).

- `combined_result`:

  A combined differential expression result table across all
  successfully analyzed cell types.

- `diagnostics`:

  A list of diagnostic summaries, including the number of pseudo-bulk
  samples, group sizes, cell counts, and library sizes.

- `errors`:

  A list of errors or skipped-analysis messages for cell types that
  could not be analyzed.

- `settings`:

  A list of analysis settings used in the run, including
  `normalize_by_ncells` and `ncells_scale`.

## Details

This function first calls
[`make_pseudobulk_from_seurat`](https://ldy0012.github.io/LPEseq2/reference/make_pseudobulk_from_seurat.md)
to generate raw pseudo-bulk count matrices. Each cell type is then
analyzed separately using
[`LPE_preprocess`](https://ldy0012.github.io/LPEseq2/reference/LPE_preprocess.md)
and
[`LPE_ANOVA`](https://ldy0012.github.io/LPEseq2/reference/LPE_ANOVA.md).

The recommended mode for differential expression analysis is
`method = "sample"`, because it preserves biological replicates. In
contrast, `method = "condition"` pools cells across all samples within
each condition and cell type, usually producing one pseudo-bulk sample
per condition. Results from condition-level pseudo-bulk analysis should
therefore be interpreted as exploratory.

When `normalize_by_ncells = TRUE`, pseudo-bulk counts are scaled to a
per-`ncells_scale`-cells basis before LPEseq2 preprocessing. This
equalizes the contribution of each sample regardless of the number of
cells, which is useful when cell numbers are substantially unequal
across samples. Because the resulting values are non-integer,
`normalize.method` is automatically forced to `"none"` when
`normalize_by_ncells = TRUE`.

## See also

[`make_pseudobulk_from_seurat`](https://ldy0012.github.io/LPEseq2/reference/make_pseudobulk_from_seurat.md),
[`LPE_preprocess`](https://ldy0012.github.io/LPEseq2/reference/LPE_preprocess.md),
[`LPE_ANOVA`](https://ldy0012.github.io/LPEseq2/reference/LPE_ANOVA.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Standard pseudo-bulk (raw integer counts, TMM normalization)
res <- LPE_pseudobulk(
  seurat_obj = seurat_obj,
  method = "sample",
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "celltype",
  normalize.method = "TMM",
  analysis.method = "auto"
)

# Cell-number-normalized pseudo-bulk
# Use when cell numbers are substantially unequal across samples
res_norm <- LPE_pseudobulk(
  seurat_obj = seurat_obj,
  method = "sample",
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "celltype",
  normalize_by_ncells = TRUE,
  ncells_scale = 1000,
  min.count = 1,
  analysis.method = "auto"
)

neuron_res <- res$results_by_celltype$Neuron
all_res <- res$combined_result
} # }
```
