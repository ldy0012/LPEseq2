# Generate pseudo-bulk RNA-seq count matrices from an annotated Seurat object

This function aggregates cell-level raw count data from a Seurat object
into pseudo-bulk RNA-seq count matrices. Aggregation can be performed
either by condition and cell type, or by sample, condition, and cell
type.

## Usage

``` r
make_pseudobulk_from_seurat(
  seurat_obj,
  method = c("sample", "condition"),
  assay = "RNA",
  counts_layer = "counts",
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "celltype",
  min_cells = 10,
  split_by_celltype = TRUE,
  id_sep = "__",
  normalize_by_ncells = FALSE,
  ncells_scale = 1000
)
```

## Arguments

- seurat_obj:

  A Seurat object containing raw counts and cell-level metadata.

- method:

  Character string specifying the pseudo-bulk aggregation method. One of
  `"sample"` or `"condition"`. `"sample"` aggregates counts by sample,
  condition, and cell type. `"condition"` aggregates counts by condition
  and cell type.

- assay:

  Character string specifying the assay to use. Default is `"RNA"`.

- counts_layer:

  Character string specifying the count layer or slot to use. Default is
  `"counts"`. For Seurat v5, this is passed as `layer`; for Seurat v4,
  it is used as `slot`.

- sample_col:

  Character string specifying the metadata column containing sample IDs.
  Required when `method = "sample"`.

- condition_col:

  Character string specifying the metadata column containing
  experimental condition or group information.

- celltype_col:

  Character string specifying the metadata column containing cell type
  annotations.

- min_cells:

  Integer specifying the minimum number of cells required for a
  pseudo-bulk sample to be retained. Default is 10.

- split_by_celltype:

  Logical. If `TRUE`, the output also includes count matrices and
  metadata split by cell type. Default is `TRUE`.

- id_sep:

  Character string used to concatenate metadata fields when creating
  pseudo-bulk sample IDs. Default is `"__"`.

- normalize_by_ncells:

  Logical. If `TRUE`, each pseudo-bulk sample's raw count vector is
  divided by its cell count and multiplied by `ncells_scale` after
  aggregation. This scales each sample to a per-`ncells_scale`-cells
  basis, reducing the effect of cell number imbalance across samples.
  Mathematically, for sample s with n_s cells, the normalized count is:
  \$\$count\_{gs}^{norm} = count\_{gs} / (n_s / scale)\$\$ This is
  equivalent to scaling each cell's raw count by `ncells_scale / n_s`
  before aggregation. Because the resulting values are non-integer,
  `normalize.method = "none"` is recommended in downstream
  [`LPE_preprocess()`](https://ldy0012.github.io/LPEseq2/reference/LPE_preprocess.md)
  calls. TMM and DESeq2 normalization require raw integer counts and
  should not be used after this normalization. Only meaningful when
  `method = "sample"`, because cell number imbalance across samples is
  the target problem. When `method = "condition"`, all cells within a
  condition are pooled regardless of sample origin, so per-sample cell
  count normalization does not apply in the same way. Default is
  `FALSE`.

- ncells_scale:

  Positive numeric. Scale factor used when `normalize_by_ncells = TRUE`.
  The pseudo-bulk count for each sample is divided by
  `n_cells / ncells_scale`. Default is `1000` (per-1000-cells
  normalization).

## Value

A list containing:

- `counts`:

  A gene-by-pseudo-sample count matrix. When
  `normalize_by_ncells = FALSE`, this is a sparse integer matrix. When
  `normalize_by_ncells = TRUE`, this is a dense numeric matrix with
  non-integer values.

- `metadata`:

  Pseudo-bulk sample metadata including condition, cell type, number of
  cells, library size, pseudo-bulk method, and normalization flags.

- `counts_by_celltype`:

  If `split_by_celltype = TRUE`, a list of gene-by-pseudo-sample count
  matrices split by cell type.

- `metadata_by_celltype`:

  If `split_by_celltype = TRUE`, a list of pseudo-bulk metadata tables
  split by cell type.

## Details

The function assumes that cell type annotation has already been
completed and stored in the Seurat object's metadata. Raw counts are
extracted from the specified assay and count layer/slot, then summed
across cells belonging to each pseudo-bulk group.

This function aggregates raw counts, not normalized or log-normalized
expression values. The resulting pseudo-bulk count matrix is intended to
be used as input for bulk RNA-seq style normalization and differential
expression analysis.

When `method = "condition"`, biological replicates are pooled within
each condition and cell type. This usually produces one pseudo-bulk
sample per condition and cell type, and is therefore mainly suitable for
exploratory comparison.

When `method = "sample"`, biological replicates are preserved by
generating one pseudo-bulk sample per sample and cell type. This is the
recommended mode for downstream differential expression analysis.

When `normalize_by_ncells = TRUE`, the pseudo-bulk count matrix is
scaled so that each sample represents a per-`ncells_scale`-cells
expression level. This removes the effect of unequal cell numbers across
samples, which can otherwise inflate library size differences and
distort normalization. The scaling is applied after aggregation and
after the `min_cells` filter. Because the scaled values are non-integer,
TMM and DESeq2 normalization are not appropriate in downstream steps.
Use `normalize.method = "none"` with `log.transform = TRUE` in
[`LPE_preprocess()`](https://ldy0012.github.io/LPEseq2/reference/LPE_preprocess.md).
This option is only meaningful when `method = "sample"`.

## Examples

``` r
if (FALSE) { # \dontrun{
# Standard pseudo-bulk (raw integer counts)
pb <- make_pseudobulk_from_seurat(
  seurat_obj = seurat_obj,
  method = "sample",
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "celltype"
)

neuron_counts <- pb$counts_by_celltype$Neuron
neuron_meta <- pb$metadata_by_celltype$Neuron

# Cell-number-normalized pseudo-bulk
# Use when cell numbers are substantially unequal across samples
pb_norm <- make_pseudobulk_from_seurat(
  seurat_obj = seurat_obj,
  method = "sample",
  sample_col = "sample_id",
  condition_col = "condition",
  celltype_col = "celltype",
  normalize_by_ncells = TRUE,
  ncells_scale = 1000
)
} # }
```
