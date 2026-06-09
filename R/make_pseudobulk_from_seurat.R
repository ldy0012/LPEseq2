#' Generate pseudo-bulk RNA-seq count matrices from an annotated Seurat object
#'
#' This function aggregates cell-level raw count data from a Seurat object into
#' pseudo-bulk RNA-seq count matrices. Aggregation can be performed either by
#' condition and cell type, or by sample, condition, and cell type.
#'
#' The function assumes that cell type annotation has already been completed and
#' stored in the Seurat object's metadata. Raw counts are extracted from the
#' specified assay and count layer/slot, then summed across cells belonging to
#' each pseudo-bulk group.
#'
#' @param seurat_obj A Seurat object containing raw counts and cell-level metadata.
#' @param method Character string specifying the pseudo-bulk aggregation method.
#'   One of \code{"sample"} or \code{"condition"}.
#'   \code{"sample"} aggregates counts by sample, condition, and cell type.
#'   \code{"condition"} aggregates counts by condition and cell type.
#' @param assay Character string specifying the assay to use. Default is \code{"RNA"}.
#' @param counts_layer Character string specifying the count layer or slot to use.
#'   Default is \code{"counts"}. For Seurat v5, this is passed as \code{layer};
#'   for Seurat v4, it is used as \code{slot}.
#' @param sample_col Character string specifying the metadata column containing
#'   sample IDs. Required when \code{method = "sample"}.
#' @param condition_col Character string specifying the metadata column containing
#'   experimental condition or group information.
#' @param celltype_col Character string specifying the metadata column containing
#'   cell type annotations.
#' @param min_cells Integer specifying the minimum number of cells required for a
#'   pseudo-bulk sample to be retained. Default is 10.
#' @param split_by_celltype Logical. If \code{TRUE}, the output also includes
#'   count matrices and metadata split by cell type. Default is \code{TRUE}.
#' @param id_sep Character string used to concatenate metadata fields when creating
#'   pseudo-bulk sample IDs. Default is \code{"__"}.
#'
#' @details
#' This function aggregates raw counts, not normalized or log-normalized
#' expression values. The resulting pseudo-bulk count matrix is intended to be
#' used as input for bulk RNA-seq style normalization and differential expression
#' analysis.
#'
#' When \code{method = "condition"}, biological replicates are pooled within each
#' condition and cell type. This usually produces one pseudo-bulk sample per
#' condition and cell type, and is therefore mainly suitable for exploratory
#' comparison.
#'
#' When \code{method = "sample"}, biological replicates are preserved by generating
#' one pseudo-bulk sample per sample and cell type. This is the recommended mode
#' for downstream differential expression analysis.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{counts}}{A gene-by-pseudo-sample sparse raw count matrix.}
#'   \item{\code{metadata}}{Pseudo-bulk sample metadata including condition,
#'   cell type, number of cells, library size, and pseudo-bulk method.}
#'   \item{\code{counts_by_celltype}}{If \code{split_by_celltype = TRUE}, a list
#'   of gene-by-pseudo-sample count matrices split by cell type.}
#'   \item{\code{metadata_by_celltype}}{If \code{split_by_celltype = TRUE}, a list
#'   of pseudo-bulk metadata tables split by cell type.}
#' }
#'
#' @examples
#' \dontrun{
#' pb <- make_pseudobulk_from_seurat(
#'   seurat_obj = seurat_obj,
#'   method = "sample",
#'   sample_col = "sample_id",
#'   condition_col = "condition",
#'   celltype_col = "celltype"
#' )
#'
#' neuron_counts <- pb$counts_by_celltype$Neuron
#' neuron_meta <- pb$metadata_by_celltype$Neuron
#' }
#'
#' @export

make_pseudobulk_from_seurat <- function(
    seurat_obj,
    method = c("sample", "condition"),
    assay = "RNA",
    counts_layer = "counts",
    sample_col = "sample_id",
    condition_col = "condition",
    celltype_col = "celltype",
    min_cells = 10,
    split_by_celltype = TRUE,
    id_sep = "__"
) {
  method <- match.arg(method)
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required.")
  }
  
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix package is required.")
  }
  
  # ------------------------------------------------------------
  # 1. Extract the raw count matrix
  #    Seurat v5 uses layer = "counts"
  #    Seurat v4 uses slot = "counts"
  # ------------------------------------------------------------
  counts <- tryCatch(
    {
      Seurat::GetAssayData(
        object = seurat_obj,
        assay = assay,
        layer = counts_layer
      )
    },
    error = function(e) {
      Seurat::GetAssayData(
        object = seurat_obj,
        assay = assay,
        slot = counts_layer
      )
    }
  )
  
  if (is.null(colnames(counts))) {
    stop("Count matrix must have cell barcodes as column names.")
  }
  
  # ------------------------------------------------------------
  # 2. Extract cell-level metadata
  # ------------------------------------------------------------
  meta <- seurat_obj[[]]
  
  required_cols <- c(condition_col, celltype_col)
  
  if (method == "sample") {
    required_cols <- c(sample_col, condition_col, celltype_col)
  }
  
  missing_cols <- setdiff(required_cols, colnames(meta))
  
  if (length(missing_cols) > 0) {
    stop(
      "These metadata columns are missing: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  # Match the cell order between the count matrix and metadata
  common_cells <- intersect(colnames(counts), rownames(meta))
  
  if (length(common_cells) == 0) {
    stop("No overlapping cells between count matrix and metadata.")
  }
  
  counts <- counts[, common_cells, drop = FALSE]
  meta <- meta[common_cells, , drop = FALSE]
  
  # ------------------------------------------------------------
  # 3. Remove cells with missing metadata values
  # ------------------------------------------------------------
  keep_cells <- complete.cases(meta[, required_cols, drop = FALSE])
  
  counts <- counts[, keep_cells, drop = FALSE]
  meta <- meta[keep_cells, , drop = FALSE]
  
  if (ncol(counts) == 0) {
    stop("No cells remain after removing cells with missing metadata.")
  }
  
  # ------------------------------------------------------------
  # 4. Define pseudo-bulk groups
  #    method = "condition":
  #      aggregate raw counts by condition and cell type
  #    method = "sample":
  #      aggregate raw counts by sample, condition, and cell type
  # ------------------------------------------------------------
  if (method == "condition") {
    group_df <- data.frame(
      condition = as.character(meta[[condition_col]]),
      celltype = as.character(meta[[celltype_col]]),
      stringsAsFactors = FALSE
    )
    
    group_id <- paste(
      group_df$condition,
      group_df$celltype,
      sep = id_sep
    )
  }
  
  if (method == "sample") {
    group_df <- data.frame(
      sample_id = as.character(meta[[sample_col]]),
      condition = as.character(meta[[condition_col]]),
      celltype = as.character(meta[[celltype_col]]),
      stringsAsFactors = FALSE
    )
    
    group_id <- paste(
      group_df$sample_id,
      group_df$condition,
      group_df$celltype,
      sep = id_sep
    )
  }
  
  group_factor <- factor(group_id, levels = unique(group_id))
  
  # ------------------------------------------------------------
  # 5. Create a sparse cell-to-pseudobulk design matrix
  #    Rows represent individual cells
  #    Columns represent pseudo-bulk samples
  # ------------------------------------------------------------
  design_mat <- Matrix::sparseMatrix(
    i = seq_along(group_factor),
    j = as.integer(group_factor),
    x = 1,
    dims = c(length(group_factor), nlevels(group_factor))
  )
  
  colnames(design_mat) <- levels(group_factor)
  rownames(design_mat) <- colnames(counts)
  
  # ------------------------------------------------------------
  # 6. Aggregate raw counts
  #    gene x cell  %*%  cell x pseudo-sample
  #    = gene x pseudo-sample
  # ------------------------------------------------------------
  pseudo_counts <- counts %*% design_mat
  pseudo_counts <- as(pseudo_counts, "dgCMatrix")
  
  # ------------------------------------------------------------
  # 7. Create pseudo-bulk sample metadata
  # ------------------------------------------------------------
  pseudo_meta <- group_df[!duplicated(group_id), , drop = FALSE]
  pseudo_meta$pseudo_sample_id <- unique(group_id)
  
  n_cells <- as.integer(table(group_factor)[pseudo_meta$pseudo_sample_id])
  pseudo_meta$n_cells <- n_cells
  
  pseudo_meta$library_size <- as.numeric(Matrix::colSums(pseudo_counts))
  pseudo_meta$method <- method
  
  rownames(pseudo_meta) <- pseudo_meta$pseudo_sample_id
  
  # Match metadata rows to the pseudo-count matrix columns
  pseudo_meta <- pseudo_meta[colnames(pseudo_counts), , drop = FALSE]
  
  # ------------------------------------------------------------
  # 8. Filter pseudo-bulk samples by the minimum number of cells
  # ------------------------------------------------------------
  keep_pb <- pseudo_meta$n_cells >= min_cells
  
  pseudo_counts <- pseudo_counts[, keep_pb, drop = FALSE]
  pseudo_meta <- pseudo_meta[keep_pb, , drop = FALSE]
  
  if (ncol(pseudo_counts) == 0) {
    stop("No pseudo-bulk samples remain after applying min_cells filter.")
  }
  
  # ------------------------------------------------------------
  # 9. Optionally split pseudo-bulk matrices by cell type
  # ------------------------------------------------------------
  if (split_by_celltype) {
    celltypes <- unique(pseudo_meta$celltype)
    
    counts_by_celltype <- lapply(celltypes, function(ct) {
      cols <- rownames(pseudo_meta)[pseudo_meta$celltype == ct]
      pseudo_counts[, cols, drop = FALSE]
    })
    
    meta_by_celltype <- lapply(celltypes, function(ct) {
      pseudo_meta[pseudo_meta$celltype == ct, , drop = FALSE]
    })
    
    names(counts_by_celltype) <- celltypes
    names(meta_by_celltype) <- celltypes
    
    return(
      list(
        counts = pseudo_counts,
        metadata = pseudo_meta,
        counts_by_celltype = counts_by_celltype,
        metadata_by_celltype = meta_by_celltype
      )
    )
  }
  
  return(
    list(
      counts = pseudo_counts,
      metadata = pseudo_meta
    )
  )
}