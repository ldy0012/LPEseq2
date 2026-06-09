#' Run LPEseq2 on pseudo-bulk RNA-seq data generated from a Seurat object
#'
#' This function creates pseudo-bulk RNA-seq count matrices from an annotated
#' Seurat object and runs LPEseq2 differential expression analysis for each cell
#' type. Pseudo-bulk aggregation can be performed either at the sample level or
#' at the condition level.
#'
#' @param seurat_obj A Seurat object containing raw counts and cell-level metadata.
#' @param method Character string specifying the pseudo-bulk aggregation method.
#'   One of \code{"sample"} or \code{"condition"}.
#' @param assay Character string specifying the assay to use. Default is \code{"RNA"}.
#' @param counts_layer Character string specifying the count layer or slot to use.
#'   Default is \code{"counts"}.
#' @param sample_col Character string specifying the metadata column containing
#'   sample IDs. Required when \code{method = "sample"}.
#' @param condition_col Character string specifying the metadata column containing
#'   experimental condition or group information.
#' @param celltype_col Character string specifying the metadata column containing
#'   cell type annotations.
#' @param celltypes Optional character vector specifying which cell types to
#'   analyze. If \code{NULL}, all available cell types are analyzed.
#' @param min_cells Integer specifying the minimum number of cells required for a
#'   pseudo-bulk sample to be retained. Default is 10.
#' @param normalize.method Character string specifying the normalization method
#'   passed to \code{\link{LPE_preprocess}}. One of \code{"TMM"},
#'   \code{"library_size"}, \code{"DESeq2"}, or \code{"none"}.
#' @param log.transform Logical. Whether to apply log transformation during
#'   preprocessing. Passed to \code{\link{LPE_preprocess}}.
#' @param min.count Numeric value specifying the minimum count threshold used
#'   during preprocessing. Passed to \code{\link{LPE_preprocess}}.
#' @param prior.count Numeric prior count added before log transformation.
#'   Passed to \code{\link{LPE_preprocess}}.
#' @param analysis.method Character string specifying the LPEseq2 analysis method.
#'   One of \code{"auto"}, \code{"LPE"}, or \code{"standard_anova"}.
#' @param standard.min.group.n Integer specifying the minimum group sample size
#'   required for standard ANOVA when \code{analysis.method = "auto"}.
#' @param n.bin Integer specifying the number of bins used for intensity-dependent
#'   variance trend estimation.
#' @param df Numeric value specifying the degrees of freedom for spline smoothing.
#' @param trim.method Character string specifying the trimming method used during
#'   variance trend estimation. One of \code{"iqr"} or \code{"none"}.
#' @param use_weighted_between Logical. Whether to use weighted between-group
#'   variance information in LPE-ANOVA.
#' @param p.method Character string specifying the p-value calculation method.
#'   One of \code{"chisq"} or \code{"F_inf"}.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @details
#' This function first calls \code{\link{make_pseudobulk_from_seurat}} to generate
#' raw pseudo-bulk count matrices. Each cell type is then analyzed separately
#' using \code{\link{LPE_preprocess}} and \code{\link{LPE_ANOVA}}.
#'
#' The recommended mode for differential expression analysis is
#' \code{method = "sample"}, because it preserves biological replicates. In
#' contrast, \code{method = "condition"} pools cells across all samples within
#' each condition and cell type, usually producing one pseudo-bulk sample per
#' condition. Results from condition-level pseudo-bulk analysis should therefore
#' be interpreted as exploratory.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{pseudobulk}}{The full pseudo-bulk object returned by
#'   \code{\link{make_pseudobulk_from_seurat}}.}
#'   \item{\code{counts_by_celltype}}{A list of raw pseudo-bulk count matrices
#'   used for LPEseq2 analysis, split by cell type.}
#'   \item{\code{metadata_by_celltype}}{A list of pseudo-bulk metadata tables
#'   used for LPEseq2 analysis, split by cell type.}
#'   \item{\code{prep_by_celltype}}{A list of preprocessed LPEseq2 objects
#'   returned by \code{\link{LPE_preprocess}}.}
#'   \item{\code{results_by_celltype}}{A list of cell-type-specific differential
#'   expression result tables returned by \code{\link{LPE_ANOVA}}.}
#'   \item{\code{combined_result}}{A combined differential expression result table
#'   across all successfully analyzed cell types.}
#'   \item{\code{diagnostics}}{A list of diagnostic summaries, including the number
#'   of pseudo-bulk samples, group sizes, cell counts, and library sizes.}
#'   \item{\code{errors}}{A list of errors or skipped-analysis messages for cell
#'   types that could not be analyzed.}
#'   \item{\code{settings}}{A list of analysis settings used in the run.}
#' }
#'
#' @examples
#' \dontrun{
#' res <- LPE_pseudobulk(
#'   seurat_obj = seurat_obj,
#'   method = "sample",
#'   sample_col = "sample_id",
#'   condition_col = "condition",
#'   celltype_col = "celltype",
#'   celltypes = c("Neuron"),
#'   normalize.method = "TMM",
#'   analysis.method = "auto"
#' )
#'
#' neuron_res <- res$results_by_celltype$Neuron
#' all_res <- res$combined_result
#' }
#'
#' @seealso
#' \code{\link{make_pseudobulk_from_seurat}},
#' \code{\link{LPE_preprocess}},
#' \code{\link{LPE_ANOVA}}
#'
#' @export

LPE_pseudobulk <- function(
    seurat_obj,
    method = c("sample", "condition"),
    assay = "RNA",
    counts_layer = "counts",
    sample_col = "sample_id",
    condition_col = "condition",
    celltype_col = "celltype",
    celltypes = NULL,
    min_cells = 10,
    normalize.method = c("TMM", "library_size", "DESeq2", "none"),
    log.transform = TRUE,
    min.count = 5,
    prior.count = 1,
    analysis.method = c("auto", "LPE", "standard_anova"),
    standard.min.group.n = 5,
    n.bin = 100,
    df = 10,
    trim.method = c("iqr", "none"),
    use_weighted_between = FALSE,
    p.method = c("chisq", "F_inf"),
    verbose = TRUE
) {
  method <- match.arg(method)
  normalize.method <- match.arg(normalize.method)
  analysis.method <- match.arg(analysis.method)
  trim.method <- match.arg(trim.method)
  p.method <- match.arg(p.method)
  
  # ------------------------------------------------------------
  # Check whether the required pseudo-bulk generation function exists
  # ------------------------------------------------------------
  if (!exists("make_pseudobulk_from_seurat")) {
    stop(
      "make_pseudobulk_from_seurat() function is not found. ",
      "Please define it before running LPE_pseudobulk()."
    )
  }
  
  # ------------------------------------------------------------
  # Check required packages
  # ------------------------------------------------------------
  if (!requireNamespace("LPEseq2", quietly = TRUE)) {
    stop(
      "LPEseq2 package is required. ",
      "Install it with devtools::install_github('ldy0012/LPEseq2')."
    )
  }
  
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Matrix package is required.")
  }
  
  # ------------------------------------------------------------
  # Warn users when condition-level pooled pseudo-bulk is selected
  # ------------------------------------------------------------
  if (method == "condition") {
    warning(
      "method = 'condition' creates pooled pseudo-bulk samples, usually n = 1 per condition per cell type. ",
      "This approach is useful for exploratory comparison, but DEG p-values should be interpreted cautiously."
    )
  }
  
  # ------------------------------------------------------------
  # 1. Generate pseudo-bulk raw count matrices from the Seurat object
  #    method = "sample":
  #      aggregate raw counts by sample, condition, and cell type
  #    method = "condition":
  #      aggregate raw counts by condition and cell type
  # ------------------------------------------------------------
  pb <- make_pseudobulk_from_seurat(
    seurat_obj = seurat_obj,
    method = method,
    assay = assay,
    counts_layer = counts_layer,
    sample_col = sample_col,
    condition_col = condition_col,
    celltype_col = celltype_col,
    min_cells = min_cells,
    split_by_celltype = TRUE
  )
  
  available_celltypes <- names(pb$counts_by_celltype)
  
  # ------------------------------------------------------------
  # 2. Select cell types to analyze
  #    If celltypes = NULL, all available cell types are analyzed
  # ------------------------------------------------------------
  if (is.null(celltypes)) {
    celltypes_to_run <- available_celltypes
  } else {
    missing_ct <- setdiff(celltypes, available_celltypes)
    
    if (length(missing_ct) > 0) {
      warning(
        "These cell types are not found in the pseudo-bulk object: ",
        paste(missing_ct, collapse = ", ")
      )
    }
    
    celltypes_to_run <- intersect(celltypes, available_celltypes)
  }
  
  if (length(celltypes_to_run) == 0) {
    stop("No valid cell types to analyze.")
  }
  
  # ------------------------------------------------------------
  # 3. Prepare containers for outputs
  # ------------------------------------------------------------
  result_list <- list()
  prep_list <- list()
  count_list <- list()
  metadata_list <- list()
  diagnostic_list <- list()
  error_list <- list()
  
  # ------------------------------------------------------------
  # 4. Run LPEseq2 for each cell type
  # ------------------------------------------------------------
  for (ct in celltypes_to_run) {
    if (verbose) {
      cat("\n==============================\n")
      cat("Running cell type:", ct, "\n")
      cat("==============================\n")
    }
    
    counts_ct <- pb$counts_by_celltype[[ct]]
    meta_ct <- pb$metadata_by_celltype[[ct]]
    
    # ----------------------------------------------------------
    # Convert sparse matrix to a standard numeric matrix
    # LPEseq2 expects a gene x sample count matrix
    # ----------------------------------------------------------
    counts_ct <- as.matrix(counts_ct)
    storage.mode(counts_ct) <- "numeric"
    
    # ----------------------------------------------------------
    # Match metadata rows to count matrix columns
    # ----------------------------------------------------------
    meta_ct <- meta_ct[colnames(counts_ct), , drop = FALSE]
    
    if (!all(colnames(counts_ct) == rownames(meta_ct))) {
      stop(
        "Column names of counts and row names of metadata are not matched for cell type: ",
        ct
      )
    }
    
    # ----------------------------------------------------------
    # Create a group variable for LPEseq2 design
    # ----------------------------------------------------------
    meta_ct$group <- factor(meta_ct$condition)
    
    group_table <- table(meta_ct$group)
    
    # ----------------------------------------------------------
    # Store diagnostic information for this cell type
    # ----------------------------------------------------------
    diagnostic_list[[ct]] <- list(
      celltype = ct,
      n_pseudobulk_samples = ncol(counts_ct),
      n_genes_before_filter = nrow(counts_ct),
      group_table = group_table,
      n_cells = meta_ct$n_cells,
      library_size = meta_ct$library_size
    )
    
    # ----------------------------------------------------------
    # Skip analysis if fewer than two groups are available
    # ----------------------------------------------------------
    if (nlevels(meta_ct$group) < 2) {
      error_list[[ct]] <- "Skipped: fewer than two groups."
      next
    }
    
    # ----------------------------------------------------------
    # Warn when any group has only one pseudo-bulk sample
    # ----------------------------------------------------------
    if (any(group_table == 1)) {
      warning(
        "Cell type '", ct, "' has at least one group with only one pseudo-bulk sample. ",
        "Statistical inference may be unstable."
      )
    }
    
    # ----------------------------------------------------------
    # Remove pseudo-bulk samples with zero library size
    # ----------------------------------------------------------
    lib_size <- colSums(counts_ct)
    
    if (any(lib_size == 0)) {
      keep_samples <- lib_size > 0
      
      counts_ct <- counts_ct[, keep_samples, drop = FALSE]
      meta_ct <- meta_ct[keep_samples, , drop = FALSE]
      meta_ct$group <- droplevels(meta_ct$group)
    }
    
    # ----------------------------------------------------------
    # Re-check whether sufficient samples and groups remain
    # ----------------------------------------------------------
    if (ncol(counts_ct) < 2 || nlevels(meta_ct$group) < 2) {
      error_list[[ct]] <- "Skipped: insufficient pseudo-bulk samples after zero-library filtering."
      next
    }
    
    # ----------------------------------------------------------
    # 5. Run LPEseq2 preprocessing and LPE-ANOVA
    # ----------------------------------------------------------
    run_one <- tryCatch(
      {
        # Preprocess raw counts:
        # filtering, normalization, and log transformation are handled here
        prep <- LPE_preprocess(
          counts = counts_ct,
          colData = meta_ct,
          design = ~ group,
          normalize.method = normalize.method,
          log.transform = log.transform,
          min.count = min.count,
          prior.count = prior.count,
          verbose = verbose
        )
        
        # Run LPEseq2 differential expression analysis
        res <- LPE_ANOVA(
          object = prep,
          analysis.method = analysis.method,
          standard.min.group.n = standard.min.group.n,
          n.bin = n.bin,
          df = df,
          trim.method = trim.method,
          use_weighted_between = use_weighted_between,
          p.method = p.method,
          verbose = verbose
        )
        
        # Add cell type and pseudo-bulk method information
        res$celltype <- ct
        res$requested_pseudobulk_method <- method
        
        selected_method <- attr(res, "analysis.method")
        if (is.null(selected_method)) {
          selected_method <- NA
        }
        
        res$selected_LPEseq2_method <- selected_method
        
        # Move important annotation columns to the front
        front_cols <- intersect(
          c("celltype", "gene", "requested_pseudobulk_method", "selected_LPEseq2_method"),
          colnames(res)
        )
        
        other_cols <- setdiff(colnames(res), front_cols)
        
        res <- res[, c(front_cols, other_cols), drop = FALSE]
        
        list(
          prep = prep,
          result = res
        )
      },
      error = function(e) {
        e
      }
    )
    
    # ----------------------------------------------------------
    # Store errors if the analysis failed for this cell type
    # ----------------------------------------------------------
    if (inherits(run_one, "error")) {
      error_list[[ct]] <- conditionMessage(run_one)
      next
    }
    
    # ----------------------------------------------------------
    # Store successful results
    # ----------------------------------------------------------
    prep_list[[ct]] <- run_one$prep
    result_list[[ct]] <- run_one$result
    count_list[[ct]] <- counts_ct
    metadata_list[[ct]] <- meta_ct
  }
  
  # ------------------------------------------------------------
  # 6. Combine cell-type-specific DEG results into one table
  # ------------------------------------------------------------
  combined_result <- NULL
  
  if (length(result_list) > 0) {
    combined_result <- do.call(rbind, result_list)
    rownames(combined_result) <- NULL
  }
  
  # ------------------------------------------------------------
  # 7. Return all outputs
  # ------------------------------------------------------------
  return(
    list(
      pseudobulk = pb,
      counts_by_celltype = count_list,
      metadata_by_celltype = metadata_list,
      prep_by_celltype = prep_list,
      results_by_celltype = result_list,
      combined_result = combined_result,
      diagnostics = diagnostic_list,
      errors = error_list,
      settings = list(
        pseudobulk_method = method,
        normalize.method = normalize.method,
        log.transform = log.transform,
        min.count = min.count,
        prior.count = prior.count,
        analysis.method = analysis.method,
        standard.min.group.n = standard.min.group.n,
        n.bin = n.bin,
        df = df,
        trim.method = trim.method,
        use_weighted_between = use_weighted_between,
        p.method = p.method
      )
    )
  )
}