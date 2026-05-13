#' Preprocess RNA-seq Count Data for LPE-ANOVA
#'
#' This function checks raw count data, filters low-count genes, normalizes
#' sample-level library size differences, and optionally applies log2
#' transformation.
#'
#' @param counts A numeric matrix of raw counts with genes as rows and samples as columns.
#' @param colData A data.frame containing sample-level metadata. Row names must match column names of counts.
#' @param design A one-way design formula, such as \code{~ group}.
#' @param normalize.method Normalization method. One of \code{"TMM"}, \code{"library_size"}, \code{"DESeq2"}, or \code{"none"}.
#' @param log.transform Logical. Whether to apply log2 transformation.
#' @param min.count Minimum count threshold for low-count filtering.
#' @param prior.count Prior count added before log transformation.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @return A list containing normalized expression matrix, group factor, design, colData, and preprocessing options.
#'
#' @export
LPE_preprocess <- function(counts,
                           colData,
                           design = ~ group,
                           normalize.method = c("TMM", "library_size", "DESeq2", "none"),
                           log.transform = TRUE,
                           min.count = 5,
                           prior.count = 1,
                           verbose = TRUE) {

  normalize.method <- match.arg(normalize.method)

  # -----------------------------
  # 1. input check
  # -----------------------------

  if (!is.matrix(counts)) {
    stop("counts must be a matrix with genes as rows and samples as columns")
  }

  if (!is.numeric(counts)) {
    stop("counts must be numeric")
  }

  if (anyNA(counts)) {
    stop("counts contains NA values")
  }

  if (any(!is.finite(counts))) {
    stop("counts contains NaN or infinite values")
  }

  if (any(counts < 0)) {
    stop("counts must be non-negative")
  }

  if (is.null(colnames(counts))) {
    stop("counts must have sample names as column names")
  }

  if (is.null(rownames(counts))) {
    rownames(counts) <- paste0("gene_", seq_len(nrow(counts)))
  }

  if (!is.data.frame(colData)) {
    stop("colData must be a data.frame")
  }

  if (is.null(rownames(colData))) {
    stop("colData must have sample names as rownames")
  }

  if (anyDuplicated(colnames(counts))) {
    stop("Duplicated sample names in counts")
  }

  if (anyDuplicated(rownames(colData))) {
    stop("Duplicated sample names in colData")
  }

  if (!setequal(colnames(counts), rownames(colData))) {
    stop("Sample names of counts and rownames of colData must match exactly")
  }

  if (!is.numeric(min.count) || length(min.count) != 1 || min.count < 0) {
    stop("min.count must be a single non-negative numeric value")
  }

  if (!is.numeric(prior.count) || length(prior.count) != 1 || prior.count < 0) {
    stop("prior.count must be a single non-negative numeric value")
  }

  if (!is.logical(log.transform) || length(log.transform) != 1) {
    stop("log.transform must be TRUE or FALSE")
  }

  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("verbose must be TRUE or FALSE")
  }

  colData <- colData[colnames(counts), , drop = FALSE]

  # -----------------------------
  # 2. design check
  # -----------------------------

  design_terms <- all.vars(design)

  if (length(design_terms) != 1) {
    stop("Only one-way design is currently supported")
  }

  group_var <- design_terms[1]

  if (!group_var %in% colnames(colData)) {
    stop("Group variable not found in colData")
  }

  if (anyNA(colData[[group_var]])) {
    stop("Group variable contains NA")
  }

  group <- droplevels(as.factor(colData[[group_var]]))

  if (nlevels(group) < 2) {
    stop("At least two groups are required")
  }

  if (any(table(group) == 1)) {
    warning("Some groups have only one sample. Inference may be unstable.")
  }

  # -----------------------------
  # 3. raw library size
  # -----------------------------

  raw_lib.size <- colSums(counts)

  if (any(raw_lib.size == 0)) {
    stop("At least one sample has zero total counts")
  }

  # -----------------------------
  # 4. low-count filtering
  # -----------------------------

  keep <- rowSums(counts >= min.count) >= 2
  counts <- counts[keep, , drop = FALSE]

  if (nrow(counts) < 10) {
    stop("Too few genes retained after filtering")
  }

  # -----------------------------
  # 5. normalization + transformation
  # -----------------------------

  if (normalize.method %in% c("TMM", "DESeq2")) {
    if (any(counts != round(counts))) {
      warning(
        "counts contains non-integer values. ",
        "Raw integer RNA-seq counts are recommended for TMM and DESeq2 normalization."
      )
    }
  }

  if (normalize.method == "none") {

    expr <- counts

    if (log.transform) {
      expr <- log2(expr + prior.count)
    }

  } else if (normalize.method == "library_size") {

    size.factor <- raw_lib.size / stats::median(raw_lib.size)
    expr <- sweep(counts, 2, size.factor, "/")

    if (log.transform) {
      expr <- log2(expr + prior.count)
    }

  } else if (normalize.method == "TMM") {

    if (!requireNamespace("edgeR", quietly = TRUE)) {
      stop("Package 'edgeR' is required for normalize.method = 'TMM'")
    }

    dge <- edgeR::DGEList(counts = counts)
    dge <- edgeR::calcNormFactors(dge, method = "TMM")

    if (log.transform) {
      expr <- edgeR::cpm(dge, log = TRUE, prior.count = prior.count)
    } else {
      expr <- edgeR::cpm(dge, log = FALSE)
    }

  } else if (normalize.method == "DESeq2") {

    if (!requireNamespace("DESeq2", quietly = TRUE)) {
      stop("Package 'DESeq2' is required for normalize.method = 'DESeq2'")
    }

    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = round(counts),
      colData = colData,
      design = design
    )

    dds <- DESeq2::estimateSizeFactors(dds)
    norm_counts <- DESeq2::counts(dds, normalized = TRUE)

    if (log.transform) {
      expr <- log2(norm_counts + prior.count)
    } else {
      expr <- norm_counts
    }
  }

  expr <- as.matrix(expr)

  # -----------------------------
  # 6. output
  # -----------------------------

  if (verbose) {
    cat("Groups detected:\n")
    print(table(group))
    cat("Genes retained:", nrow(expr), "\n")
    cat("Normalization method:", normalize.method, "\n")
    cat("Log transform:", log.transform, "\n")
  }

  return(list(
    expr = expr,
    group = group,
    design = design,
    colData = colData,
    normalize.method = normalize.method,
    log.transform = log.transform,
    min.count = min.count,
    prior.count = prior.count
  ))
}
