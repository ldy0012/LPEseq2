#' Preprocess RNA-seq Count Data for LPE-ANOVA
#'
#' This function checks raw count data, filters low-count genes, normalizes
#' sample-level library size differences, and optionally applies log2
#' transformation.
#'
#' @param counts A numeric matrix of raw counts with genes as rows and samples as columns.
#' @param colData A data.frame containing sample-level metadata. Row names must match column names of counts.
#' @param design A design formula with no left-hand side, e.g. \code{~ group}
#'   or \code{~ genotype + treatment}. One or more main-effect terms are
#'   supported (multi-way ANOVA without interaction); interaction terms
#'   (\code{A:B} or \code{A*B}) are not supported and will raise an error.
#'   Every term must name a column of \code{colData} that can be coerced to
#'   a factor with at least two levels.
#' @param normalize.method Normalization method. One of \code{"TMM"}, \code{"library_size"}, \code{"DESeq2"}, or \code{"none"}.
#' @param log.transform Logical. Whether to apply log2 transformation.
#' @param min.count Minimum count threshold for low-count filtering.
#' @param prior.count Prior count added before log transformation.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @return A list containing normalized expression matrix, a factor
#'   identifying each unique combination of the design factors (\code{group},
#'   kept for backward compatibility with one-way designs), a named list of
#'   the individual design factors (\code{factors}), design, colData, and
#'   preprocessing options.
#'
#' @export
LPE_preprocess <- function(counts,
                           colData,
                           design = ~ group,
                           normalize.method = c("library_size", "TMM", "DESeq2", "none"),
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
  # Multi-way (main-effects-only) designs are supported: ~ A, ~ A + B,
  # ~ A + B + C, ... Interaction terms (A:B, A*B) are rejected because the
  # LPE-ANOVA variance-trend / test-statistic machinery below assumes an
  # additive model.

  if (!inherits(design, "formula")) {
    stop("design must be a formula, e.g. ~ group or ~ A + B")
  }

  design_tt <- stats::terms(design)

  if (!is.null(attr(design_tt, "response")) && attr(design_tt, "response") != 0) {
    stop("design must not have a left-hand side, e.g. use ~ group, not y ~ group")
  }

  if (attr(design_tt, "intercept") == 0) {
    stop("design must include an intercept (do not use ~ . - 1 or ~ . + 0)")
  }

  term_labels <- attr(design_tt, "term.labels")
  term_order  <- attr(design_tt, "order")

  if (length(term_labels) < 1) {
    stop("design must contain at least one factor term, e.g. ~ group or ~ A + B")
  }

  if (any(term_order > 1)) {
    stop(
      "Interaction terms are not supported (e.g. 'A:B' or 'A*B'). ",
      "LPE_ANOVA()/LPE_preprocess() only support main-effect (additive) ",
      "multi-way designs, e.g. ~ A + B."
    )
  }

  # terms() backquotes non-syntactic names (e.g. a column called "my group"
  # comes back as "`my group`"); strip them so they match colData columns.
  design_terms <- gsub("^`|`$", "", term_labels)

  missing_terms <- setdiff(design_terms, colnames(colData))
  if (length(missing_terms) > 0) {
    stop(
      "The following design variable(s) were not found in colData: ",
      paste(missing_terms, collapse = ", ")
    )
  }

  factors <- vector("list", length(design_terms))
  names(factors) <- design_terms

  for (v in design_terms) {
    if (anyNA(colData[[v]])) {
      stop("Design variable '", v, "' contains NA")
    }

    f <- droplevels(as.factor(colData[[v]]))

    if (nlevels(f) < 2) {
      stop("Design variable '", v, "' must have at least two levels")
    }

    factors[[v]] <- f
  }

  # `group`: a single factor identifying each unique combination of the
  # design factors (i.e. each experimental "cell"). For a one-way design
  # (length(design_terms) == 1) this is exactly the original factor, so
  # existing one-way callers of LPE_ANOVA()/LPE_ANOVA_var() that only look
  # at object$group are unaffected. For a multi-way design, it is the
  # crossing of all factors, and is what LPE_ANOVA_var() uses to define
  # "true replicates" (samples sharing the same combination of factor
  # levels) when estimating the intensity-dependent variance trend.
  if (length(design_terms) == 1) {
    group <- factors[[1]]
  } else {
    group <- droplevels(interaction(factors, drop = TRUE, sep = "."))
  }

  if (any(table(group) == 1)) {
    warning(
      if (length(design_terms) == 1) {
        "Some groups have only one sample. Inference may be unstable."
      } else {
        "Some factor-level combinations have only one sample. Inference may be unstable."
      }
    )
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

  filtered_lib.size <- colSums(counts)

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

    size.factor <- filtered_lib.size / stats::median(filtered_lib.size)
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
    if (length(design_terms) == 1) {
      cat("Groups detected:\n")
      print(table(group))
    } else {
      cat("Design factors (no interaction):", paste(design_terms, collapse = " + "), "\n")
      for (v in design_terms) {
        cat(" -", v, ":\n")
        print(table(factors[[v]]))
      }
      cat("Factor-level combinations (cells) detected:\n")
      print(table(group))
    }
    cat("Genes retained:", nrow(expr), "\n")
    cat("Normalization method:", normalize.method, "\n")
    cat("Log transform:", log.transform, "\n")
  }

  return(list(
    expr = expr,
    group = group,
    factors = factors,
    design = design,
    colData = colData,
    normalize.method = normalize.method,
    log.transform = log.transform,
    min.count = min.count,
    prior.count = prior.count
  ))
}
