#' Run LPE-ANOVA
#'
#' Performs multi-group differential expression analysis using an
#' intensity-dependent local pooled error variance trend.
#'
#' @param object A list returned by \code{LPE_preprocess()}.
#' @param n.bin Number of quantile bins used for variance trend estimation.
#' @param df Degrees of freedom for smoothing spline.
#' @param trim.method Outlier trimming method. One of \code{"mad"}, \code{"quantile"}, or \code{"fixed"}.
#' @param d Fixed trimming threshold used when \code{trim.method = "fixed"}.
#' @param use_weighted_between Logical. Whether to include weighted between-group differences in variance trend estimation.
#' @param verbose Logical. Whether to print progress messages.
#' @param p.method P-value calculation method. One of \code{"chisq"} or \code{"F_inf"}.
#'
#' @return A data.frame containing gene-level LPE-ANOVA statistics.
#'
#' @export
LPE_ANOVA <- function(object,
                      n.bin = 100,
                      df = 10,
                      trim.method = c("mad", "quantile", "fixed"),
                      d = 1.2,
                      use_weighted_between = FALSE,
                      verbose = TRUE,
                      p.method = c("chisq", "F_inf")) {

  trim.method <- match.arg(trim.method)
  p.method <- match.arg(p.method)

  # -----------------------------
  # 1. input check
  # -----------------------------

  if (!is.list(object)) {
    stop("object must be a list returned by LPE_preprocess()")
  }

  if (is.null(object$expr)) {
    stop("object must contain expr")
  }

  if (is.null(object$group)) {
    stop("object must contain group")
  }

  expr <- as.matrix(object$expr)
  group <- droplevels(as.factor(object$group))

  if (!is.numeric(expr)) {
    stop("object$expr must be numeric")
  }

  if (anyNA(expr) || any(!is.finite(expr))) {
    stop("object$expr contains NA, NaN, or infinite values")
  }

  if (ncol(expr) != length(group)) {
    stop("The number of columns in object$expr must match the length of object$group")
  }

  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("verbose must be TRUE or FALSE")
  }

  if (!is.logical(use_weighted_between) || length(use_weighted_between) != 1) {
    stop("use_weighted_between must be TRUE or FALSE")
  }

  k <- nlevels(group)
  n_i <- table(group)

  if (k < 2) {
    stop("At least two groups are required")
  }

  if (verbose) {
    cat("Running LPE-ANOVA\n")
    print(n_i)
  }

  # -----------------------------
  # 2. variance trend estimation
  # -----------------------------

  var.spline <- LPE_ANOVA_var(
    expr = expr,
    group = group,
    n.bin = n.bin,
    df = df,
    trim.method = trim.method,
    d = d,
    use_weighted_between = use_weighted_between
  )

  gene.mean <- rowMeans(expr, na.rm = TRUE)

  pred.var <- fixbounds.predict.smooth.spline(
    var.spline,
    gene.mean
  )$y

  positive_y <- var.spline$y[is.finite(var.spline$y) & var.spline$y > 0]

  if (length(positive_y) == 0) {
    stop("Variance spline produced no positive fitted values")
  }

  var_floor <- min(positive_y, na.rm = TRUE)

  pred.var[!is.finite(pred.var)] <- var_floor
  pred.var <- pmax(pred.var, var_floor)

  # -----------------------------
  # 3. between-group mean square
  # -----------------------------

  MS_between <- apply(expr, 1, function(x) {
    group_means <- tapply(x, group, mean)
    grand_mean <- mean(x)

    sum(as.numeric(n_i) * (group_means - grand_mean)^2) / (k - 1)
  })

  Fstat <- MS_between / pred.var

  # -----------------------------
  # 4. p-value calculation
  # -----------------------------

  if (p.method == "chisq") {
    p.val <- stats::pchisq(
      (k - 1) * Fstat,
      df = k - 1,
      lower.tail = FALSE
    )
  } else {
    p.val <- stats::pf(
      Fstat,
      df1 = k - 1,
      df2 = 1e6,
      lower.tail = FALSE
    )
  }

  p.val[!is.finite(p.val)] <- 1
  p.val <- pmax(p.val, .Machine$double.xmin)

  adj.p <- stats::p.adjust(p.val, method = "BH")

  # -----------------------------
  # 5. output
  # -----------------------------

  if (verbose) {
    cat("Finished.\n")
  }

  gene_id <- rownames(expr)

  if (is.null(gene_id)) {
    gene_id <- paste0("gene_", seq_len(nrow(expr)))
  }

  res <- data.frame(
    gene = gene_id,
    mean = gene.mean,
    var = pred.var,
    MS_between = MS_between,
    F = Fstat,
    p.value = p.val,
    q.value = adj.p,
    row.names = NULL
  )

  return(res)
}
