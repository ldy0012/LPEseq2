#' Run LPE-ANOVA
#'
#' Performs multi-group differential expression analysis using an
#' intensity-dependent local pooled error variance trend.
#'
#' @param object A list returned by \code{LPE_preprocess()}.
#' @param n.bin Number of quantile bins used for variance trend estimation.
#' @param df Degrees of freedom for smoothing spline.
#' @param trim.method Outlier trimming method applied to pairwise values
#'   used for variance trend estimation. One of \code{"iqr"} or \code{"none"}.
#'   \code{"iqr"} applies the conventional 1.5*IQR boxplot rule within
#'   expression-intensity A-bins after pooling within-group and between-group
#'   values. \code{"none"} performs no outlier trimming.
#' @param trend.method Method used to fit the variance trend. One of
#'   \code{"mean_spline"} or \code{"quantile_regression"}.
#'   \code{"mean_spline"} fits a smoothing spline to bin-level variance
#'   estimates. \code{"quantile_regression"} fits an upper-quantile
#'   regression trend to reduce potential variance underestimation.
#' @param tau Quantile level used when
#' @param use_weighted_between Logical. Whether to include weighted between-group
#'   differences in variance trend estimation.
#' @param analysis.method Analysis method. One of \code{"LPE"},
#'   \code{"standard_anova"}, or \code{"auto"}. \code{"LPE"} uses the
#'   local pooled error-based ANOVA. \code{"standard_anova"} uses conventional
#'   gene-wise one-way ANOVA. \code{"auto"} selects standard ANOVA when all
#'   groups have at least \code{standard.min.group.n} samples; otherwise,
#'   LPE-ANOVA is used.
#' @param standard.min.group.n Minimum per-group sample size required to use
#'   standard one-way ANOVA when \code{analysis.method = "auto"}.
#' @param verbose Logical. Whether to print progress messages.
#' @param p.method P-value calculation method. One of \code{"chisq"} or \code{"F_inf"}.
#'
#' @return A data.frame containing gene-level test statistics. The selected
#'   analysis method is stored in \code{attr(result, "analysis.method")}.
#'   For LPE-ANOVA, trimming information, variance trend information, and
#'   bin-level variance points are stored in \code{attr(result, "trim.info")},
#'   \code{attr(result, "trend.info")}, and \code{attr(result, "base.var")}.
#'
#' @export
LPE_ANOVA <- function(object,
                      n.bin = 100,
                      df = 10,
                      trim.method = c("iqr", "none"),
                      trend.method = c("mean_spline", "quantile_regression"),
                      tau = 0.75,
                      use_weighted_between = FALSE,
                      analysis.method = c("LPE", "standard_anova", "auto"),
                      standard.min.group.n = 5,
                      verbose = TRUE,
                      p.method = c("chisq", "F_inf")) {

  trim.method <- match.arg(trim.method)
  trend.method <- match.arg(trend.method)
  analysis.method <- match.arg(analysis.method)
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

  if (!is.numeric(standard.min.group.n) ||
      length(standard.min.group.n) != 1 ||
      standard.min.group.n < 2) {
    stop("standard.min.group.n must be a single numeric value >= 2")
  }

  if (!is.numeric(tau) || length(tau) != 1 || tau <= 0 || tau >= 1) {
    stop("tau must be a single numeric value between 0 and 1")
  }

  k <- nlevels(group)
  n_i <- table(group)

  if (k < 2) {
    stop("At least two groups are required")
  }

  selected.method <- analysis.method

  if (analysis.method == "auto") {
    if (min(n_i) >= standard.min.group.n) {
      selected.method <- "standard_anova"
    } else {
      selected.method <- "LPE"
    }
  }

  if (verbose) {
    cat("Selected analysis method:", selected.method, "\n")
    print(n_i)
  }

  if (selected.method == "standard_anova") {
    if (verbose) {
      cat("Running standard one-way ANOVA\n")
    }

    res <- standard_ANOVA_expr(
      expr = expr,
      group = group,
      p.adjust.method = "BH"
    )

    attr(res, "analysis.method") <- "standard_anova"
    attr(res, "requested.analysis.method") <- analysis.method
    attr(res, "standard.min.group.n") <- standard.min.group.n
    attr(res, "trim.info") <- NULL
    attr(res, "trend.info") <- NULL

    return(res)
  }

  if (verbose) {
    cat("Running LPE-ANOVA\n")
  }

  if (trend.method == "quantile_regression" &&
      !requireNamespace("quantreg", quietly = TRUE)) {
    stop(
      "Package 'quantreg' is required when trend.method = 'quantile_regression'. ",
      "Install it with install.packages('quantreg')."
    )
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
    trend.method = trend.method,
    tau = tau,
    use_weighted_between = use_weighted_between
  )

  trim.info <- attr(var.spline, "trim.info")
  trend.info <- attr(var.spline, "trend.info")
  base.var <- attr(var.spline, "base.var")

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
    method = "LPE",
    row.names = NULL
  )

  attr(res, "analysis.method") <- "LPE"
  attr(res, "requested.analysis.method") <- analysis.method
  attr(res, "standard.min.group.n") <- standard.min.group.n
  attr(res, "trim.info") <- trim.info
  attr(res, "trend.info") <- trend.info
  attr(res, "base.var") <- base.var
  attr(res, "base.var") <- NULL

  return(res)
}
