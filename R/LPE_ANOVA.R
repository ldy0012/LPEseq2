#' Run LPE-ANOVA
#'
#' Performs multi-group differential expression analysis using an
#' intensity-dependent local pooled error variance trend.
#'
#' @param object A list returned by \code{LPE_preprocess()}.
#' @param n.bin Number of quantile bins used for variance trend estimation.
#' @param df Degrees of freedom for smoothing spline.
#' @param trim.method Outlier trimming method applied to pairwise values
#'   used for variance trend estimation. One of \code{"iqr"}, \code{"dvalue"},
#'   or \code{"none"}. See \code{\link{LPE_ANOVA_var}} for details.
#' @param use_weighted_between Logical. Whether to include weighted between-group
#'   differences in variance trend estimation.
#' @param d.threshold Numeric. Fixed threshold applied to the raw pairwise
#'   difference D when \code{trim.method = "dvalue"}. Default 1.2, matching
#'   the default reported in LPEseq1 (Gim et al. 2016). Ignored for other
#'   \code{trim.method} values.
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
#' @param variance.eval Character string specifying how the fitted variance
#'   trend is evaluated for the test statistic. One of \code{"grand_mean"}
#'   or \code{"per_group"}. \code{"grand_mean"} (default) evaluates the
#'   variance trend once at each gene's overall (grand) mean expression and
#'   uses it as a single pooled variance in a standard one-way ANOVA-style
#'   F statistic. \code{"per_group"} evaluates the variance trend
#'   separately at each group's mean expression and combines the
#'   group-specific variances via inverse-variance weighting into a
#'   Welch-type unequal-variance test statistic (returned with
#'   \code{method = "LPE_Welch"} and a different \code{var}/\code{MS_between}
#'   definition; see Details).
#' @return A data.frame containing gene-level test statistics. The selected
#'   analysis method is stored in \code{attr(result, "analysis.method")}.
#'   analysis method is stored in \code{attr(result, "analysis.method")}.
#'   For LPE-ANOVA, trimming information, variance trend information,
#'   bin-level variance points, and the fitted variance trend object are
#'   stored in \code{attr(result, "trim.info")},
#'   \code{attr(result, "trend.info")}, \code{attr(result, "base.var")},
#'   and \code{attr(result, "var.spline")}.
#'
#' @export
LPE_ANOVA <- function(object, n.bin = 100, df = 10,
                      trim.method = c("iqr", "dvalue", "none"),
                      use_weighted_between = FALSE, d.threshold = 1.2,
                      analysis.method = c("LPE", "standard_anova", "auto"),
                      standard.min.group.n = 5, verbose = TRUE,
                      p.method = c("chisq", "F_inf"),
                      variance.eval = c("grand_mean", "per_group")) {

  variance.eval <- match.arg(variance.eval)
  trim.method <- match.arg(trim.method)
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
    attr(res, "base.var") <- NULL
    attr(res, "var.spline") <- NULL

    return(res)
  }

  if (verbose) {
    cat("Running LPE-ANOVA\n")
  }

  # -----------------------------
  # 2. variance trend estimation
  # -----------------------------

  var.result <- LPE_ANOVA_var(
    expr = expr,
    group = group,
    n.bin = n.bin,
    df = df,
    trim.method = trim.method,
    use_weighted_between = use_weighted_between,
    d.threshold = d.threshold
  )

  trim.info <- attr(var.result, "trim.info")
  trend.info <- attr(var.result, "trend.info")
  base.var <- attr(var.result, "base.var")

  gene.mean <- rowMeans(expr, na.rm = TRUE)

  if (variance.eval == "grand_mean") {
    pred.var <- fixbounds.predict.smooth.spline(var.result$object, gene.mean)$y
  } else {
    group.levels <- levels(group)
    group.means <- sapply(group.levels, function(g) {
      rowMeans(expr[, group == g, drop = FALSE], na.rm = TRUE)
    })
    group.vars <- sapply(seq_len(k), function(i) {
      fixbounds.predict.smooth.spline(var.result$object, group.means[, i])$y
    })
  }

  positive_y <- var.result$object$y[is.finite(var.result$object$y) & var.result$object$y > 0]
  var_floor <- min(positive_y, na.rm = TRUE)

  if (variance.eval == "grand_mean") {
    pred.var[!is.finite(pred.var)] <- var_floor
    pred.var <- pmax(pred.var, var_floor)
  } else {
    group.vars[!is.finite(group.vars)] <- var_floor
    group.vars <- pmax(group.vars, var_floor)
  }

  # -----------------------------
  # 3. between-group mean square
  # -----------------------------

  if (variance.eval == "grand_mean") {
    MS_between <- apply(expr, 1, function(x) {
      group_means <- tapply(x, group, mean)
      grand_mean <- mean(x)
      sum(as.numeric(n_i) * (group_means - grand_mean)^2) / (k - 1)
    })
    Fstat <- MS_between / pred.var
  } else {
    n_i_vec <- as.numeric(n_i)
    w <- sweep(1 / group.vars, 2, n_i_vec, "*")
    weighted.mean <- rowSums(w * group.means) / rowSums(w)
    T.stat <- rowSums(w * (group.means - weighted.mean)^2)
  }

  # -----------------------------
  # 4. p-value calculation
  # -----------------------------

  if (variance.eval == "grand_mean") {
    if (p.method == "chisq") {
      p.val <- stats::pchisq((k - 1) * Fstat, df = k - 1, lower.tail = FALSE)
    } else {
      p.val <- stats::pf(Fstat, df1 = k - 1, df2 = 1e6, lower.tail = FALSE)
    }
  } else {
    p.val <- stats::pchisq(T.stat, df = k - 1, lower.tail = FALSE)
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

  if (variance.eval == "grand_mean") {
    res <- data.frame(
      gene = gene_id, mean = gene.mean, var = pred.var,
      MS_between = MS_between, F = Fstat,
      p.value = p.val, q.value = adj.p, method = "LPE", row.names = NULL
    )
  } else {
    res <- data.frame(
      gene = gene_id, mean = gene.mean,
      var = 1 / rowSums(w),
      MS_between = T.stat,
      F = T.stat / (k - 1),
      p.value = p.val, q.value = adj.p, method = "LPE_Welch", row.names = NULL
    )
  }

  attr(res, "analysis.method") <- "LPE"
  attr(res, "requested.analysis.method") <- analysis.method
  attr(res, "standard.min.group.n") <- standard.min.group.n
  attr(res, "trim.info") <- trim.info
  attr(res, "trend.info") <- trend.info
  attr(res, "base.var") <- base.var
  attr(res, "var.spline") <- var.result

  return(res)
}
