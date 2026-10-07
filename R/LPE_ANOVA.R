#' Run LPE-ANOVA
#'
#' Performs multi-group (and multi-way, interaction-free) differential
#' expression analysis using an intensity-dependent local pooled error
#' variance trend.
#'
#' @param object A list returned by \code{LPE_preprocess()}.
#' @param n.bin Number of quantile bins used for variance trend estimation.
#' @param df Degrees of freedom for smoothing spline.
#' @param trim.method Outlier trimming method applied to pairwise values
#'   used for variance trend estimation. One of \code{"iqr"}, \code{"dvalue"},
#'   or \code{"none"}. See \code{\link{LPE_ANOVA_var}} for details.
#' @param use_weighted_between Logical. Whether to include weighted between-group
#'   differences in variance trend estimation.
#' @param d.threshold Numeric. Threshold on the LPEseq1 raw-D scale, used to
#'   derive the fixed M-scale cutoff (\code{d.threshold/sqrt(2)}) applied
#'   when \code{trim.method = "dvalue"}. Default 1.2, matching the default
#'   reported in LPEseq1 (Gim et al. 2016). Ignored for other
#'   \code{trim.method} values. See \code{\link{LPE_ANOVA_var}} for details.
#' @param analysis.method Analysis method. One of \code{"LPE"},
#'   \code{"standard_anova"}, or \code{"auto"}. \code{"LPE"} uses the
#'   local pooled error-based ANOVA. \code{"standard_anova"} uses conventional
#'   gene-wise ANOVA (Type II sums of squares for multi-way designs).
#'   \code{"auto"} selects standard ANOVA when every factor-level
#'   combination (design "cell") has at least \code{standard.min.group.n}
#'   samples; otherwise, LPE-ANOVA is used.
#' @param standard.min.group.n Minimum per-cell sample size required to use
#'   standard ANOVA when \code{analysis.method = "auto"}. A "cell" is a
#'   unique combination of the design factor levels (for a one-way design,
#'   this is simply a group).
#' @param verbose Logical. Whether to print progress messages.
#' @param p.method P-value calculation method. One of \code{"chisq"} or \code{"F_inf"}.
#' @param variance.eval Character string specifying how the fitted variance
#'   trend is evaluated for the test statistic. One of \code{"grand_mean"}
#'   or \code{"per_group"}. \code{"grand_mean"} (default) evaluates the
#'   variance trend once at each gene's overall (grand) mean expression and
#'   uses it as a single pooled variance for the test statistic; this mode
#'   supports both one-way and multi-way (main-effects-only) designs.
#'   \code{"per_group"} evaluates the variance trend separately at each
#'   group's mean expression and combines the group-specific variances via
#'   inverse-variance weighting into a Welch-type unequal-variance test
#'   statistic (returned with \code{method = "LPE_Welch"} and a different
#'   \code{var}/\code{MS_between} definition; see Details). \code{"per_group"}
#'   is only supported for one-way designs (a single design term); it raises
#'   an error for multi-way designs.
#' @return A data.frame containing gene-level test statistics. For a
#'   one-way design, columns are unsuffixed (\code{MS_between}, \code{F},
#'   \code{p.value}, \code{q.value}, ...), matching previous versions. For a
#'   multi-way design (e.g. \code{design = ~ A + B}), one main-effect test is
#'   reported per design term, with columns suffixed by the term name (e.g.
#'   \code{F_A}, \code{p.value_A}, \code{q.value_A}, \code{F_B},
#'   \code{p.value_B}, \code{q.value_B}, ...); \code{gene}, \code{mean},
#'   \code{var}, and \code{method} are shared across terms. The selected
#'   analysis method is stored in \code{attr(result, "analysis.method")},
#'   and the design term names are stored in
#'   \code{attr(result, "design.terms")}. For LPE-ANOVA, trimming
#'   information, variance trend information, bin-level variance points, and
#'   the fitted variance trend object are stored in
#'   \code{attr(result, "trim.info")}, \code{attr(result, "trend.info")},
#'   \code{attr(result, "base.var")}, and \code{attr(result, "var.spline")}.
#'
#' @details
#' For a multi-way design \code{~ A + B + ...} (no interaction terms
#' allowed), each design term's main effect is tested by comparing, per
#' gene, the full additive model against the model with that term dropped:
#' \code{SS_term = RSS(reduced) - RSS(full)}, with
#' \code{df_term = rank(X_full) - rank(X_reduced)}. This is the standard
#' Type II sum-of-squares decomposition, which for a main-effects-only model
#' does not depend on term order. The intensity-dependent pooled variance
#' (from the same LPE variance-trend spline used for one-way designs) is
#' used as the common denominator for every term's test statistic. This
#' reduces exactly to the original one-way LPE-ANOVA when the design has a
#' single term.
#'
#' @export
LPE_ANOVA <- function(object, n.bin = 100, df = 10,
                      trim.method = c("dvalue", "iqr", "none"),
                      use_weighted_between = TRUE, d.threshold = 1.2,
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

  # `factors`: a named list of the individual design factors. LPE_preprocess()
  # always sets this; a manually-constructed object (only expr + group) is
  # still supported by falling back to a one-way design named "group".
  if (is.null(object$factors)) {
    factors <- stats::setNames(list(group), "group")
  } else {
    factors <- lapply(object$factors, function(f) droplevels(as.factor(f)))
  }

  design.terms <- names(factors)

  if (!is.numeric(expr)) {
    stop("object$expr must be numeric")
  }

  if (anyNA(expr) || any(!is.finite(expr))) {
    stop("object$expr contains NA, NaN, or infinite values")
  }

  if (ncol(expr) != length(group)) {
    stop("The number of columns in object$expr must match the length of object$group")
  }

  for (nm in design.terms) {
    if (length(factors[[nm]]) != ncol(expr)) {
      stop("The number of columns in object$expr must match the length of each design factor")
    }
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

  if (variance.eval == "per_group" && length(design.terms) > 1) {
    stop(
      "variance.eval = 'per_group' (Welch-type) is only supported for ",
      "one-way designs (a single design term). Use variance.eval = ",
      "'grand_mean' for multi-way designs (", paste(design.terms, collapse = " + "), ")."
    )
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
    if (length(design.terms) > 1) {
      cat("Design terms (no interaction):", paste(design.terms, collapse = " + "), "\n")
    }
    print(n_i)
  }

  if (selected.method == "standard_anova") {
    if (verbose) {
      cat("Running standard ANOVA\n")
    }

    res <- standard_ANOVA_expr(
      expr = expr,
      factors = factors,
      p.adjust.method = "BH"
    )

    attr(res, "analysis.method") <- "standard_anova"
    attr(res, "requested.analysis.method") <- analysis.method
    attr(res, "standard.min.group.n") <- standard.min.group.n
    attr(res, "design.terms") <- design.terms
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
  # Uses `group` (the crossing of all design factors, i.e. one level per
  # experimental "cell") so that "within-group" pairwise differences are
  # always genuine same-condition replicates, regardless of how many
  # design terms there are. LPE_ANOVA_var() itself is unchanged.

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
  # 3. between-group signal + p-value
  # -----------------------------

  if (variance.eval == "grand_mean") {

    # Type II sum-of-squares per design term (order-independent because the
    # design is additive/main-effects-only): for each term, compare the
    # full model against the model with only that term dropped.
    # SS_term = RSS(reduced) - RSS(full); numerically identical to the
    # classical sum(n_i * (group_means - grand_mean)^2) group-mean formula
    # when there is a single design term (one-way ANOVA).

    factor_df <- as.data.frame(factors)
    colnames(factor_df) <- design.terms

    full_formula <- stats::reformulate(paste0("`", design.terms, "`"))  # backquote: allows non-syntactic names
    X_full <- stats::model.matrix(full_formula, data = factor_df)

    tt_full <- stats::terms(full_formula)
    term_labels <- attr(tt_full, "term.labels")

    fit_full <- stats::lm.fit(X_full, t(expr))
    RSS_full <- colSums(fit_full$residuals^2)
    rank_full <- fit_full$rank

    per_factor <- list()

    for (term_label in term_labels) {
      # stats::drop.terms() cannot drop down to a bare intercept (it errors
      # internally when the result would have zero terms), so the one-term
      # case (dropping the only design factor -> reduced model = ~1) is
      # special-cased directly instead of going through drop.terms().
      if (length(term_labels) == 1) {
        reduced_labels <- character(0)
      } else {
        drop_idx <- which(term_labels == term_label)
        tt_reduced <- stats::drop.terms(tt_full, dropx = drop_idx, keep.response = FALSE)
        reduced_labels <- attr(stats::terms(tt_reduced), "term.labels")
      }

      X_reduced <- if (length(reduced_labels) == 0) {
        matrix(1, nrow = nrow(X_full), ncol = 1, dimnames = list(NULL, "(Intercept)"))
      } else {
        stats::model.matrix(stats::reformulate(reduced_labels), data = factor_df)
      }

      fit_reduced <- stats::lm.fit(X_reduced, t(expr))
      RSS_reduced <- colSums(fit_reduced$residuals^2)
      rank_reduced <- fit_reduced$rank

      df_i <- rank_full - rank_reduced

      if (df_i <= 0) {
        stop(
          "Design term '", term_label, "' is fully confounded with the ",
          "other design terms (e.g. missing factor-level combinations). ",
          "Its main effect cannot be tested."
        )
      }

      SS_i <- pmax(RSS_reduced - RSS_full, 0)
      MS_i <- SS_i / df_i
      Fstat_i <- MS_i / pred.var

      if (p.method == "chisq") {
        p_i <- stats::pchisq(df_i * Fstat_i, df = df_i, lower.tail = FALSE)
      } else {
        p_i <- stats::pf(Fstat_i, df1 = df_i, df2 = 1e6, lower.tail = FALSE)
      }

      p_i[!is.finite(p_i)] <- 1
      p_i <- pmax(p_i, .Machine$double.xmin)
      q_i <- stats::p.adjust(p_i, method = "BH")

      per_factor[[term_label]] <- list(
        df = df_i, MS_between = MS_i, F = Fstat_i, p.value = p_i, q.value = q_i
      )
    }

  } else {

    n_i_vec <- as.numeric(n_i)
    w <- sweep(1 / group.vars, 2, n_i_vec, "*")
    weighted.mean <- rowSums(w * group.means) / rowSums(w)
    T.stat <- rowSums(w * (group.means - weighted.mean)^2)

    if (p.method == "chisq") {
      p.val <- stats::pchisq(T.stat, df = k - 1, lower.tail = FALSE)
    } else {
      p.val <- stats::pf(T.stat / (k - 1), df1 = k - 1, df2 = 1e6, lower.tail = FALSE)
    }

    p.val[!is.finite(p.val)] <- 1
    p.val <- pmax(p.val, .Machine$double.xmin)

    adj.p <- stats::p.adjust(p.val, method = "BH")
  }

  # -----------------------------
  # 4. output
  # -----------------------------

  if (verbose) {
    cat("Finished.\n")
  }

  gene_id <- rownames(expr)

  if (is.null(gene_id)) {
    gene_id <- paste0("gene_", seq_len(nrow(expr)))
  }

  if (variance.eval == "grand_mean") {

    if (length(term_labels) == 1) {

      fac <- per_factor[[1]]

      res <- data.frame(
        gene = gene_id, mean = gene.mean, var = pred.var,
        MS_between = fac$MS_between, F = fac$F,
        p.value = fac$p.value, q.value = fac$q.value,
        method = "LPE", row.names = NULL
      )

    } else {

      res <- data.frame(gene = gene_id, mean = gene.mean, var = pred.var, row.names = NULL)

      for (term_label in term_labels) {
        fac <- per_factor[[term_label]]
        term_name <- gsub("^`|`$", "", term_label)  # drop backquotes of non-syntactic names
        res[[paste0("MS_between_", term_name)]] <- fac$MS_between
        res[[paste0("F_", term_name)]]          <- fac$F
        res[[paste0("p.value_", term_name)]]    <- fac$p.value
        res[[paste0("q.value_", term_name)]]    <- fac$q.value
      }

      res$method <- "LPE"
    }

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
  attr(res, "design.terms") <- design.terms
  attr(res, "trim.info") <- trim.info
  attr(res, "trend.info") <- trend.info
  attr(res, "base.var") <- base.var
  attr(res, "var.spline") <- var.result

  return(res)
}
