#' @keywords internal

standard_ANOVA_expr <- function(expr,
                                factors,
                                p.adjust.method = "BH") {

  if (!is.matrix(expr)) {
    expr <- as.matrix(expr)
  }

  if (!is.numeric(expr)) {
    stop("expr must be a numeric matrix")
  }

  if (anyNA(expr) || any(!is.finite(expr))) {
    stop("expr contains NA, NaN, or infinite values")
  }

  if (!is.list(factors) || length(factors) < 1 || is.null(names(factors))) {
    stop("factors must be a named list of factors (one per design term)")
  }

  factor_names <- names(factors)
  n <- ncol(expr)

  for (nm in factor_names) {
    factors[[nm]] <- droplevels(as.factor(factors[[nm]]))
    if (length(factors[[nm]]) != n) {
      stop("The number of columns in expr must match the length of each factor")
    }
  }

  # -----------------------------
  # design matrix (main effects only, no interaction)
  # -----------------------------

  factor_df <- as.data.frame(factors)
  colnames(factor_df) <- factor_names

  full_formula <- stats::reformulate(paste0("`", factor_names, "`"))  # backquote: allows non-syntactic names
  X_full <- stats::model.matrix(full_formula, data = factor_df)

  tt_full <- stats::terms(full_formula)
  term_labels <- attr(tt_full, "term.labels")

  gene.mean <- rowMeans(expr, na.rm = TRUE)

  fit_full <- stats::lm.fit(X_full, t(expr))
  RSS_full <- colSums(fit_full$residuals^2)
  rank_full <- fit_full$rank
  df_within <- n - rank_full

  if (df_within <= 0) {
    stop("Residual degrees of freedom must be positive for standard ANOVA")
  }

  MS_within <- RSS_full / df_within
  MS_within[!is.finite(MS_within) | MS_within <= 0] <- .Machine$double.xmin

  # -----------------------------
  # Type II sum-of-squares per design term:
  # SS_term = RSS(model without term) - RSS(full model)
  # (order-independent because the design is additive/main-effects-only)
  # -----------------------------

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

    df1 <- rank_full - rank_reduced

    if (df1 <= 0) {
      stop(
        "Design term '", term_label, "' is fully confounded with the other ",
        "design terms (e.g. missing factor-level combinations). Its main ",
        "effect cannot be tested."
      )
    }

    SS_i <- pmax(RSS_reduced - RSS_full, 0)
    MS_i <- SS_i / df1
    F_i  <- MS_i / MS_within

    p_i <- stats::pf(F_i, df1 = df1, df2 = df_within, lower.tail = FALSE)
    p_i[!is.finite(p_i)] <- 1
    p_i <- pmax(p_i, .Machine$double.xmin)

    q_i <- stats::p.adjust(p_i, method = p.adjust.method)

    per_factor[[term_label]] <- list(
      df1 = df1, MS_between = MS_i, F = F_i, p.value = p_i, q.value = q_i
    )
  }

  gene_id <- rownames(expr)

  if (is.null(gene_id)) {
    gene_id <- paste0("gene_", seq_len(nrow(expr)))
  }

  # -----------------------------
  # output: one-way designs keep the original (unsuffixed) column names for
  # backward compatibility; multi-way designs get one suffixed set of
  # columns per design term ("wide" format).
  # -----------------------------

  if (length(term_labels) == 1) {

    fac <- per_factor[[1]]

    res <- data.frame(
      gene = gene_id,
      mean = gene.mean,
      var = MS_within,
      MS_between = fac$MS_between,
      MS_within = MS_within,
      F = fac$F,
      df1 = fac$df1,
      df2 = df_within,
      p.value = fac$p.value,
      q.value = fac$q.value,
      method = "standard_anova",
      row.names = NULL
    )

  } else {

    res <- data.frame(
      gene = gene_id,
      mean = gene.mean,
      var = MS_within,
      MS_within = MS_within,
      df2 = df_within,
      row.names = NULL
    )

    for (term_label in term_labels) {
      fac <- per_factor[[term_label]]
      term_name <- gsub("^`|`$", "", term_label)  # drop backquotes of non-syntactic names
      res[[paste0("MS_between_", term_name)]] <- fac$MS_between
      res[[paste0("F_", term_name)]]          <- fac$F
      res[[paste0("df1_", term_name)]]        <- fac$df1
      res[[paste0("p.value_", term_name)]]    <- fac$p.value
      res[[paste0("q.value_", term_name)]]    <- fac$q.value
    }

    res$method <- "standard_anova"
  }

  attr(res, "design.terms") <- term_labels

  return(res)
}
