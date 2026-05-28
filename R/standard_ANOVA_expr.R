#' @keywords internal

standard_ANOVA_expr <- function(expr,
                                group,
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

  group <- droplevels(as.factor(group))

  if (ncol(expr) != length(group)) {
    stop("The number of columns in expr must match the length of group")
  }

  k <- nlevels(group)
  n <- length(group)
  n_i <- table(group)

  if (k < 2) {
    stop("At least two groups are required")
  }

  if (any(n_i < 2)) {
    stop("Standard one-way ANOVA requires at least two samples per group")
  }

  df_between <- k - 1
  df_within <- n - k

  if (df_within <= 0) {
    stop("Residual degrees of freedom must be positive for standard ANOVA")
  }

  gene.mean <- rowMeans(expr, na.rm = TRUE)

  stat_mat <- t(apply(expr, 1, function(x) {

    group_means <- tapply(x, group, mean)
    grand_mean <- sum(as.numeric(n_i) * group_means) / sum(as.numeric(n_i))

    ss_between <- sum(as.numeric(n_i) * (group_means - grand_mean)^2)

    ss_within <- sum(
      tapply(x, group, function(z) {
        sum((z - mean(z))^2)
      })
    )

    ms_between <- ss_between / df_between
    ms_within <- ss_within / df_within

    if (!is.finite(ms_within) || ms_within <= 0) {
      ms_within <- .Machine$double.xmin
    }

    f_stat <- ms_between / ms_within

    p_val <- stats::pf(
      f_stat,
      df1 = df_between,
      df2 = df_within,
      lower.tail = FALSE
    )

    c(
      var = ms_within,
      MS_between = ms_between,
      MS_within = ms_within,
      F = f_stat,
      df1 = df_between,
      df2 = df_within,
      p.value = p_val
    )
  }))

  p.val <- stat_mat[, "p.value"]
  p.val[!is.finite(p.val)] <- 1
  p.val <- pmax(p.val, .Machine$double.xmin)

  adj.p <- stats::p.adjust(p.val, method = p.adjust.method)

  gene_id <- rownames(expr)

  if (is.null(gene_id)) {
    gene_id <- paste0("gene_", seq_len(nrow(expr)))
  }

  res <- data.frame(
    gene = gene_id,
    mean = gene.mean,
    var = stat_mat[, "var"],
    MS_between = stat_mat[, "MS_between"],
    MS_within = stat_mat[, "MS_within"],
    F = stat_mat[, "F"],
    df1 = stat_mat[, "df1"],
    df2 = stat_mat[, "df2"],
    p.value = p.val,
    q.value = adj.p,
    method = "standard_anova",
    row.names = NULL
  )

  return(res)
}
