#' Estimate Intensity-Dependent Variance Trend for LPE-ANOVA
#'
#' This function estimates an intensity-dependent variance trend using
#' within-group pairwise differences. Optionally, weighted between-group
#' differences can be included when replicate information is limited.
#'
#' @param expr A numeric expression matrix with genes as rows and samples as columns.
#' @param group A factor or vector indicating sample group labels.
#' @param n.bin Number of quantile bins used for variance trend estimation.
#' @param df Degrees of freedom for smoothing spline.
#' @param trim.method Outlier trimming method. One of \code{"mad"}, \code{"quantile"}, or \code{"fixed"}.
#' @param d Fixed trimming threshold used when \code{trim.method = "fixed"}.
#' @param use_weighted_between Logical. Whether to include weighted between-group differences.
#'
#' @return A \code{smooth.spline} object representing the estimated variance trend.
#'
#' @export
LPE_ANOVA_var <- function(expr,
                          group,
                          n.bin = 100,
                          df = 10,
                          trim.method = c("mad", "quantile", "fixed"),
                          d = 1.2,
                          use_weighted_between = FALSE) {

  trim.method <- match.arg(trim.method)

  # -----------------------------
  # 1. input check
  # -----------------------------

  if (!is.matrix(expr)) {
    expr <- as.matrix(expr)
  }

  if (!is.numeric(expr)) {
    stop("expr must be a numeric matrix")
  }

  if (anyNA(expr) || any(!is.finite(expr))) {
    stop("expr contains NA, NaN, or infinite values")
  }

  if (ncol(expr) != length(group)) {
    stop("The number of columns in expr must match the length of group")
  }

  if (!is.numeric(n.bin) || length(n.bin) != 1 || n.bin < 5) {
    stop("n.bin must be a single numeric value >= 5")
  }

  if (!is.numeric(df) || length(df) != 1 || df < 2) {
    stop("df must be a single numeric value >= 2")
  }

  if (!is.numeric(d) || length(d) != 1 || d <= 0) {
    stop("d must be a single positive numeric value")
  }

  if (!is.logical(use_weighted_between) || length(use_weighted_between) != 1) {
    stop("use_weighted_between must be TRUE or FALSE")
  }

  group <- droplevels(as.factor(group))
  k <- nlevels(group)
  n_i <- table(group)

  if (k < 2) {
    stop("At least two groups are required")
  }

  split_index <- split(seq_along(group), group)

  M_all <- numeric(0)
  A_all <- numeric(0)
  W_all <- numeric(0)

  # -----------------------------
  # 2. within-group pairwise differences
  # -----------------------------

  for (g in seq_len(nrow(expr))) {
    y <- expr[g, ]

    for (i in seq_len(k)) {
      idx <- split_index[[i]]
      ni <- length(idx)

      if (ni >= 2) {
        comb <- utils::combn(idx, 2)

        yi <- y[comb[1, ]]
        yj <- y[comb[2, ]]

        M_all <- c(M_all, (yi - yj) / sqrt(2))
        A_all <- c(A_all, (yi + yj) / 2)
        W_all <- c(W_all, rep(1, ncol(comb)))
      }
    }
  }

  # -----------------------------
  # 3. optional weighted between-group information
  # -----------------------------

  if (use_weighted_between) {
    warning(
      "Between-group differences are used for variance training. ",
      "This may inflate variance if many genes are truly differentially expressed."
    )

    for (g in seq_len(nrow(expr))) {
      y <- expr[g, ]
      group_means <- sapply(split_index, function(idx) mean(y[idx]))
      comb_groups <- utils::combn(seq_len(k), 2)

      for (j in seq_len(ncol(comb_groups))) {
        i1 <- comb_groups[1, j]
        i2 <- comb_groups[2, j]

        n1 <- as.numeric(n_i[i1])
        n2 <- as.numeric(n_i[i2])

        if (n1 == 1 && n2 == 1) {
          next
        }

        m_star <- (group_means[i1] - group_means[i2]) /
          sqrt(1 / n1 + 1 / n2)

        a_val <- (group_means[i1] + group_means[i2]) / 2

        alpha <- min(n1, n2) / (n1 + n2)

        M_all <- c(M_all, m_star)
        A_all <- c(A_all, a_val)
        W_all <- c(W_all, alpha)
      }
    }
  }

  # -----------------------------
  # 4. full non-replicate fallback
  # -----------------------------

  if (length(M_all) == 0 && all(n_i == 1)) {
    warning(
      "All groups have only one sample. ",
      "Variance estimation relies entirely on between-group differences. ",
      "P-values should be interpreted cautiously."
    )

    for (g in seq_len(nrow(expr))) {
      y <- expr[g, ]
      comb_g <- utils::combn(seq_len(k), 2)

      yi <- y[comb_g[1, ]]
      yj <- y[comb_g[2, ]]

      M_all <- c(M_all, (yi - yj) / sqrt(2))
      A_all <- c(A_all, (yi + yj) / 2)
      W_all <- c(W_all, rep(1, ncol(comb_g)))
    }
  }

  if (length(M_all) < 10) {
    stop("Too few pairwise differences to estimate variance trend")
  }

  valid <- is.finite(M_all) & is.finite(A_all) & is.finite(W_all) & W_all > 0
  M_all <- M_all[valid]
  A_all <- A_all[valid]
  W_all <- W_all[valid]

  if (length(M_all) < 10) {
    stop("Too few valid pairwise differences to estimate variance trend")
  }

  # -----------------------------
  # 5. robust outlier trimming
  # -----------------------------

  if (trim.method == "mad") {

    med <- stats::median(M_all, na.rm = TRUE)
    s <- stats::mad(M_all, center = med, na.rm = TRUE)

    if (is.finite(s) && s > 0) {
      keep <- abs(M_all - med) <= 3 * s
    } else {
      keep <- rep(TRUE, length(M_all))
    }

  } else if (trim.method == "quantile") {

    lo <- stats::quantile(M_all, 0.01, na.rm = TRUE)
    hi <- stats::quantile(M_all, 0.99, na.rm = TRUE)
    keep <- M_all >= lo & M_all <= hi

  } else {

    keep <- abs(M_all) < d
  }

  M_all <- M_all[keep]
  A_all <- A_all[keep]
  W_all <- W_all[keep]

  if (length(M_all) < 10) {
    stop("Too few observations after trimming")
  }

  # -----------------------------
  # 6. quantile binning
  # -----------------------------

  n.bin <- min(n.bin, floor(length(A_all) / 5))

  if (n.bin < 5) {
    stop("Too few observations for quantile binning")
  }

  probs <- seq(0, 1, length.out = n.bin + 1)
  quantile.A <- unique(as.numeric(stats::quantile(A_all, probs = probs, na.rm = TRUE)))

  if (length(quantile.A) < 5) {
    stop("Too few unique A bins. Check input data or filtering")
  }

  var.M <- numeric(length(quantile.A) - 1)
  medianAs <- numeric(length(quantile.A) - 1)

  var.M[] <- NA_real_
  medianAs[] <- NA_real_

  for (i in 2:length(quantile.A)) {

    if (i == length(quantile.A)) {
      idx <- which(A_all >= quantile.A[i - 1] & A_all <= quantile.A[i])
    } else {
      idx <- which(A_all >= quantile.A[i - 1] & A_all < quantile.A[i])
    }

    if (length(idx) > 1) {
      w <- W_all[idx]
      x <- M_all[idx]

      w_mean <- sum(w * x) / sum(w)
      var.M[i - 1] <- sum(w * (x - w_mean)^2) / sum(w)
      medianAs[i - 1] <- stats::median(A_all[idx])
    }
  }

  base.var <- data.frame(
    A = medianAs,
    var.M = var.M
  )

  base.var <- base.var[
    is.finite(base.var$A) &
      is.finite(base.var$var.M) &
      base.var$var.M > 0,
  ]

  if (nrow(base.var) < 5) {
    stop("Too few valid bins to fit variance trend")
  }

  # -----------------------------
  # 7. low-expression variance floor
  # -----------------------------

  max_var <- max(base.var$var.M, na.rm = TRUE)
  idx_max <- which(base.var$var.M == max_var)[1]

  base.var$var.M[seq_len(idx_max)] <- max_var

  # -----------------------------
  # 8. smoothing spline
  # -----------------------------

  df_use <- min(df, nrow(base.var) - 1)

  sm.spline <- stats::smooth.spline(
    x = base.var$A,
    y = base.var$var.M,
    df = df_use
  )

  return(sm.spline)
}
