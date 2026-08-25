#' Estimate Intensity-Dependent Variance Trend for LPE-ANOVA
#'
#' This function estimates an intensity-dependent variance trend using
#' within-group pairwise differences. Following LPEseq1 (Gim et al. 2016),
#' each raw pairwise difference D is symmetrized to \eqn{\pm D} before
#' quantile binning, so that the pooled difference distribution used for
#' variance estimation is centered at zero by construction rather than by
#' subtracting an empirically estimated mean. Optionally, weighted
#' between-group differences can be included when replicate information is
#' limited. When trimming is enabled, the conventional 1.5*IQR rule is
#' applied within expression-intensity bins to pairwise values used for
#' variance trend estimation.
#'
#'
#' @param expr A numeric expression matrix with genes as rows and samples as columns.
#' @param group A factor or vector indicating sample group labels.
#' @param n.bin Number of quantile bins used for variance trend estimation.
#' @param df Degrees of freedom for smoothing spline.
#' @param trim.method Outlier trimming method applied to pairwise values
#'   used for variance trend estimation. One of \code{"iqr"}, \code{"dvalue"},
#'   or \code{"none"}. \code{"iqr"} applies the conventional 1.5*IQR boxplot
#'   rule within expression-intensity A-bins after pooling within-group and
#'   between-group values. \code{"dvalue"} applies a fixed global threshold
#'   to the raw pairwise difference D (as in LPEseq1's non-replicate outlier
#'   procedure, Gim et al. 2016): any pairwise value with |D| > d.threshold
#'   is removed, independent of bin. \code{"none"} performs no outlier
#'   trimming.
#' @param use_weighted_between Logical. Whether to include weighted between-group
#'   differences in variance trend estimation.
#' @param d.threshold Numeric. Fixed threshold applied to the raw pairwise
#'   difference D when \code{trim.method = "dvalue"}. Default is 1.2,
#'   matching the default reported in LPEseq1 (Gim et al. 2016). Ignored for
#'   other trim.method values.
#' @note Pairwise D/M values are symmetrized (\eqn{\pm D}) before trimming
#'   and quantile binning, following LPEseq1's approach of fixing the
#'   difference distribution's center at zero. Consequently, the pairwise
#'   counts reported in the returned \code{trim.info} (e.g.
#'   \code{n_total_before}, \code{n_within_before}, \code{n_between_before},
#'   and their \verb{_after}/\verb{_removed} counterparts) are twice the
#'   number of underlying replicate/group pairs.
#'
#' @export
LPE_ANOVA_var <- function(expr,
                          group,
                          n.bin = 100,
                          df = 10,
                          trim.method = c("iqr", "dvalue", "none"),
                          use_weighted_between = FALSE,
                          d.threshold = 1.2) {

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

  if (!is.logical(use_weighted_between) || length(use_weighted_between) != 1) {
    stop("use_weighted_between must be TRUE or FALSE")
  }

  group <- droplevels(as.factor(group))
  k <- nlevels(group)
  n_i <- table(group)

  if (k < 2) {
    stop("At least two groups are required")
  }

  if (trim.method == "dvalue") {
    if (!is.numeric(d.threshold) || length(d.threshold) != 1 || d.threshold <= 0) {
      stop("d.threshold must be a single positive numeric value")
    }
  }

  split_index <- split(seq_along(group), group)

  # -----------------------------
  # Initialize containers by source
  # -----------------------------

  M_within <- numeric(0)
  A_within <- numeric(0)
  W_within <- numeric(0)
  D_within <- numeric(0)

  M_between <- numeric(0)
  A_between <- numeric(0)
  W_between <- numeric(0)
  D_between <- numeric(0)

  # -----------------------------
  # 2. within-group pairwise differences
  #    These are genuine replicate-based residual differences.
  #    Following LPEseq1 (Gim et al. 2016), each raw difference D is
  #    symmetrized to +-D (mirrored with the same A and weight) so that
  #    the pooled M distribution used for variance estimation is centered
  #    at zero by construction. As a result, all pairwise counts recorded
  #    below (and in trim.info) are twice the number of underlying pairs.
  #    If trim.method = "iqr", they will be pooled with between-derived
  #    values and trimmed within expression-intensity A-bins later.
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

        d_raw <- yi - yj
        a_val <- (yi + yj) / 2
        w_val <- rep(1, ncol(comb))

        D_within <- c(D_within, d_raw, -d_raw)
        M_within <- c(M_within, d_raw / sqrt(2), -d_raw / sqrt(2))
        A_within <- c(A_within, a_val, a_val)
        W_within <- c(W_within, w_val, w_val)
      }
    }
  }

  # Keep only valid within-group values
  valid_within <- is.finite(M_within) &
    is.finite(A_within) &
    is.finite(W_within) &
    is.finite(D_within) &
    W_within > 0

  M_within <- M_within[valid_within]
  A_within <- A_within[valid_within]
  W_within <- W_within[valid_within]
  D_within <- D_within[valid_within]

  # -----------------------------
  # 3. optional weighted between-group information
  #    These values may contain true DE signal.
  #    If trim.method = "iqr", they will be pooled with within-group
  #    pairwise values and trimmed within expression-intensity A-bins later.
  # -----------------------------

  if (use_weighted_between) {
    warning(
      "Between-group differences are used for variance training. ",
      "If trim.method = 'iqr', within-group and between-group pairwise values ",
      "will be pooled and trimmed within expression-intensity A-bins using ",
      "the conventional 1.5*IQR rule."
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

        d_raw <- group_means[i1] - group_means[i2]

        m_star <- d_raw / sqrt(1 / n1 + 1 / n2)

        a_val <- (group_means[i1] + group_means[i2]) / 2

        alpha <- min(n1, n2) / (n1 + n2)

        D_between <- c(D_between, d_raw, -d_raw)
        M_between <- c(M_between, m_star, -m_star)
        A_between <- c(A_between, a_val, a_val)
        W_between <- c(W_between, alpha, alpha)
      }
    }
  }

  # -----------------------------
  # 4. full non-replicate fallback
  #    If all groups have n = 1, variance estimation must rely
  #    on between-group differences; trim these values.
  # -----------------------------

  if (length(M_within) == 0 && all(n_i == 1)) {
    warning(
      "All groups have only one sample. ",
      "Variance estimation relies entirely on between-group differences. ",
      if (trim.method != "none") {
        "Outlier trimming will be applied to these between-group differences. "
      } else {
        "No outlier trimming will be applied because trim.method = 'none'. "
      },
      "P-values should be interpreted cautiously."
    )

    for (g in seq_len(nrow(expr))) {
      y <- expr[g, ]
      comb_g <- utils::combn(seq_len(k), 2)

      yi <- y[comb_g[1, ]]
      yj <- y[comb_g[2, ]]

      d_raw <- yi - yj
      a_val <- (yi + yj) / 2
      w_val <- rep(1, ncol(comb_g))

      D_between <- c(D_between, d_raw, -d_raw)
      M_between <- c(M_between, d_raw / sqrt(2), -d_raw / sqrt(2))
      A_between <- c(A_between, a_val, a_val)
      W_between <- c(W_between, w_val, w_val)
    }
  }

  # -----------------------------
  # 5. outlier trimming (pooled bin-wise IQR, or fixed D-value threshold)
  #    For trim.method = "iqr": pool within- and between-derived values,
  #    then apply the conventional 1.5*IQR rule within A-bins (on the M scale).
  #    For trim.method = "dvalue": apply a single fixed global threshold
  #    to the raw D values, independent of bin (see block below).
  # -----------------------------

  M_all <- c(M_within, M_between)
  A_all <- c(A_within, A_between)
  W_all <- c(W_within, W_between)
  D_all <- c(D_within, D_between)

  source_all <- c(
    rep("within", length(M_within)),
    rep("between", length(M_between))
  )

  valid_all0 <- is.finite(M_all) &
    is.finite(A_all) &
    is.finite(W_all) &
    is.finite(D_all) &
    W_all > 0

  M_all <- M_all[valid_all0]
  A_all <- A_all[valid_all0]
  W_all <- W_all[valid_all0]
  D_all <- D_all[valid_all0]
  source_all <- source_all[valid_all0]

  n_within_before <- sum(source_all == "within")
  n_between_before <- sum(source_all == "between")

  trim.rule <- switch(
    trim.method,
    iqr    = "pooled_bin_wise_boxplot_IQR_rule",
    dvalue = "fixed_D_threshold",
    none   = "none"
  )

  trim.scale <- switch(
    trim.method,
    iqr    = "M",
    dvalue = "D",
    none   = NA_character_
  )

  trim.info <- list(
    method = trim.method,
    rule = trim.rule,
    trim.scale = trim.scale,
    d.threshold = if (trim.method == "dvalue") d.threshold else NA_real_,
    n_total_before = length(M_all),
    n_total_after = length(M_all),
    n_total_removed = 0,
    n_within_before = n_within_before,
    n_within_after = n_within_before,
    n_within_removed = 0,
    n_between_before = n_between_before,
    n_between_after = n_between_before,
    n_between_removed = 0,
    threshold.table = data.frame()
  )

  if (trim.method == "iqr" && length(M_all) > 0) {

    keep_all <- rep(TRUE, length(M_all))
    threshold.list <- list()

    n.iqr.bin <- min(n.bin, floor(length(A_all) / 5))

    if (n.iqr.bin < 1) {
      n.iqr.bin <- 1
    }

    if (n.iqr.bin == 1) {
      breaks <- c(min(A_all, na.rm = TRUE), max(A_all, na.rm = TRUE))
    } else {
      probs.iqr <- seq(0, 1, length.out = n.iqr.bin + 1)

      breaks <- unique(as.numeric(stats::quantile(
        A_all,
        probs = probs.iqr,
        na.rm = TRUE
      )))
    }

    if (length(breaks) < 2) {
      breaks <- c(min(A_all, na.rm = TRUE), max(A_all, na.rm = TRUE))
    }

    for (i in 2:length(breaks)) {

      if (i == length(breaks)) {
        idx_bin <- which(A_all >= breaks[i - 1] & A_all <= breaks[i])
      } else {
        idx_bin <- which(A_all >= breaks[i - 1] & A_all < breaks[i])
      }

      if (length(idx_bin) == 0) {
        next
      }

      # Use M values for outlier detection because M is the value used
      # to estimate local pooled variance.
      X_bin <- M_all[idx_bin]

      if (length(X_bin) >= 5) {

        q1 <- as.numeric(stats::quantile(X_bin, 0.25, na.rm = TRUE))
        q3 <- as.numeric(stats::quantile(X_bin, 0.75, na.rm = TRUE))
        iqr_value <- q3 - q1

        if (is.finite(iqr_value) && iqr_value > 0) {
          lower <- q1 - 1.5 * iqr_value
          upper <- q3 + 1.5 * iqr_value

          keep_bin <- X_bin >= lower & X_bin <= upper
        } else {
          lower <- q1
          upper <- q3
          keep_bin <- rep(TRUE, length(X_bin))
        }

      } else {

        q1 <- NA_real_
        q3 <- NA_real_
        iqr_value <- NA_real_
        lower <- NA_real_
        upper <- NA_real_
        keep_bin <- rep(TRUE, length(X_bin))
      }

      keep_all[idx_bin] <- keep_bin

      source_bin <- source_all[idx_bin]

      threshold.list[[length(threshold.list) + 1]] <- data.frame(
        bin = i - 1,
        A_low = breaks[i - 1],
        A_high = breaks[i],
        Q1 = q1,
        Q3 = q3,
        IQR = iqr_value,
        lower_bound = lower,
        upper_bound = upper,
        n_values = length(idx_bin),
        n_within = sum(source_bin == "within"),
        n_between = sum(source_bin == "between"),
        n_removed = sum(!keep_bin),
        n_within_removed = sum(!keep_bin & source_bin == "within"),
        n_between_removed = sum(!keep_bin & source_bin == "between")
      )
    }

    M_all <- M_all[keep_all]
    A_all <- A_all[keep_all]
    W_all <- W_all[keep_all]
    D_all <- D_all[keep_all]
    source_all <- source_all[keep_all]

    if (length(threshold.list) > 0) {
      trim.info$threshold.table <- do.call(rbind, threshold.list)
    }

    trim.info$n_within_after <- sum(source_all == "within")
    trim.info$n_between_after <- sum(source_all == "between")

    trim.info$n_within_removed <- trim.info$n_within_before - trim.info$n_within_after
    trim.info$n_between_removed <- trim.info$n_between_before - trim.info$n_between_after

    trim.info$n_total_after <- length(M_all)
    trim.info$n_total_removed <- trim.info$n_total_before - trim.info$n_total_after
  }

  if (trim.method == "dvalue" && length(D_all) > 0) {

    # LPEseq1-style: a single fixed global threshold on the raw
    # pairwise difference D, applied uniformly regardless of bin.
    keep_all <- abs(D_all) <= d.threshold

    n_removed_within  <- sum(!keep_all & source_all == "within")
    n_removed_between <- sum(!keep_all & source_all == "between")

    M_all <- M_all[keep_all]
    A_all <- A_all[keep_all]
    W_all <- W_all[keep_all]
    D_all <- D_all[keep_all]
    source_all <- source_all[keep_all]

    trim.info$n_within_after  <- sum(source_all == "within")
    trim.info$n_between_after <- sum(source_all == "between")
    trim.info$n_within_removed  <- n_removed_within
    trim.info$n_between_removed <- n_removed_between
    trim.info$n_total_after   <- length(M_all)
    trim.info$n_total_removed <- trim.info$n_total_before - trim.info$n_total_after

    trim.info$threshold.table <- data.frame(
      d.threshold = d.threshold,
      n_total_before = trim.info$n_total_before,
      n_total_removed = trim.info$n_total_removed
    )
  }

  # -----------------------------
  # 6. validate variance-training information
  # -----------------------------

  if (length(M_all) < 10) {
    stop("Too few pairwise differences to estimate variance trend")
  }

  valid_all <- is.finite(M_all) &
    is.finite(A_all) &
    is.finite(W_all) &
    W_all > 0

  M_all <- M_all[valid_all]
  A_all <- A_all[valid_all]
  W_all <- W_all[valid_all]

  if (length(M_all) < 10) {
    stop("Too few valid pairwise differences to estimate variance trend")
  }

  # -----------------------------
  # 7. quantile binning
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
      var.M[i - 1] <- sum(w * x^2) / sum(w)
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

  if (length(unique(base.var$A)) < 4) {
    stop("Too few unique A values to fit variance trend")
  }

  # -----------------------------
  # 8. low-expression variance floor
  # -----------------------------

  max_var <- max(base.var$var.M, na.rm = TRUE)
  idx_max <- which(base.var$var.M == max_var)[1]

  base.var$var.M[seq_len(idx_max)] <- max_var

  # -----------------------------
  # 9. smoothing spline
  # -----------------------------

  df_use <- min(df, nrow(base.var) - 1)

  sm.spline <- stats::smooth.spline(
    x = base.var$A,
    y = base.var$var.M,
    df = df_use
  )

  trend.info <- list(
    method = "mean_spline",
    spline.df = df_use
  )

  result <- list(
    type   = "smooth.spline",
    object = sm.spline,
    x_min  = min(base.var$A),
    x_max  = max(base.var$A)
  )

  attr(result, "trim.info")  <- trim.info
  attr(result, "trend.info") <- trend.info
  attr(result, "base.var")   <- base.var

  return(result)
}
