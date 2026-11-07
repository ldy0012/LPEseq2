## variance.trend = "per_group" (group-specific variance trends, LPEseq1-style)

make_obj_pg <- function(ns, size = rep(20, length(ns)), seed = 1, n_gene = 1500) {
  set.seed(seed)
  g <- factor(rep(paste0("G", seq_along(ns)), ns))
  mu <- exp(stats::rnorm(n_gene, 4, 1.5))
  # `size` per group: smaller size = more overdispersed (noisier) group
  cnt <- vapply(seq_along(g), function(j) {
    stats::rnbinom(n_gene, mu = mu, size = size[as.integer(g[j])])
  }, numeric(n_gene))
  dimnames(cnt) <- list(paste0("g", seq_len(n_gene)), paste0("s", seq_along(g)))
  cd <- data.frame(group = g, row.names = colnames(cnt))
  suppressWarnings(LPE_preprocess(cnt, cd, ~ group, normalize.method = "library_size",
                                  min.count = 0, verbose = FALSE))
}
run_pg <- function(obj, ...) suppressMessages(LPE_ANOVA(obj, verbose = FALSE, ...))

test_that("default variance.trend = 'pooled' reproduces the previous per_group result", {
  obj <- make_obj_pg(c(3, 3))
  a <- run_pg(obj, variance.eval = "per_group")
  b <- run_pg(obj, variance.eval = "per_group", variance.trend = "pooled")
  expect_identical(a$p.value, b$p.value)
  expect_identical(attr(a, "variance.settings")$variance.trend, "pooled")
  expect_true(all(attr(a, "variance.settings")$trend.source == "pooled"))
  expect_null(attr(a, "var.spline.group"))
})

test_that("grand_mean results are unchanged and carry no group trends", {
  obj <- make_obj_pg(c(3, 3, 3))
  a <- run_pg(obj)
  expect_identical(attr(a, "variance.settings")$variance.trend, "pooled")
  expect_null(attr(a, "variance.settings")$trend.source)
  expect_null(attr(a, "var.spline.group"))
})

test_that("replicated groups get their own trend, fitted from their own pairs", {
  obj <- make_obj_pg(c(3, 3))
  a <- run_pg(obj, variance.eval = "per_group", variance.trend = "per_group")
  src <- attr(a, "variance.settings")$trend.source
  expect_identical(unname(src), c("own", "own"))
  sp <- attr(a, "var.spline.group")
  expect_named(sp, c("G1", "G2"))

  # the stored trend of G1 equals LPE_ANOVA_var() on G1's columns alone
  idx <- which(obj$group == "G1")
  ref <- suppressMessages(LPE_ANOVA_var(obj$expr[, idx], factor(rep("G1", length(idx))),
                                        trim.method = "none",
                                        use_weighted_between = FALSE))
  expect_equal(sp$G1$object$y, ref$object$y)
  expect_identical(attr(sp$G1, "trim.info")$method, "none")
})

test_that("single-sample groups fall back to the pooled trend", {
  obj <- make_obj_pg(c(3, 1, 3))
  a <- run_pg(obj, variance.eval = "per_group", variance.trend = "per_group")
  src <- attr(a, "variance.settings")$trend.source
  expect_identical(unname(src), c("own", "pooled", "own"))
  expect_null(attr(a, "var.spline.group")$G2)
  expect_true(all(is.finite(a$p.value)))
})

test_that("all single-sample groups give the same result as the pooled trend", {
  obj <- make_obj_pg(c(1, 1, 1))
  a <- suppressWarnings(run_pg(obj, variance.eval = "per_group", variance.trend = "per_group"))
  b <- suppressWarnings(run_pg(obj, variance.eval = "per_group", variance.trend = "pooled"))
  expect_identical(a$p.value, b$p.value)
  expect_true(all(attr(a, "variance.settings")$trend.source == "pooled"))
})

test_that("a noisier group gets a higher own trend", {
  obj <- make_obj_pg(c(4, 4), size = c(50, 3))   # G2 much more overdispersed
  a <- run_pg(obj, variance.eval = "per_group", variance.trend = "per_group")
  sp <- attr(a, "var.spline.group")
  grid <- seq(4, 8, length.out = 20)
  v1 <- predict(sp$G1$object, grid)$y
  v2 <- predict(sp$G2$object, grid)$y
  expect_gt(mean(v2), mean(v1))
})

test_that("explicit trim.method is applied to the group trends", {
  obj <- make_obj_pg(c(3, 3))
  a <- run_pg(obj, variance.eval = "per_group", variance.trend = "per_group",
              trim.method = "iqr")
  expect_identical(attr(attr(a, "var.spline.group")$G1, "trim.info")$method, "iqr")
})

test_that("invalid combinations are rejected", {
  obj <- make_obj_pg(c(3, 3))
  expect_error(run_pg(obj, variance.trend = "per_group"),
               "only available with variance.eval = 'per_group'")
})

test_that("LPE_ANOVA_var accepts a single replicated group only without between pairs", {
  obj <- make_obj_pg(c(3, 3))
  idx <- which(obj$group == "G1")
  g1 <- factor(rep("G1", length(idx)))
  expect_error(LPE_ANOVA_var(obj$expr[, idx], g1, use_weighted_between = TRUE),
               "At least two groups")
  expect_error(LPE_ANOVA_var(obj$expr[, idx[1], drop = FALSE], g1[1],
                             use_weighted_between = FALSE),
               "at least two samples")
  expect_s3_class(
    suppressMessages(LPE_ANOVA_var(obj$expr[, idx], g1, trim.method = "none",
                                   use_weighted_between = FALSE))$object,
    "smooth.spline")
})
