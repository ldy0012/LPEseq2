## trim.method = "auto" / use_weighted_between = "auto" (LPEseq2 >= 0.1.0)

make_obj <- function(ns, seed = 1, n_gene = 1500) {
  set.seed(seed)
  g <- factor(rep(paste0("G", seq_along(ns)), ns))
  mu <- exp(stats::rnorm(n_gene, 4, 1.5))
  cnt <- vapply(seq_along(g), function(j) stats::rnbinom(n_gene, mu = mu, size = 20),
                numeric(n_gene))
  dimnames(cnt) <- list(paste0("g", seq_len(n_gene)), paste0("s", seq_along(g)))
  cd <- data.frame(group = g, row.names = colnames(cnt))
  suppressWarnings(LPE_preprocess(cnt, cd, ~ group, normalize.method = "library_size",
                                  min.count = 0, verbose = FALSE))
}
run <- function(obj, ...) suppressMessages(LPE_ANOVA(obj, verbose = FALSE, ...))

test_that("all cells replicated -> within pairs, no trimming", {
  obj <- make_obj(c(3, 3))
  a <- run(obj)
  m <- run(obj, trim.method = "none", use_weighted_between = FALSE)
  vs <- attr(a, "variance.settings")
  expect_identical(vs$trim.method, "none")
  expect_false(vs$use_weighted_between)
  expect_true(vs$auto)
  expect_identical(a$p.value, m$p.value)
})

test_that("a single-sample cell -> between pairs with dvalue trimming", {
  for (ns in list(c(1, 1), c(3, 1))) {
    obj <- make_obj(ns)
    a <- run(obj)
    m <- run(obj, trim.method = "dvalue", use_weighted_between = TRUE)
    expect_identical(attr(a, "variance.settings")$trim.method, "dvalue")
    expect_true(attr(a, "variance.settings")$use_weighted_between)
    expect_identical(a$p.value, m$p.value)
  }
})

test_that("explicit values are used unchanged", {
  obj <- make_obj(c(3, 3))
  vs <- attr(run(obj, trim.method = "dvalue", use_weighted_between = TRUE), "variance.settings")
  expect_false(vs$auto)
  expect_identical(vs$trim.method, "dvalue")
  expect_true(vs$use_weighted_between)
})

test_that("invalid use_weighted_between is rejected", {
  expect_error(run(make_obj(c(2, 2)), use_weighted_between = "yes"), "TRUE, FALSE or")
})
