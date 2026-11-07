## Results on the bundled test data must match the shipped expected results
## (regenerate with inst/testdat/make_expected_results.R after an intended change).

testdat <- system.file("testdat", package = "LPEseq2")
read_counts <- function(f) as.matrix(read.csv(file.path(testdat, f), row.names = 1, check.names = FALSE))
read_meta <- function(f) read.csv(file.path(testdat, f), row.names = 1)

test_that("two-group test data reproduces the expected results", {
  prep <- LPE_preprocess(read_counts("LPEseq2_2group_test_counts.csv"),
                         read_meta("LPEseq2_2group_test_metadata.csv"),
                         design = ~ group, verbose = FALSE)
  res <- suppressMessages(LPE_ANOVA(prep, verbose = FALSE))
  expected <- read.csv(file.path(testdat, "LPEseq2_2group_expected_results.csv"))
  expect_identical(res$gene, expected$gene)
  expect_equal(res$p.value, expected$p.value, tolerance = 1e-8)
  expect_equal(res$q.value, expected$q.value, tolerance = 1e-8)
})

test_that("multi-way test data reproduces the expected results", {
  prep <- LPE_preprocess(read_counts("LPEseq2_test_counts.csv"),
                         read_meta("LPEseq2_test_metadata.csv"),
                         design = ~ group + batch, verbose = FALSE)
  res <- suppressMessages(LPE_ANOVA(prep, verbose = FALSE))
  expected <- read.csv(file.path(testdat, "LPEseq2_multiway_expected_results.csv"))
  expect_identical(res$gene, expected$gene)
  expect_equal(res$p.value_group, expected$p.value_group, tolerance = 1e-8)
  expect_equal(res$p.value_batch, expected$p.value_batch, tolerance = 1e-8)
})
