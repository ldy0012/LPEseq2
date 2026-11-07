## Reproduces the expected results shipped in inst/testdat/.
## Run from the package root after installing LPEseq2:
##   Rscript inst/testdat/make_expected_results.R

library(LPEseq2)

testdat <- system.file("testdat", package = "LPEseq2")
read_counts <- function(f) as.matrix(read.csv(file.path(testdat, f), row.names = 1, check.names = FALSE))
read_meta <- function(f) read.csv(file.path(testdat, f), row.names = 1)

## 1. Two-group design (6 Control vs 6 Treatment_A), ~ group
counts2 <- read_counts("LPEseq2_2group_test_counts.csv")
meta2 <- read_meta("LPEseq2_2group_test_metadata.csv")
prep2 <- LPE_preprocess(counts2, meta2, design = ~ group, verbose = FALSE)
res2 <- LPE_ANOVA(prep2, verbose = FALSE)

## 2. Three groups with a batch factor (6 samples per group), ~ group + batch
counts3 <- read_counts("LPEseq2_test_counts.csv")
meta3 <- read_meta("LPEseq2_test_metadata.csv")
prep3 <- LPE_preprocess(counts3, meta3, design = ~ group + batch, verbose = FALSE)
res3 <- LPE_ANOVA(prep3, verbose = FALSE)

out_dir <- if (dir.exists("inst/testdat")) "inst/testdat" else "."
write.csv(res2, file.path(out_dir, "LPEseq2_2group_expected_results.csv"), row.names = FALSE)
write.csv(res3, file.path(out_dir, "LPEseq2_multiway_expected_results.csv"), row.names = FALSE)

cat("Two-group: ", sum(res2$q.value < 0.05), "of", nrow(res2), "genes with q < 0.05\n")
cat("Multi-way, group: ", sum(res3$q.value_group < 0.05), "of", nrow(res3), "genes with q < 0.05\n")
cat("Multi-way, batch: ", sum(res3$q.value_batch < 0.05), "of", nrow(res3), "genes with q < 0.05\n")
