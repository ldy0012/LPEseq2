# Supplementary Note S1 / Figure 2B
# GSE142647 SUM149PT: Rep1 vs Rep2 within DMSO and within TSA.
# Methods: LPEseq2 (pseudo-bulk), MAST with CDR adjustment, and Wilcoxon rank-sum.

library(Seurat)
library(LPEseq2)

input_file <- file.path("S1_SUM149PT", "data", "sum149pt_qc.rds")
output_dir <- file.path("S1_SUM149PT", "results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Load the QC filtered Seurat objects.
prepared <- readRDS(input_file)
summary_results <- data.frame()

for (condition in c("DMSO", "TSA")) {
  message("Analyzing: ", condition)
  obj <- prepared[[condition]]

  # LPEseq2: one pseudo-bulk sample per biological replicate.
  pb <- LPE_pseudobulk(
    seurat_obj = obj,
    method = "sample",
    sample_col = "sample",
    condition_col = "replicate",
    celltype_col = "celltype",
    min_cells = 10,
    normalize_by_ncells = TRUE,
    ncells_scale = 1000,
    normalize.method = "library_size",
    analysis.method = "LPE",
    n.bin = 100,
    df = 10,
    trim.method = "dvalue",
    d.threshold = 1.2,
    use_weighted_between = TRUE,
    p.method = "chisq",
    verbose = FALSE
  )

  lpe <- pb$results_by_celltype[["all"]]
  lpe_p <- data.frame(gene = lpe$gene, LPEseq2 = lpe$p.value)

  # MAST: cell-level test with cellular detection rate (CDR) adjustment.
  obj <- NormalizeData(obj, verbose = FALSE)
  Idents(obj) <- "replicate"
  obj$cngeneson <- as.numeric(scale(obj$nFeature_RNA))

  mast <- FindMarkers(
    obj,
    ident.1 = "Rep1",
    ident.2 = "Rep2",
    test.use = "MAST",
    latent.vars = "cngeneson",
    logfc.threshold = 0,
    min.pct = 0,
    verbose = FALSE
  )
  mast_p <- data.frame(gene = rownames(mast), MAST = mast$p_val)
  
  # Wilcoxon
  wilcox <- FindMarkers(
    obj, ident.1 = "Rep1", ident.2 = "Rep2", test.use = "wilcox",
    logfc.threshold = 0, min.pct = 0, verbose = FALSE
  )
  wilcox_p <- data.frame(gene = rownames(wilcox), Wilcoxon = wilcox$p_val)

  # Use the genes tested by all methods for the BH adjustment.
  p_values <- merge(merge(lpe_p, mast_p, by = "gene"), wilcox_p, by = "gene")
  p_values <- p_values[complete.cases(p_values), , drop = FALSE]

  write.csv(
    p_values,
    file.path(output_dir, paste0("pvalues_", condition, "_Rep1vRep2.csv")),
    row.names = FALSE
  )

  condition_summary <- data.frame(
    condition = condition,
    method = c("LPEseq2", "MAST", "Wilcoxon"),
    n_shared = nrow(p_values),
    common_BH = c(sum(p.adjust(p_values$LPEseq2,  "BH") < 0.05),
                  sum(p.adjust(p_values$MAST,     "BH") < 0.05),
                  sum(p.adjust(p_values$Wilcoxon, "BH") < 0.05))
  )
  condition_summary$pct_common_BH <- round(
    100 * condition_summary$common_BH / condition_summary$n_shared, 2
  )

  summary_results <- rbind(summary_results, condition_summary)
}

# Save the table for Supplementary Note S1 and Figure 2B.
write.csv(
  summary_results,
  file.path(output_dir, "negctrl_summary_commonBH.csv"),
  row.names = FALSE
)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
print(summary_results, row.names = FALSE)
message("Results saved in ", output_dir)
