# Supplementary Note S2 / Figure 2C
# U2OS DRUG-seq (GSE176150): replicate-subsampling benchmark

library(edgeR)
library(LPEseq2)

input_file <- file.path("S2_U2OS", "data", "u2os_counts_meta.rds")
output_dir <- file.path("S2_U2OS", "results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Load the merged count matrix and sample metadata.
data <- readRDS(input_file)
counts <- data$counts
meta <- data$meta

# Exclude untreated controls and samples with missing compound names.
keep_samples <- !is.na(meta$compound) & meta$compound != "DMSO"
meta <- meta[keep_samples, , drop = FALSE]

# 2. Statistical tests: dose effect adjusted for batch.
run_edger <- function(gene_counts, sample_info) {
  design <- model.matrix(~ dose + batch, data = sample_info)

  y <- DGEList(counts = gene_counts)
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)

  dose_columns <- grep("^dose", colnames(design))
  test <- glmQLFTest(fit, coef = dose_columns)
  return(test$table$PValue)
}

run_lpe <- function(gene_counts, sample_info) {
  processed <- LPE_preprocess(
    gene_counts, sample_info,
    design = ~ dose + batch,
    normalize.method = "TMM",
    min.count = 0,
    verbose = FALSE
  )

  result <- LPE_ANOVA(
    object = processed,
    analysis.method = "LPE",
    n.bin = 100,
    df = 10,
    trim.method = "dvalue",
    d.threshold = 1.2,
    use_weighted_between = TRUE,
    variance.eval = "grand_mean",
    p.method = "chisq",
    verbose = FALSE
  )

  # Keep the p-values in the same gene order as gene_counts.
  gene_positions <- match(rownames(gene_counts), result$gene)
  return(result$p.value_dose[gene_positions])
}

# 3. Compare significant calls with the full-data edgeR reference.
score_results <- function(p_values, reference) {
  p_values[is.na(p_values)] <- 1
  significant <- p.adjust(p_values, method = "BH") < 0.05

  true_positives <- sum(significant & reference)
  total_calls <- sum(significant)
  total_reference <- sum(reference)

  precision <- if (total_calls > 0) true_positives / total_calls else NA_real_
  recall <- true_positives / total_reference
  F1 <- if (true_positives == 0) 0 else 2 * precision * recall / (precision + recall)

  return(data.frame(calls = total_calls, precision = precision,
                    recall = recall, F1 = F1))
}

# 4. Select compounds with 8 doses, 3 batches, and >=6 wells in every cell.
compound_names <- sort(unique(meta$compound))
compounds <- character(0)

for (compound in compound_names) {
  compound_meta <- meta[meta$compound == compound, , drop = FALSE]
  well_table <- table(as.character(compound_meta$dose),
                      as.character(compound_meta$batch))

  if (nrow(well_table) == 8 && ncol(well_table) == 3 &&
      all(well_table >= 6)) {
    compounds <- c(compounds, compound)
  }
}

# 5. Full-data reference and 10 matched subsamples for each sample size.
reference_summary <- data.frame()
metrics_by_run <- data.frame()

for (compound in compounds) {
  message("Analyzing ", compound)

  sample_info <- droplevels(meta[meta$compound == compound, c("dose", "batch"), drop = FALSE])
  sample_info$batch <- factor(as.character(sample_info$batch))
  gene_counts <- counts[, rownames(sample_info), drop = FALSE]

  design <- model.matrix(~ dose + batch, data = sample_info)
  keep_genes <- filterByExpr(DGEList(gene_counts), design = design)
  gene_counts <- gene_counts[keep_genes, , drop = FALSE]

  reference_p <- run_edger(gene_counts, sample_info)
  reference <- p.adjust(reference_p, method = "BH") < 0.05

  reference_summary <- rbind(reference_summary, data.frame(
    compound = compound,
    wells = ncol(gene_counts),
    genes = nrow(gene_counts),
    ref_DEG = sum(reference),
    ref_pct = 100 * mean(reference)
  ))

  # Each dose-by-batch combination is sampled independently.
  well_groups <- interaction(sample_info$dose, sample_info$batch, drop = TRUE)
  wells_by_group <- split(rownames(sample_info), well_groups)

  for (n in 1:3) {
    for (repeat_id in 1:10) {
      set.seed(5000 * n + repeat_id)
      selected_wells <- character(0)

      for (group_wells in wells_by_group) {
        selected_wells <- c(selected_wells, sample(group_wells, size = n))
      }

      sub_counts <- gene_counts[, selected_wells, drop = FALSE]
      sub_info <- sample_info[selected_wells, , drop = FALSE]

      # Both methods use exactly the same selected wells.
      for (method in c("LPEseq2", "edgeR")) {
        if (method == "LPEseq2") {
          p_values <- run_lpe(sub_counts, sub_info)
        } else {
          p_values <- run_edger(sub_counts, sub_info)
        }

        metrics <- score_results(p_values, reference)
        metrics_by_run <- rbind(metrics_by_run, data.frame(
          compound = compound,
          method = method,
          n = n,
          rep = repeat_id,
          calls = metrics$calls,
          precision = metrics$precision,
          recall = metrics$recall,
          F1 = metrics$F1
        ))
      }
    }
  }
}

write.csv(reference_summary, file.path(output_dir, "reference_summary.csv"),
          row.names = FALSE)
write.csv(metrics_by_run, file.path(output_dir, "metrics_by_run.csv"),
          row.names = FALSE)

# 6. Average the 10 subsampling runs for each compound, method, and n.
metrics_summary <- data.frame()

for (compound in compounds) {
  for (method in c("LPEseq2", "edgeR")) {
    for (n in 1:3) {
      selected <- metrics_by_run[
        metrics_by_run$compound == compound &
          metrics_by_run$method == method &
          metrics_by_run$n == n, , drop = FALSE
      ]

      metrics_summary <- rbind(metrics_summary, data.frame(
        compound = compound, method = method, n = n,
        calls = mean(selected$calls),
        precision = mean(selected$precision),
        recall = mean(selected$recall),
        F1 = mean(selected$F1)
      ))
    }
  }
}

metrics_summary <- metrics_summary[
  order(metrics_summary$compound, metrics_summary$method, metrics_summary$n),
  , drop = FALSE
]

write.csv(metrics_summary, file.path(output_dir, "metrics_summary.csv"),
          row.names = FALSE)

# 7. Supplementary Table S3: median of the per-compound means.
table_S3 <- data.frame()

for (method in c("LPEseq2", "edgeR")) {
  for (n in 1:3) {
    selected <- metrics_summary[
      metrics_summary$method == method & metrics_summary$n == n, , drop = FALSE
    ]

    table_S3 <- rbind(table_S3, data.frame(
      method = method, n = n,
      F1 = median(selected$F1),
      calls = median(selected$calls),
      precision = median(selected$precision),
      recall = median(selected$recall)
    ))
  }
}

table_S3 <- table_S3[order(table_S3$n, table_S3$method), , drop = FALSE]
table_S3$F1 <- round(table_S3$F1, 3)
table_S3$precision <- round(table_S3$precision, 3)
table_S3$recall <- round(table_S3$recall, 3)
table_S3$calls <- round(table_S3$calls)

write.csv(table_S3, file.path(output_dir, "tableS3_median_over_compounds.csv"),
          row.names = FALSE)
print(table_S3, row.names = FALSE)

# 8. Figure 2C: paired mean and SD of F1 (edgeR vs LPEseq2).
figure_2C <- data.frame()

for (compound in compounds) {
  for (n in 1:3) {
    edge <- metrics_by_run[
      metrics_by_run$compound == compound & metrics_by_run$method == "edgeR" &
        metrics_by_run$n == n, , drop = FALSE
    ]
    lpe <- metrics_by_run[
      metrics_by_run$compound == compound & metrics_by_run$method == "LPEseq2" &
        metrics_by_run$n == n, , drop = FALSE
    ]

    ref_percent <- reference_summary$ref_pct[
      reference_summary$compound == compound
    ]

    figure_2C <- rbind(figure_2C, data.frame(
      compound = compound, n = n,
      mean_edgeR = mean(edge$F1), sd_edgeR = sd(edge$F1),
      mean_LPE = mean(lpe$F1), sd_LPE = sd(lpe$F1),
      pct = ref_percent
    ))
  }
}

write.csv(figure_2C, file.path(output_dir, "paired_F1_mean_SD_edgeR.csv"),
          row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("Results saved in ", output_dir)
