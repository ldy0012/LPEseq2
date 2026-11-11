# Figure 2C: per-compound F1, LPEseq2 vs edgeR (U2OS DRUG-seq, GSE176150)
# Plot only (no analysis). Run from the repository root: Rscript S2_U2OS/plot_figure2C.R
#
# Input : S2_U2OS/results/metrics_by_run.csv     F1 of every subsample run
#         S2_U2OS/results/reference_summary.csv  reference-positive % per compound
# Output: S2_U2OS/results/Figure_2C.png / .pdf / .svg  (11 x 9.5 cm)
#
# Point = mean F1 over the 10 subsamples (one compound x wells per condition);
# bars  = SD over the same 10 subsamples; colour = reference-positive genes (%).

library(ggplot2)

result_dir <- file.path("S2_U2OS", "results")
runs <- read.csv(file.path(result_dir, "metrics_by_run.csv"))
reference <- read.csv(file.path(result_dir, "reference_summary.csv"))

# mean and SD of F1 over the 10 subsamples, per compound x method x n
f1_mean <- aggregate(F1 ~ compound + method + n, data = runs, FUN = mean)
f1_sd   <- aggregate(F1 ~ compound + method + n, data = runs, FUN = sd)
f1 <- merge(f1_mean, f1_sd, by = c("compound", "method", "n"), suffixes = c("_mean", "_sd"))

plot_data <- merge(
  subset(f1, method == "edgeR",   c(compound, n, F1_mean, F1_sd)),
  subset(f1, method == "LPEseq2", c(compound, n, F1_mean, F1_sd)),
  by = c("compound", "n"), suffixes = c("_edgeR", "_LPE")
)
plot_data$pct   <- reference$ref_pct[match(plot_data$compound, reference$compound)]
plot_data$wells <- factor(plot_data$n, levels = 1:3)
plot_data$label <- ifelse(plot_data$n == 3, plot_data$compound, "")   # names at 3 wells

clamp <- function(v) pmin(1, pmax(0, v))

p <- ggplot(plot_data, aes(F1_mean_edgeR, F1_mean_LPE, colour = pct, shape = wells)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(ymin = clamp(F1_mean_LPE - F1_sd_LPE),
                    ymax = clamp(F1_mean_LPE + F1_sd_LPE)),
                width = 0, linewidth = 0.35, alpha = 0.7) +
  geom_errorbar(aes(xmin = clamp(F1_mean_edgeR - F1_sd_edgeR),
                    xmax = clamp(F1_mean_edgeR + F1_sd_edgeR)),
                width = 0, linewidth = 0.35, alpha = 0.7, orientation = "y") +
  geom_point(size = 2) +
  ggrepel::geom_text_repel(aes(label = label), colour = "grey25", size = 6 / .pt,
                           show.legend = FALSE, max.overlaps = Inf, min.segment.length = 0,
                           segment.size = 0.25, segment.colour = "grey45",
                           box.padding = 0.35, point.padding = 0.25, force = 2, seed = 2026) +
  scale_colour_viridis_c(name = "Reference\nDEGs (%)") +
  scale_shape_manual(values = c(16, 17, 15), name = "Wells per\ncondition") +
  scale_x_continuous(breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
  labs(x = "edgeR F1 (reference: full-data edgeR)",
       y = "LPEseq2 F1 (reference: full-data edgeR)") +
  guides(shape = guide_legend(order = 1),
         colour = guide_colourbar(order = 2, barheight = grid::unit(2.2, "cm"),
                                  barwidth = grid::unit(0.3, "cm"))) +
  theme_classic(base_size = 9) +
  theme(axis.text = element_text(size = 8, colour = "black"),
        legend.text = element_text(size = 8),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA),
        legend.key = element_rect(fill = "transparent", colour = NA),
        plot.margin = margin(6, 8, 6, 6))

for (ext in c("png", "pdf", "svg")) {
  ggsave(file.path(result_dir, paste0("Figure_2C.", ext)), p,
         width = 11, height = 9.5, units = "cm", dpi = 600, bg = "transparent",
         device = switch(ext, png = "png", pdf = grDevices::cairo_pdf, svg = svglite::svglite))
}
message("Saved Figure_2C in ", result_dir)
