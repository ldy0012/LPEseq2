# Figure 2B: false positives in same-treatment null comparisons (SUM149PT, GSE142647)
# Input : S1_SUM149PT/results/negctrl_summary_commonBH.csv
#         (columns condition, method, n_shared, common_BH; methods LPEseq2, MAST, Wilcoxon)
# Output: S1_SUM149PT/results/Figure_2B.png / .pdf / .svg  (11 x 7 cm)
# Bar = genes with common-BH adjusted p < 0.05 in Rep1 vs Rep2 (all false positives);
# label = count (% of genes tested by all methods).

library(ggplot2)

result_dir <- file.path("S1_SUM149PT", "results")
fp <- read.csv(file.path(result_dir, "negctrl_summary_commonBH.csv"))

methods <- c("LPEseq2", "MAST", "Wilcoxon")
stopifnot(all(methods %in% fp$method))
fp <- fp[fp$method %in% methods, ]

fp$method <- factor(fp$method, levels = methods,
                    labels = c("LPEseq2", "MAST (CDR)", "Wilcoxon"))
fp$label <- sprintf("%s\n(%.1f%%)", format(fp$common_BH, big.mark = ",", trim = TRUE),
                    100 * fp$common_BH / fp$n_shared)
fp$panel <- sprintf("%s (%s genes)", fp$condition, format(fp$n_shared, big.mark = ",", trim = TRUE))
fp$panel <- factor(fp$panel, levels = unique(fp$panel[order(fp$condition != "DMSO")]))   # DMSO first

colours <- c("LPEseq2" = "#009E73", "MAST (CDR)" = "#8C6D4F", "Wilcoxon" = "#F4A582")

p <- ggplot(fp, aes(method, common_BH, fill = method)) +
  geom_col(width = 0.65, show.legend = FALSE) +
  geom_text(aes(label = label), vjust = -0.25, size = 7 / .pt, lineheight = 0.9) +
  facet_wrap(~ panel, nrow = 1) +
  scale_fill_manual(values = colours) +
  scale_y_continuous(labels = scales::label_comma(), expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = "False-positive genes") +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 9) +
  theme(axis.text = element_text(size = 8, colour = "black"),
        axis.text.x = element_text(angle = 30, hjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(size = 9, face = "bold"),
        panel.spacing = grid::unit(1, "lines"),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA))

for (ext in c("png", "pdf", "svg")) {
  ggsave(file.path(result_dir, paste0("Figure_2B.", ext)), p,
         width = 11, height = 7, units = "cm", dpi = 600, bg = "transparent",
         device = switch(ext, png = "png", pdf = grDevices::cairo_pdf, svg = svglite::svglite))
}
message("Saved Figure_2B in ", result_dir)
