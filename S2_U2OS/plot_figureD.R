# Figure D: medians across the 14 U2OS compounds (F1, selected genes, precision, recall)
# Input : S2_U2OS/results/metrics_summary.csv
#         (mean of the 10 subsamples per compound x method x n, written by run_analysis.R)
# Output: S2_U2OS/results/Figure_D.png / .pdf / .svg  (18 x 6 cm)
#         S2_U2OS/results/tableS3_median_over_compounds.csv (the plotted medians)
# Value = median over the 14 compounds of the per-compound mean (10 subsamples);
library(ggplot2)

result_dir <- file.path("S2_U2OS", "results")
summ <- read.csv(file.path(result_dir, "metrics_summary.csv"))

# median over compounds, per method x n
med <- aggregate(cbind(F1, calls, precision, recall) ~ method + n, data = summ, FUN = median)
med <- med[order(med$n, med$method), ]

table_S3 <- med
table_S3[c("F1", "precision", "recall")] <- round(table_S3[c("F1", "precision", "recall")], 3)
table_S3$calls <- round(table_S3$calls)
write.csv(table_S3, file.path(result_dir, "tableS3_median_over_compounds.csv"), row.names = FALSE)

# long format: one row per method x n x metric
metrics <- c(F1 = "F1", calls = "Selected genes", precision = "Precision", recall = "Recall")
long <- do.call(rbind, lapply(names(metrics), function(m)
  data.frame(method = med$method, n = med$n, metric = metrics[[m]], value = med[[m]])))
long$metric <- factor(long$metric, levels = metrics)
long$method <- factor(long$method, levels = c("LPEseq2", "edgeR"))
long$label  <- ifelse(long$metric == "Selected genes",
                      format(round(long$value), big.mark = ",", trim = TRUE),
                      sprintf("%.3f", long$value))

# label above the higher method and below the lower one at each n
higher <- ave(long$value, long$metric, long$n, FUN = max) == long$value
long$vjust <- ifelse(higher, -0.9, 1.9)

# y range per panel: 0-1 for proportions, 0-3,000 for selected genes
limits <- data.frame(metric = factor(rep(metrics, each = 2), levels = metrics),
                     value  = c(0, 1, 0, 3000, 0, 1, 0, 1), n = 1, method = "LPEseq2")
y_breaks <- function(lim) if (max(lim) > 10) c(0, 1500, 3000) else c(0, 0.5, 1)
y_labels <- function(b) if (max(b, na.rm = TRUE) > 10) format(b, big.mark = ",", trim = TRUE) else as.character(b)

colours <- c(LPEseq2 = "#1F6F9A", edgeR = "#C0622D")

p <- ggplot(long, aes(n, value, colour = method)) +
  geom_blank(data = limits) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.4) +
  geom_text(aes(label = label, vjust = vjust), size = 6 / .pt, show.legend = FALSE) +
  facet_wrap(~ metric, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = colours, name = NULL) +
  scale_x_continuous(breaks = 1:3, expand = expansion(add = 0.3)) +
  scale_y_continuous(breaks = y_breaks, labels = y_labels,
                     expand = expansion(mult = c(0.02, 0.12))) +
  labs(x = "Wells per dose × batch combination", y = NULL) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 9) +
  theme(axis.text = element_text(size = 8, colour = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size = 9, face = "bold"),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.3),
        panel.spacing = grid::unit(1.2, "lines"),
        legend.position = "top", legend.justification = "right",
        legend.text = element_text(size = 8),
        legend.margin = margin(0, 0, 0, 0),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        legend.background = element_rect(fill = "transparent", colour = NA))

for (ext in c("png", "pdf", "svg")) {
  ggsave(file.path(result_dir, paste0("Figure_D.", ext)), p,
         width = 18, height = 6, units = "cm", dpi = 600, bg = "transparent",
         device = switch(ext, png = "png", pdf = grDevices::cairo_pdf, svg = svglite::svglite))
}
print(table_S3, row.names = FALSE)
message("Saved Figure_D in ", result_dir)
