scMetaTraj_plot_trend_by_cluster <- function(trend_by_cluster, palette) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required.")
  }
  library(ggplot2)
  
  ggplot(
    trend_by_cluster,
    aes(x = mPT_bin, y = score_smooth, color = cluster)
  ) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ module, scales = "free_y", ncol = 2) +
    scale_color_manual(values = palette) +
    theme_classic() +
    labs(
      x = "Metabolic pseudotime (mPT)",
      y = "Module score (smoothed)",
      color = "Metabolic cluster"
    )
}
