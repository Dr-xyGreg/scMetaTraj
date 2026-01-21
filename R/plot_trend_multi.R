#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_vline facet_wrap labs
#' @importFrom ggplot2 theme_classic theme_minimal

scMetaTraj_plot_trend_multi <- function(trend_long, switchpoints) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required.")
  }
  library(ggplot2)
  
  ggplot(trend_long, aes(x = mPT_bin, y = score_smooth)) +
    geom_line(linewidth = 1.1, color = "black") +
    geom_vline(
      data = switchpoints,
      aes(xintercept = mPT_switch),
      linetype = "dashed",
      linewidth = 0.5,
      color = "red"
    ) +
    facet_wrap(~ module, scales = "free_y", ncol = 2) +
    theme_classic() +
    labs(
      x = "Metabolic pseudotime (mPT)",
      y = "Module score (smoothed)"
    )
}
