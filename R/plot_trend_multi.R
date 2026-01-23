#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_vline facet_wrap labs
#' @importFrom ggplot2 theme_classic theme_minimal
#' @importFrom scMetaTraj scMetaTraj_trend scMetaTraj_switchpoint

scMetaTraj_plot_trend_multi <- function(trend_long, switchpoints) {
  
  ggplot2::ggplot(
    trend_long,
    ggplot2::aes(
      x = .data$mPT_bin,
      y = .data$score_smooth
    )
  ) +
    ggplot2::geom_line(linewidth = 1.0, color = "black") +
    ggplot2::geom_vline(
      data = switchpoints,
      ggplot2::aes(xintercept = .data$mPT_switch),
      linetype = "dashed",
      linewidth = 0.4
    ) +
    ggplot2::facet_wrap(~ module, scales = "free_y", ncol = 2) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = "Metabolic pseudotime (mPT)",
      y = "Module score (smoothed)"
    )
}
#' @keywords internal
NULL
