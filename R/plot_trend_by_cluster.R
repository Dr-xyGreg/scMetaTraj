#' Plot module trends by metabolic cluster
#'
#' @param trend_by_cluster Data.frame of trends by cluster.
#' @param palette Named vector of colors for clusters.
#'
#' @return ggplot object.
#' @export
scMetaTraj_plot_trend_by_cluster <- function(trend_by_cluster, palette) {
  
  ggplot2::ggplot(
    trend_by_cluster,
    ggplot2::aes(
      x = .data$mPT_bin,
      y = .data$score_smooth,
      color = .data$cluster
    )
  ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::facet_wrap(
      ~ module,
      scales = "free_y",
      ncol = 2
    ) +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = "Metabolic pseudotime (mPT)",
      y = "Module score (smoothed)",
      color = "Metabolic cluster"
    )
}

