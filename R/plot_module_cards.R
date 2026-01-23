# =========================================================
# scMetaTraj_plot_module_cards
# =========================================================

#' Plot metabolic module cards
#'
#' @param cluster_profile Cluster-by-module matrix
#' @param module_axis_map Data.frame with columns: module, axis
#' @param axis_palette Named vector of colors for axes
#'
#' @return ggplot object
#' @export
scMetaTraj_plot_module_cards <- function(
    cluster_profile,
    module_axis_map,
    axis_palette
) {
  
  df_long <- cluster_profile |>
    as.data.frame() |>
    dplyr::mutate(cluster = rownames(cluster_profile)) |>
    tidyr::pivot_longer(
      cols = -cluster,
      names_to = "module",
      values_to = "score"
    ) |>
    dplyr::left_join(module_axis_map, by = "module") |>
    dplyr::group_by(.data$module) |>
    dplyr::arrange(dplyr::desc(.data$score), .by_group = TRUE) |>
    dplyr::mutate(
      cluster_rank = factor(.data$cluster, levels = .data$cluster)
    )
  
  ggplot2::ggplot(
    df_long,
    ggplot2::aes(
      x = .data$cluster_rank,
      y = .data$score,
      group = .data$module
    )
  ) +
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = -Inf,
        xmax = Inf,
        ymin = -Inf,
        ymax = Inf,
        fill = .data$axis
      ),
      alpha = 0.12,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey70",
      linewidth = 0.3
    ) +
    ggplot2::geom_line(color = "white", linewidth = 0.8) +
    ggplot2::geom_point(color = "white", size = 2) +
    ggplot2::facet_wrap(~ module, scales = "free_x", ncol = 4) +
    ggplot2::scale_fill_manual(values = axis_palette) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        size = 8
      ),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "none"
    ) +
    ggplot2::labs(
      x = "Metabolic clusters (ranked within module)",
      y = "Module activity (Z-score)"
    )
}
