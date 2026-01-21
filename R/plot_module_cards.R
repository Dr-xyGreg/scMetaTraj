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
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  df_long <- cluster_profile %>%
    as.data.frame() %>%
    mutate(cluster = rownames(cluster_profile)) %>%
    pivot_longer(
      cols = -cluster,
      names_to = "module",
      values_to = "score"
    ) %>%
    left_join(module_axis_map, by = "module") %>%
    group_by(module) %>%
    arrange(desc(score), .by_group = TRUE) %>%
    mutate(cluster_rank = factor(cluster, levels = cluster))
  
  ggplot(df_long, aes(cluster_rank, score, group = module)) +
    geom_rect(
      aes(
        xmin = -Inf, xmax = Inf,
        ymin = -Inf, ymax = Inf,
        fill = axis
      ),
      alpha = 0.12,
      inherit.aes = FALSE
    ) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey70", linewidth = 0.3) +
    geom_line(color = "white", linewidth = 0.8) +
    geom_point(color = "white", size = 2) +
    facet_wrap(~ module, scales = "free_x", ncol = 4) +
    scale_fill_manual(values = axis_palette) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      strip.text = element_text(face = "bold"),
      legend.position = "none"
    ) +
    labs(
      x = "Metabolic clusters (ranked within module)",
      y = "Module activity (Z-score)"
    )
}
