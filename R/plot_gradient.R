# =========================================================
# scMetaTraj_plot_gradient
# =========================================================

#' Plot metabolic gradient on UMAP
#'
#' @param embedding UMAP coordinates (data.frame or matrix)
#' @param score Vector of metabolic axis score
#' @param title Plot title
#' @param palette Continuous color palette
#'
#' @return ggplot object
#' @export
scMetaTraj_plot_gradient <- function(
    embedding,
    score,
    title = "Metabolic gradient",
    palette
) {
  df <- data.frame(
    UMAP_1 = embedding[, 1],
    UMAP_2 = embedding[, 2],
    score  = score
  )
  
  ggplot(df, aes(UMAP_1, UMAP_2)) +
    geom_point(
      aes(color = score),
      size = 3.5,
      alpha = 0.85
    ) +
    scale_color_gradientn(colors = palette) +
    theme_classic() +
    labs(
      title = title,
      color = "Metabolic\nactivity"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5)
    )
}
