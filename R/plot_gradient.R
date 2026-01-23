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
  
  ggplot2::ggplot(
    df,
    ggplot2::aes(x = UMAP_1, y = UMAP_2)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(color = score),
      size = 3.5,
      alpha = 0.85
    ) +
    ggplot2::scale_color_gradientn(colors = palette) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = title,
      color = "Metabolic\nactivity"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
}
