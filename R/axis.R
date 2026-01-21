# =========================================================
# scMetaTraj_axis_score
# Construct metabolic axis scores from module scores
# =========================================================

#' Calculate metabolic axis scores
#'
#' @param scores Cell-by-module metabolic score matrix
#' @param axis_list Named list defining metabolic axes
#' @param scale Logical, whether to scale axis scores (Z-score)
#'
#' @return A matrix of cell-by-axis scores
#' @export
scMetaTraj_axis_score <- function(
    scores,
    axis_list,
    scale = TRUE
) {
  stopifnot(is.matrix(scores) || is.data.frame(scores))
  
  axis_scores <- sapply(names(axis_list), function(ax) {
    modules <- axis_list[[ax]]
    modules <- intersect(modules, colnames(scores))
    if (length(modules) == 0) {
      rep(NA, nrow(scores))
    } else {
      rowMeans(scores[, modules, drop = FALSE], na.rm = TRUE)
    }
  })
  
  axis_scores <- as.matrix(axis_scores)
  rownames(axis_scores) <- rownames(scores)
  
  if (scale) {
    axis_scores <- scale(axis_scores)
  }
  
  axis_scores
}
