#' Summarize metabolic profiles of clusters
#'
#' @description
#' scMetaTraj_cluster_profile() computes representative metabolic
#' pathway activities for each metabolic cluster.
#'
#' @param scores Numeric matrix, cells x pathways.
#' @param metabolic_cluster Factor or character vector, length = nrow(scores).
#' @param stat Character. "median" (default) or "mean".
#' @param scale Logical. Whether to z-score pathways across clusters.
#'
#' @return A data.frame: clusters x pathways.
#'
scMetaTraj_cluster_profile <- function(
    scores,
    metabolic_cluster,
    stat = c("median", "mean"),
    scale = TRUE
) {
  
  stat <- match.arg(stat)
  
  if (!is.matrix(scores)) {
    stop("scores must be a numeric matrix (cells x pathways).")
  }
  
  if (length(metabolic_cluster) != nrow(scores)) {
    stop("Length of metabolic_cluster must equal number of rows in scores.")
  }
  
  metabolic_cluster <- as.factor(metabolic_cluster)
  clusters <- levels(metabolic_cluster)
  
  prof <- base::sapply(clusters, function(cl) {
    idx <- which(metabolic_cluster == cl)
    sub <- scores[idx, , drop = FALSE]
    
    if (stat == "median") {
      base::apply(sub, 2, stats::median, na.rm = TRUE)
    } else {
      base::colMeans(sub, na.rm = TRUE)
    }
  })
  
  prof <- base::t(prof)
  rownames(prof) <- clusters
  
  if (scale) {
    prof <- base::scale(prof)
  }
  
  as.data.frame(prof)
}
