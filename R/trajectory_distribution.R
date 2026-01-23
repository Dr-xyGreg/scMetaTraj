#' Prepare mPT distribution data for metabolic subclusters
#'
#' @description
#' Organize metabolic pseudotime distributions by cluster and
#' automatically order clusters along median mPT.
#'
#' @param mPT Numeric vector of metabolic pseudotime.
#' @param cluster Factor or character vector of metabolic cluster labels.
#'
#' @return Data.frame with columns mPT and cluster (ordered factor).
#'
#' @export
scMetaTraj_mPT_distribution <- function(mPT, cluster) {
  df <- data.frame(
    mPT = mPT,
    cluster = as.factor(cluster)
  )
  
  ord <- stats::aggregate(mPT ~ cluster, df, stats::median, na.rm = TRUE)
  ord <- ord[order(ord$mPT), "cluster"]
  
  df$cluster <- factor(df$cluster, levels = as.character(ord))
  df
}
