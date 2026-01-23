#' Compute module trends along mPT stratified by metabolic cluster
#'
#' @param score_mat Matrix/data.frame (cells x modules).
#' @param mPT Numeric vector.
#' @param cluster Factor/character vector of cluster labels.
#' @param modules Character vector of module names.
#' @param n_bins Integer. Number of bins.
#' @param smooth Logical. Whether to loess smooth.
#' @param span Numeric. Loess span.
#' @param min_cells Integer. Minimum cells per cluster to compute trends.
#'
#' @return Long-format data.frame with columns:
#'   cluster, module, mPT_bin, score, score_smooth, n_cells
#'
#' @export
scMetaTraj_trend_by_cluster <- function(
    score_mat,
    mPT,
    cluster,
    modules,
    n_bins = 30,
    smooth = TRUE,
    span = 0.3,
    min_cells = 50
) {
  score_mat <- as.data.frame(score_mat)
  cluster <- as.factor(cluster)
  
  if (nrow(score_mat) != length(mPT) || length(mPT) != length(cluster)) {
    stop("Lengths mismatch: score_mat rows, mPT, cluster must match.")
  }
  if (!all(modules %in% colnames(score_mat))) {
    miss <- modules[!modules %in% colnames(score_mat)]
    stop("Modules not found: ", paste(miss, collapse = ", "))
  }
  
  out <- list()
  for (cl in levels(cluster)) {
    idx <- which(cluster == cl)
    if (length(idx) < min_cells) next
    
    for (m in modules) {
      td <- scMetaTraj_trend(
        scores = score_mat[idx, m],
        mPT = mPT[idx],
        n_bins = n_bins,
        smooth = smooth,
        span = span
      )
      td$cluster <- cl
      td$module <- m
      td$n_cells <- length(idx)
      out[[length(out) + 1]] <- td
    }
  }
  
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}
#' Plot module trends along mPT stratified by cluster
#'
#' @param trend_by_cluster Output from scMetaTraj_trend_by_cluster().
#' @param palette Named vector of cluster colors (e.g., scMetaTraj_palette_discrete).
#'
#' @return ggplot object.
#' @export
scMetaTraj_plot_trend_by_cluster <- function(trend_by_cluster, palette) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required.")
  
  ggplot(trend_by_cluster, aes(x = mPT_bin, y = score_smooth, color = cluster)) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ module, scales = "free_y", ncol = 2) +
    ggplot2::scale_color_manual(values = palette) +
    theme_classic() +
    labs(x = "Metabolic pseudotime (mPT)", y = "Module score (smoothed)", color = "Cluster")
}
#' @keywords internal
NULL
