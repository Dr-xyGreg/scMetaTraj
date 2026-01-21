#' Compute trends and switchpoints for multiple modules along mPT
#'#' @importFrom stats aggregate median dist

#' @param score_mat Matrix/data.frame (cells x modules). Row order must match mPT.
#' @param mPT Numeric vector (length = n_cells).
#' @param modules Character vector of module names (columns of score_mat).
#' @param n_bins Integer. Number of mPT bins.
#' @param smooth Logical. Whether to loess smooth.
#' @param span Numeric. Loess span.
#'
#' @return A list with:
#'   - trend_long: long-format data.frame for plotting
#'   - switchpoints: data.frame of module-wise switchpoints
#'
#' @export
scMetaTraj_trend_multi <- function(
    score_mat,
    mPT,
    modules,
    n_bins = 30,
    smooth = TRUE,
    span = 0.3
) {
  score_mat <- as.data.frame(score_mat)
  
  if (!all(modules %in% colnames(score_mat))) {
    miss <- modules[!modules %in% colnames(score_mat)]
    stop("Modules not found in score_mat: ", paste(miss, collapse = ", "))
  }
  if (nrow(score_mat) != length(mPT)) stop("nrow(score_mat) must equal length(mPT).")
  
  trend_list <- lapply(modules, function(m) {
    td <- scMetaTraj_trend(scores = score_mat[[m]], mPT = mPT, n_bins = n_bins, smooth = smooth, span = span)
    td$module <- m
    td
  })
  
  trend_long <- do.call(rbind, trend_list)
  rownames(trend_long) <- NULL
  
  # switchpoint per module
  sp <- lapply(modules, function(m) {
    td <- trend_long[trend_long$module == m, , drop = FALSE]
    sw <- scMetaTraj_switchpoint(td)
    data.frame(module = m, mPT_switch = sw$mPT_switch, index = sw$index)
  })
  switchpoints <- do.call(rbind, sp)
  rownames(switchpoints) <- NULL
  
  list(trend_long = trend_long, switchpoints = switchpoints)
}
#' Plot module trends along mPT with switchpoints
#'
#' @param trend_long Output $trend_long from scMetaTraj_trend_multi().
#' @param switchpoints Output $switchpoints from scMetaTraj_trend_multi().
#'
#' @return ggplot object.
#' @export
scMetaTraj_plot_trend_multi <- function(trend_long, switchpoints) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required.")
  
  library(ggplot2)
  
  ggplot(trend_long, aes(x = mPT_bin, y = score_smooth)) +
    geom_line(linewidth = 1.0, color = "black") +
    geom_vline(
      data = switchpoints,
      aes(xintercept = mPT_switch),
      linetype = "dashed",
      linewidth = 0.4
    ) +
    facet_wrap(~ module, scales = "free_y", ncol = 2) +
    theme_classic() +
    labs(x = "Metabolic pseudotime (mPT)", y = "Module score (smoothed)")
}
