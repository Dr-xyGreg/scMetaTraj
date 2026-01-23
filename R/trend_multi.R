#' Compute trends and switchpoints for multiple modules along mPT
#'
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
  if (nrow(score_mat) != length(mPT)) {
    stop("nrow(score_mat) must equal length(mPT).")
  }
  
  trend_list <- base::lapply(modules, function(m) {
    td <- scMetaTraj_trend(
      scores = score_mat[[m]],
      mPT = mPT,
      n_bins = n_bins,
      smooth = smooth,
      span = span
    )
    td$module <- m
    td
  })
  
  trend_long <- base::do.call(rbind, trend_list)
  rownames(trend_long) <- NULL
  
  sp <- base::lapply(modules, function(m) {
    td <- trend_long[trend_long$module == m, , drop = FALSE]
    sw <- scMetaTraj_switchpoint(td)
    data.frame(
      module = m,
      mPT_switch = sw$mPT_switch,
      index = sw$index
    )
  })
  
  switchpoints <- base::do.call(rbind, sp)
  rownames(switchpoints) <- NULL
  
  list(
    trend_long = trend_long,
    switchpoints = switchpoints
  )
}
