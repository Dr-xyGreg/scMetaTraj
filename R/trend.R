#' Compute metabolic module trend along pseudotime
#'
#' @description
#' Bins cells along metabolic pseudotime (mPT) and computes mean module 
#' scores per bin, with optional loess smoothing.
#'
#' @param scores Numeric vector of module scores (length = n_cells).
#' @param mPT Numeric vector of metabolic pseudotime values (length = n_cells).
#' @param n_bins Integer. Number of bins along mPT for trend computation.
#' @param smooth Logical. Whether to apply loess smoothing to the binned trend.
#' @param span Numeric. Loess span parameter (only used if smooth = TRUE).
#'
#' @return A data frame with columns:
#'   \item{mPT_bin}{Mid-point of each mPT bin}
#'   \item{score}{Mean score per bin}
#'   \item{score_smooth}{Smoothed score (if smooth = TRUE, otherwise same as score)}
#'
#' @examples
#' # Create example data
#' set.seed(123)
#' n_cells <- 200
#' mPT <- runif(n_cells, 0, 1)
#' scores <- sin(mPT * 2 * pi) + rnorm(n_cells, 0, 0.1)
#' 
#' # Compute trend
#' trend <- scMetaTraj_trend(
#'   scores = scores,
#'   mPT = mPT,
#'   n_bins = 30,
#'   smooth = TRUE,
#'   span = 0.3
#' )
#' 
#' # Plot trend
#' plot(trend$mPT_bin, trend$score_smooth, type = "l",
#'      xlab = "Metabolic pseudotime", ylab = "Module score")
#'
#' @export
scMetaTraj_trend <- function(
    scores,
    mPT,
    n_bins = 30,
    smooth = TRUE,
    span = 0.3
) {
  
  # Input validation
  if (length(scores) != length(mPT)) {
    stop("scores and mPT must have the same length.")
  }
  
  if (any(is.na(scores)) || any(is.na(mPT))) {
    stop("scores and mPT cannot contain NA values.")
  }
  
  # Create bins
  mPT_breaks <- seq(min(mPT), max(mPT), length.out = n_bins + 1)
  mPT_bins <- cut(mPT, breaks = mPT_breaks, include.lowest = TRUE)
  
  # Compute mean score per bin
  score_mean <- tapply(scores, mPT_bins, mean, na.rm = TRUE)
  
  # Get bin midpoints
  mPT_bin <- (mPT_breaks[-length(mPT_breaks)] + mPT_breaks[-1]) / 2
  
  # Create result data frame
  result <- data.frame(
    mPT_bin = mPT_bin,
    score = as.numeric(score_mean)
  )
  
  # Remove bins with no data
  result <- result[!is.na(result$score), ]
  
  # Apply loess smoothing if requested
  if (smooth && nrow(result) > 3) {
    tryCatch({
      loess_fit <- stats::loess(
        score ~ mPT_bin,
        data = result,
        span = span
      )
      result$score_smooth <- stats::predict(loess_fit, newdata = result$mPT_bin)
    }, error = function(e) {
      warning("Loess smoothing failed, using raw scores.")
      result$score_smooth <- result$score
    })
  } else {
    result$score_smooth <- result$score
  }
  
  return(result)
}


#' Identify metabolic trajectory switchpoint
#'
#' @description
#' Identifies the point along metabolic pseudotime where a module shows
#' maximum change in trend (inflection point).
#'
#' @param trend_df Data frame with columns: mPT_bin and score_smooth.
#'   Typically output from \code{\link{scMetaTraj_trend}}.
#'
#' @return A list with:
#'   \item{mPT_switch}{Numeric. The mPT value at the switchpoint}
#'   \item{index}{Integer. The index (row number) of the switchpoint in trend_df}
#'
#' @examples
#' # Create example trend data
#' set.seed(456)
#' n_cells <- 200
#' mPT <- runif(n_cells, 0, 1)
#' 
#' # Simulate trend with switchpoint at mPT = 0.5
#' scores <- ifelse(mPT < 0.5, 
#'                  0.3 + rnorm(n_cells, 0, 0.05),
#'                  0.7 + rnorm(n_cells, 0, 0.05))
#' 
#' # Compute trend
#' trend <- scMetaTraj_trend(scores, mPT, n_bins = 30, smooth = TRUE)
#' 
#' # Find switchpoint
#' switchpoint <- scMetaTraj_switchpoint(trend)
#' print(switchpoint$mPT_switch)
#' 
#' # Visualize
#' plot(trend$mPT_bin, trend$score_smooth, type = "l",
#'      xlab = "Metabolic pseudotime", ylab = "Module score")
#' abline(v = switchpoint$mPT_switch, col = "red", lty = 2)
#'
#' @export
scMetaTraj_switchpoint <- function(trend_df) {
  
  # Input validation
  if (!is.data.frame(trend_df)) {
    stop("trend_df must be a data frame.")
  }
  
  required_cols <- c("mPT_bin", "score_smooth")
  if (!all(required_cols %in% colnames(trend_df))) {
    stop("trend_df must contain columns: mPT_bin and score_smooth.")
  }
  
  if (nrow(trend_df) < 3) {
    warning("Too few points to identify switchpoint. Returning midpoint.")
    mid_idx <- ceiling(nrow(trend_df) / 2)
    return(list(
      mPT_switch = trend_df$mPT_bin[mid_idx],
      index = mid_idx
    ))
  }
  
  # Compute first derivative (rate of change)
  score_diff <- diff(trend_df$score_smooth)
  
  # Compute second derivative (change in rate of change)
  # Switchpoint = maximum absolute second derivative
  if (length(score_diff) > 1) {
    score_diff2 <- diff(score_diff)
    
    # Find index of maximum absolute change
    switch_idx <- which.max(abs(score_diff2)) + 1  # +1 to account for diff
    
  } else {
    # Fallback: use maximum absolute first derivative
    switch_idx <- which.max(abs(score_diff)) + 1
  }
  
  # Ensure valid index
  if (switch_idx < 1) switch_idx <- 1
  if (switch_idx > nrow(trend_df)) switch_idx <- nrow(trend_df)
  
  return(list(
    mPT_switch = trend_df$mPT_bin[switch_idx],
    index = switch_idx
  ))
}
