#' Identify candidate transition zone along mPT
#'
#' @description
#' Identify a candidate transition zone along mPT where the composition
#' of metabolic subclusters changes most rapidly. Intended as a
#' hypothesis-generating indicator.
#'
#' @param mPT Numeric vector of metabolic pseudotime.
#' @param cluster Metabolic cluster labels.
#' @param n_bins Integer. Number of bins along mPT.
#' @param top_frac Fraction of bins with highest composition change to define zone.
#'
#' @return Named numeric vector with xmin and xmax.
#'
#' @export
scMetaTraj_transition_zone <- function(
    mPT,
    cluster,
    n_bins = 30,
    top_frac = 0.2
) {
  df <- data.frame(
    mPT = mPT,
    cluster = as.factor(cluster)
  )
  
  bins <- cut(df$mPT, breaks = seq(0, 1, length.out = n_bins + 1),
              include.lowest = TRUE)
  
  tab <- table(bins, df$cluster)
  prop <- prop.table(tab, margin = 1)
  
  # L1 distance between consecutive bins
  delta <- apply(prop[-1, ] - prop[-nrow(prop), ], 1,
                 function(x) sum(abs(x)))
  
  idx <- order(delta, decreasing = TRUE)[
    seq_len(max(1, ceiling(length(delta) * top_frac)))
  ]
  
  bin_levels <- levels(bins)
  left_edges <- seq(0, 1, length.out = n_bins + 1)[-length(seq(0, 1, length.out = n_bins + 1))]
  
  xmin <- min(left_edges[idx])
  xmax <- max(left_edges[idx] + diff(left_edges)[1])
  
  c(xmin = xmin, xmax = xmax)
}
