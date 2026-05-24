#' Plot multiple module trends with switchpoints
#'
#' @description
#' Visualize metabolic module trends along pseudotime with identified
#' switchpoints marked by vertical dashed lines.
#'
#' @param trend_long Data frame with columns: module, mPT_bin, score_smooth.
#'   Output from scMetaTraj_trend_multi()$trend_long.
#' @param switchpoints Data frame with columns: module, mPT_switch.
#'   Output from scMetaTraj_trend_multi()$switchpoints.
#'
#' @return A ggplot2 object showing trends faceted by module.
#'
#' @examples
#' # Create example trend data
#' trend_long <- data.frame(
#'   module = rep(c("Glycolysis", "OXPHOS"), each = 30),
#'   mPT_bin = rep(seq(0, 1, length.out = 30), 2),
#'   score_smooth = c(sin(seq(0, pi, length.out = 30)),
#'                    cos(seq(0, pi, length.out = 30)))
#' )
#' 
#' # Create example switchpoint data
#' switchpoints <- data.frame(
#'   module = c("Glycolysis", "OXPHOS"),
#'   mPT_switch = c(0.5, 0.3)
#' )
#' 
#' # Plot
#' p <- scMetaTraj_plot_trend_multi(trend_long, switchpoints)
#' print(p)
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_vline facet_wrap labs
#' @importFrom ggplot2 theme_classic theme_minimal
#' @export
scMetaTraj_plot_trend_multi <- function(trend_long, switchpoints) {
  
  ggplot2::ggplot(
    trend_long,
    ggplot2::aes(
      x = .data$mPT_bin,
      y = .data$score_smooth
    )
  ) +
    ggplot2::geom_line(linewidth = 1.0, color = "black") +
    ggplot2::geom_vline(
      data = switchpoints,
      ggplot2::aes(xintercept = .data$mPT_switch),
      linetype = "dashed",
      linewidth = 0.4
    ) +
    ggplot2::facet_wrap(~ module, scales = "free_y", ncol = 2) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = "Metabolic pseudotime (mPT)",
      y = "Module score (smoothed)"
    )
}
