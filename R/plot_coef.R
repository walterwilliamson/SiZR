#' Plot Coefficient Estimates with Confidence Bands
#'
#' @param coef_df Data frame containing coefficient estimates and confidence bands.
#' @param f_X True coefficient function for comparison.
#'
#' @return A ggplot object showing the coefficient estimates, confidence bands, and true function.
#' @import ggplot2
#' @export
#'
#' @examples
#' # Plot data


plot_coef <- function(coef_df, f_X) {

  a <- data.frame(x = coef_df$xind, y = f_X(coef_df$xind))

  ggplot() +
    geom_line(data = coef_df, aes(x = .data$xind, y = .data$estimate)) +
    geom_line(data = coef_df, aes(x = .data$xind, y = .data$UB, group = .data$coefficient),
              color = "blue", linetype = "dashed") +
    geom_line(data = coef_df, aes(x = .data$xind,y = .data$LB, group = .data$coefficient),
              color = "blue", linetype = "dashed") +
    geom_line(data = a, aes(x = .data$x, y = .data$y), color = "red") +
    facet_wrap(~.data$coefficient, scales = "free_y") +
    theme_classic()
}
