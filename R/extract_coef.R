#' Extract Coefficients and Confidence Intervals from mgcv Fits
#'
#' @param fit A fitted model object from mgcv::gam or mgcv::bam.
#' @param xind A numeric vector of values at which to evaluate the smooth terms. Default is seq(0, 1, len = 20).
#' @param conf_level Confidence level for the intervals. Default is 0.95.
#' @param type Type of model fitted. Options are "SiZR" (default) or "mean".
#'
#' @return A data frame containing the estimates, standard errors, confidence intervals, and corresponding x values.
#' @import dplyr
#' @importFrom stats qnorm
#' @export
#'
#' @examples
#' # Assuming fit_SiZR is a fitted model from the SiZR approach
#' # coef_df <- extract_coef(fit_SiZR, xind = seq(0, 1, len = 100), conf_level = 0.95, type = "SiZR")



extract_coef <- function(fit, xind = seq(0, 1, len = 20), conf_level = 0.95,
                                 type = "SiZR") {
  # Build prediction grid
  if(type == "SiZR"){
    xmat_pred <- data.frame(X = xind, Lmat = 1)
  } else if (type == "mean"){
    vars <- names(fit$model)[-1]  # drop response var
    xmat_pred <- expand.grid(rep(list(xind), length(vars)))
    names(xmat_pred) <- vars
  }

  # Predict terms
  fhat_mat <- stats::predict(fit, type = "terms", newdata = xmat_pred, se.fit = TRUE)

  # Critical value
  alpha <- 1 - conf_level
  crit <- qnorm(1 - alpha/2)

  # Match the "first" style: flatten predictor matrix
  plt_mat <- data.frame(
    "estimate"    = as.vector(fhat_mat$fit),
    "se"          = as.vector(fhat_mat$se.fit),
    "xind"        = as.vector(as.matrix(xmat_pred)),
    "coefficient" = rep(paste0("f", seq_len(ncol(fhat_mat$fit))),
                      each = nrow(fhat_mat$fit))
  )

  # Add confidence intervals
  plt_mat <- plt_mat %>%
    dplyr::mutate(
      LB = .data$estimate - crit * .data$se,
      UB = .data$estimate + crit * .data$se
    ) |> dplyr::arrange(.data$coefficient, .data$xind)

  return(plt_mat)
}
