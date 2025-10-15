#' Fit a model to the simulated data
#'
#' @param sim_result A list containing the simulated data frame and the coefficient function f_X, returned from make_data()
#' @param k Number of basis functions to use in the smooth term. Default is 20.
#' @param method Method for smoothing parameter estimation. Default is "REML".
#' @param family Family object specifying the error distribution and link function. Default is gaussian().
#' @param fit Type of fitting function to use. Options are "gam" (default) or "bam".
#' @param model Type of model to fit. Options are "SiZR" (default) or "mean".
#'
#'
#' @return A fitted model object.
#' @export
#'
#' @importFrom mgcv gam bam s
#'
#' @examples
#' # Simulate data
#' sim_data <- sim_data(N = 100, J = 5, s2_e = 1, seed = 123)
#' # Fit SiZR model
#' fit_SiZR <- fitting(sim_data, model = "SiZR")


fitting <- function(sim_result, k = 20, method = "REML", family = gaussian(),
                    fit = "gam", model = "SiZR"){

  df <- sim_result$df
  f_X <- sim_result$f_X

  fit_func <- if(fit == "bam") mgcv::bam else mgcv::gam

  if (model == "SiZR"){
    fit <- fit_func(Y ~ s(X, by=Lmat,k=k), method=method, family=family, data=df)
  } else if (model == "mean"){
    fit <- fit_func(Y ~ s(Xbar,k=k), method=method, family=family, data=df)
  }

  return(fit = fit)
}
