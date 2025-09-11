fitting <- function(sim_result, k = 20, method = "REML", family = gaussian(),
                    fit = "gam", model = "SiZR", data_type = c("wide", "long")){

  df <- sim_result$df
  f_X <- sim_result$f_X

  if (fit == "bam") {
    fit_func <- function(...) mgcv::bam(...)
  } else {
    fit_func <- function(...) mgcv::gam(...)
  }

  if (model == "SiZR"){
    fit <- fit_func(Y ~ s(X, by=Lmat,k=k), method=method, family=family, data=df)
  } else if (model == "mean"){
    fit <- fit_func(Y ~ s(Xbar,k=k), method=method, family=family, data=df)
  }

  return(fit = fit)
}
