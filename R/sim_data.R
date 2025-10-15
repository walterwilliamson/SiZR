#' Simulate data for SiZR model
#'
#' @param N Sample size (number of subjects).
#' @param J Number of covariates.
#' @param s2_e Variance of the residual error term.
#' @param seed Random seed for reproducibility.
#' @param p_X Function to generate the covariates X. Default is uniform distribution on (0, 1).
#' @param f_X Coefficient function defining the association between X and Y. Default is the identity function.
#' @param family Family of the outcome variable. Options are "gaussian", "binomial", or "poisson". Default is "gaussian".
#'
#' @return A list containing: data frame with simulated data and the coefficient function f_X.
#' @export
#'
#' @examples
#' # Simulate data with default parameters
#' sim_data <- sim_data(N = 100, J = 5, s2_e = 1, seed = 123)
#' head(sim_data$df)
#'
#' @importFrom stats gaussian rnorm runif rbinom rpois
#'

sim_data <- function(N = integer(), J = integer(), s2_e = numeric(), seed = integer(),
                      p_X = function(j) runif(j, 0, 1), f_X = function(x) x,
                      family = c("gaussian", "binomial", "poisson")) {
  # Set seed
  set.seed(seed)
  family <- match.arg(family)

  # Simulate the data
  # Features
  X <- t(vapply(1:N, function(x,j) p_X(j), j=J, numeric(J)))

  # Outcomes
  eta <- apply(X, 1, function(x, f) mean(f(x)), f = f_X)

  if (family == "gaussian"){
    Y <- eta + rnorm(N, 0, sqrt(s2_e))
  }

  if (family == "poisson"){
    lambda <- exp(eta)
    Y <- rpois(N, lambda=lambda)
  }

  if (family == "binomial"){
    expit <- function(x) 1/(1+exp(-x))
    p <- expit(eta)
    Y <- rbinom(N, size=1, prob=p)
  }

  Lmat <- matrix(1/J, N, J)
  df <- data.frame(Y=Y, X=I(X), Lmat=I(Lmat))
  df$Xbar <- rowMeans(X)

  return(list(df = df, f_X = f_X))
}
