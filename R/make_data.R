#' Simulate data for SiZR model
#'
#' @param N Sample size (number of subjects).
#' @param J Number of covariates.
#' @param s2_e Variance of the residual error term.
#' @param seed Random seed for reproducibility.
#' @param p_X Function to generate the covariates X. Default is uniform distribution on (0, 1).
#' @param f_X Coefficient function defining the association between X and Y. Default is the identity function.
#'
#' @return A list containing: data frame with simulated data and the coefficient function f_X.
#' @export
#'
#' @examples
#' # Simulate data with default parameters
#' sim_data <- make_data(N = 100, J = 5, s2_e = 1, seed = 123)
#' head(sim_data$df)
#'
#' @importFrom stats gaussian rnorm runif
#'

make_data <- function(N = integer(), J = integer(), s2_e = numeric(), seed = integer(),
                      p_X = function(j) runif(j, 0, 1), f_X = function(x) x){
  # Set seed
  set.seed(seed)

  # Simulate the data
  # Features
  X <- t(vapply(1:N, function(x,j) p_X(j), j=J, numeric(J)))
  # Outcomes
  Y <- apply(X, 1, function(x, f, s2) mean(f(x)) + rnorm(1, 0, sd=sqrt(s2)),
             f=f_X, s2=s2_e)

  Lmat <- matrix(1/J, N, J)
  df <- data.frame(Y=Y, X=I(X), Lmat=I(Lmat))
  df$Xbar <- rowMeans(X)

  return(list(df = df, f_X = f_X))
}
