make_data <- function(N = integer(), J = integer(), f_X = function(x) x,
                      s2_e = numeric(), p_X = function(j) runif(j,0,1),
                      seed = integer()){
  # Set seed
  set.seed(seed)

  #Simulate the data
  X <- t(vapply(1:N, function(x,j) p_X(j), j=J, numeric(J)))
  Y <- apply(X, 1, function(x, f, s2) mean(f(x)) + rnorm(1, 0, sd=sqrt(s2)),
             f=f_X, s2=s2_e)

  Lmat <- matrix(1/J, N, J)
  df <- data.frame(Y=Y, X=I(X), Lmat=I(Lmat))
  df$Xbar <- rowMeans(X)

  return(list(df = df, f_X = f_X))
}
