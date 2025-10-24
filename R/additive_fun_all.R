### Library
require(fields)
require(mgcv)
library(doParallel)
library(foreach)
library(patchwork)
library(ggplot2)

### Functions of Real Coefficients
g_fun_collection <- list(
  fun1 = function(t, x) sin(11 + 0.3 * (t) / (x)),
  fun2 = function(t, x) sin(11 + t + x * 3 / 10),
  fun3 = function(t, x) sin(11 * t + t * x^2 / 10),
  fun4 = function(t, x) sin(11 * t + x^2 / 10),
  fun5 = function(t, x) sin(11 + t + x / 4),
  fun6 = function(t, x) sin(pi * t/10 + x / 4),
  fun7 = function(t, x) sin(11 * t * x / 100),
  fun8 = function(t, x) sin(pi * t/10) * tanh(x),
  fun9 = function(t, x) exp(-((x-7)/5)^2-(((t/10)^2-0.5)/0.3)^2)
)
## selected_fun <- g_fun_collection$fun2 ##select the function you want
# N is the number of participants
# J is the number of observations
# tlen is number of observations per function per subject
########################################
####          Sizer Method          ####
########################################
sizer <- function(N, J, tlen, func){
  tind <- seq(0, 10, length=tlen)
  X.arr <- array(dim=c(N,J,tlen))
  ## simulate functiona data ##
  for (i in 1:N) {
    # subject-specific coefficients
    u_i1 <- rnorm(1, mean = 0, sd = 5)       # u_i1 ~ N(0, 25)
    u_i2 <- rnorm(1, mean = 0, sd = 0.2)     # u_i2 ~ N(0, 0.04)
    # repeated noisy observations (exchangeable)
    for (j in 1:J) {
      # latent true function X_i(t)
      X_true <- u_i1 + u_i2 * tind
      for (k in 1:10) {
        v_ik1 <- rnorm(1, mean = 0, sd = 1/k)     
        v_ik2 <- rnorm(1, mean = 0, sd = 1/k)    
        X_true <- X_true + v_ik1 * sin(2 * pi * k * tind / 10) +
          v_ik2 * cos(2 * pi * k * tind / 10) ### fourier expansion ###
      }
      X.arr[i, j, ] <- X_true + rnorm(tlen, mean = 0, sd = 0.1)
    }
  }
  ## true coefficient
  g_fun <- g_fun_collection[[func]]
  ## set up matrices for FGAM
  trow <- matrix(tind, nrow = J, ncol = tlen, byrow = TRUE)
  tmat.wide <- matrix(rep(tind, N*J), ncol=tlen*J, byrow=TRUE)
  Xmat.wide <- t(apply(X.arr, 1, function(x) as.vector(t(x))))
  Lmat.wide <- matrix(1/(tlen*J), ncol=tlen*J, nrow=N)
  
  ### Gaussian ###
  Y <- vapply(1:N, function(i) mean(g_fun(trow, X.arr[i, , ])), numeric(1))
  fit_fGAM_G <- gam(Y ~ ti(tmat.wide, Xmat.wide, by = Lmat.wide,
                           bs = c('cr', 'cr'), k = c(10, 10), mc = c(F,T)), method = "REML")
  
  ### Binary ###
  eta <- 10*vapply(1:N, function(i) mean(g_fun(trow, X.arr[i, , ])), numeric(1))
  p <- 1/(1 + exp(-eta))
  Y <- rbinom(N, size = 1, prob = p)
  fit_fGAM_B <- gam(Y ~ ti(tmat.wide, Xmat.wide, by = Lmat.wide,
                           bs = c('cr', 'cr'), k = c(10, 10), mc = c(F,T)),
                    family = binomial(), method = "REML")
  
  ### Poisson ###
  eta <- 5*vapply(1:N, function(i) mean(g_fun(trow, X.arr[i, , ])), numeric(1))
  lambda <- exp(1 + eta)
  Y <- rpois(N, lambda)
  fit_fGAM_P <- gam(Y ~ ti(tmat.wide, Xmat.wide, by = Lmat.wide,
                           bs = c('cr', 'cr'), k = c(10, 10), mc = c(F,T)),
                    family = poisson(), method = "REML")
  
  # prediction grid in [0, 1] × [0, 1]
  tind <- seq(0, 10, length.out = 101)
  xind <- seq(-10, 10, length.out = 100)
  
  # estimated coefficient
  get_coefMat <- function(fit_model, xind, tind) {
    coef <- predict(
      fit_model, type = "terms", se.fit = TRUE,
      newdata = data.frame(
        Xmat.wide = rep(xind, each = length(tind)),
        tmat.wide = rep(tind,      length(xind)),
        Lmat.wide = 1
      )
    )
    matrix(coef$fit, nrow = length(xind), ncol = length(tind), byrow = TRUE)
  }
  
  models <- list(G = fit_fGAM_G, B = fit_fGAM_B, P = fit_fGAM_P)
  coefMats <- lapply(models, get_coefMat, xind = xind, tind = tind)
  return(coefMats)
}
test = sizer(1000, 10, 100, 'fun6')

# prediction grid
tind <- seq(0, 10, length.out = 101)
xind <- seq(-10, 10, length.out = 100)

par(mfrow = c(2,2))

coef_list <- list(test$G, test$B, test$P)
titles <- c("Gaussian", "Binary", "Poisson")

# plot estimated coefficients
for (i in 1:3) {
  image.plot(
    tind, xind, t(coef_list[[i]]),
    main = paste("Estimated Coefficient Function -", titles[i]),
    xlab = "t", ylab = "x"
  )
}

# true surface
t <- seq(0, 10, len = 101)
x <- seq(-10, 10, len = 100)
z <- outer(t, x, function(t, x) exp(-((x-7)/5)^2-(((t/10)^2-0.5)/0.3)^2))
image.plot(t, x, z, main = "True Coefficient Function", xlab = "t", ylab = "x")


#########################################
####           Mean Method          ####
#########################################
mean_original <- function(N, J, tlen, func){
  tind <- seq(0, 10, length=tlen)
  X.arr <- array(dim=c(N,J,tlen))
  for (i in 1:N) {
    # subject-specific coefficients
    u_i1 <- rnorm(1, mean = 0, sd = 5)       # u_i1 ~ N(0, 25)
    u_i2 <- rnorm(1, mean = 0, sd = 0.2)     # u_i2 ~ N(0, 0.04)
    # repeated noisy observations (exchangeable)
    for (j in 1:J) {
      # latent true function X_i(t)
      X_true <- u_i1 + u_i2 * tind
      for (k in 1:10) {
        v_ik1 <- rnorm(1, mean = 0, sd = 1/k)     
        v_ik2 <- rnorm(1, mean = 0, sd = 1/k)    
        X_true <- X_true + v_ik1 * sin(2 * pi * k * tind / 10) +
          v_ik2 * cos(2 * pi * k * tind / 10)
      }
      X.arr[i, j, ] <- X_true + rnorm(tlen, mean = 0, sd = 0.1)
    }
  }
  g_fun <- g_fun_collection[[func]]
  ## set up matrices for FGAM
  Xbar <- apply(X.arr, c(1,3), mean)
  trow_N <- matrix(tind, nrow = N, ncol = tlen, byrow = TRUE)
  tmat <- matrix(tind, nrow = N, ncol = tlen, byrow = TRUE)
  Lmat <- matrix(1/tlen, nrow = N, ncol = tlen)
  
  ### Gaussian ###
  Y <- rowMeans(g_fun(trow_N, Xbar))
  fit_fGAM_G <- gam(Y ~ ti(tmat, Xbar, by = Lmat,
                           bs = c('cr', 'cr'), k = c(10, 10), mc = c(F,T)), 
                    method = "REML")
  
  ### Binary ###
  eta <- 10*rowMeans(g_fun(trow_N, Xbar))
  p <- 1/(1 + exp(-eta))
  Y <- rbinom(N, size = 1, prob = p)
  fit_fGAM_B <- gam(Y ~ ti(tmat, Xbar, by = Lmat,
                           bs = c('cr', 'cr'), k = c(10, 10), mc = c(F,T)),
                    family = binomial(), method = "REML")
  
  ### Poisson ###
  eta <- 5*rowMeans(g_fun(trow_N, Xbar))
  lambda <- exp(1 + eta)
  Y <- rpois(N, lambda)
  fit_fGAM_P <- gam(Y ~ ti(tmat, Xbar, by = Lmat,
                           bs = c('cr', 'cr'), k = c(10, 10), mc = c(F,T)),
                    family = poisson(), method = "REML")
  
  # prediction grid in [0, 1] × [0, 1]
  tind <- seq(0, 10, length.out = 101)
  xind <- seq(-10, 10, length.out = 100)
  
  # estimated coefficient
  get_coefMat <- function(fit_model, xind, tind) {
    coef <- predict(
      fit_model, type = "terms", se.fit = TRUE,
      newdata = data.frame(
        Xbar = rep(xind, each = length(tind)),
        tmat = rep(tind, length(xind)),
        Lmat = rep(1, length(xind) * length(tind))
      )
    )
    matrix(coef$fit, nrow = length(xind), ncol = length(tind), byrow = TRUE)
  }
  
  models <- list(G = fit_fGAM_G, B = fit_fGAM_B, P = fit_fGAM_P)
  coefMats <- lapply(models, get_coefMat, xind = xind, tind = tind)
  return(coefMats)
}
test = mean_original(5000, 10, 100, 'fun7')

# prediction grid
tind <- seq(0, 10, length.out = 101)
xind <- seq(-10, 10, length.out = 100)

par(mfrow = c(2,2))

coef_list <- list(test$G, test$B, test$P)
titles <- c("Gaussian", "Binary", "Poisson")

## plot of estimated coefficients
for (i in 1:3) {
  image.plot(
    tind, xind, t(coef_list[[i]]),
    main = paste("Estimated Coefficient Function -", titles[i]),
    xlab = "t", ylab = "x"
  )
}

# true surface
t <- seq(0, 1, len = 100)
x <- seq(-10, 10, len = 100)
z <- outer(t, x, function(t, x) sin(11 * t * x / 10))
image.plot(t, x, z, main = "True Coefficient Function", xlab = "t", ylab = "x")


#########################################
####         Multiple Results        ####
#########################################
# parallel computation
n.cores <- parallel::detectCores() - 11
cl <- makeCluster(n.cores)
registerDoParallel(cl)
# initial
MSE_GM <- MSE_PM <- MSE_BM <- Bias_GM <- Bias_PM <- Bias_BM <- c()
MSE_GS <- MSE_PS <- MSE_BS <- Bias_GS <- Bias_PS <- Bias_BS <- c()
G_TotalM <- P_TotalM <- B_TotalM <- matrix(0, nrow = 100, ncol = 101)
G_TotalS <- P_TotalS <- B_TotalS <- matrix(0, nrow = 100, ncol = 101)
results <- foreach(i = 1:20, .combine = 'rbind', .packages = c("mgcv")) %dopar% {
  # estimated coefficients
  test_M <- mean_original(5000, 10, 100, 'fun9')
  test_S <- sizer(5000, 10, 100, 'fun9')
  
  # true functions
  t <- seq(0, 10, len = 101)
  x <- seq(-10, 10, len = 100)
  z <- outer(t, x, function(t, x) exp(-((x-7)/5)^2-(((t/10)^2-0.5)/0.3)^2))
  
  # MSE & Bias
  MSE_GM <- mean((t(test_M$G) - z)^2)
  MSE_PM <- mean((t(test_M$P)/5 - z)^2)
  MSE_BM <- mean((t(test_M$B)/10 - z)^2)
  Bias_GM <- mean(abs(t(test_M$G) - z))
  Bias_PM <- mean(abs(t(test_M$P)/5 - z))
  Bias_BM <- mean(abs(t(test_M$B)/10 - z))
  
  MSE_GS <- mean((t(test_S$G) - z)^2)
  MSE_PS <- mean((t(test_S$P)/5 - z)^2)
  MSE_BS <- mean((t(test_S$B)/10 - z)^2)
  Bias_GS <- mean(abs(t(test_S$G) - z))
  Bias_PS <- mean(abs(t(test_S$P)/5 - z))
  Bias_BS <- mean(abs(t(test_S$B)/10 - z))
  
  list(
    MSE_GM = MSE_GM, MSE_PM = MSE_PM, MSE_BM = MSE_BM,
    Bias_GM = Bias_GM, Bias_PM = Bias_PM, Bias_BM = Bias_BM,
    MSE_GS = MSE_GS, MSE_PS = MSE_PS, MSE_BS = MSE_BS,
    Bias_GS = Bias_GS, Bias_PS = Bias_PS, Bias_BS = Bias_BS,
    G_M = test_M$G, P_M = test_M$P/5, B_M = test_M$B/10,
    G_S = test_S$G, P_S = test_S$P/5, B_S = test_S$B/10
  )
}
for (i in 1:20) {
  G_TotalM <- G_TotalM + results[i,13][[1]]
  P_TotalM <- P_TotalM + results[i,14][[1]]
  B_TotalM <- B_TotalM + results[i,15][[1]]
  G_TotalS <- G_TotalS + results[i,16][[1]]
  P_TotalS <- P_TotalS + results[i,17][[1]]
  B_TotalS <- B_TotalS + results[i,18][[1]]
}
stopCluster(cl)
#2D heatmap
par(mfrow = c(2,3))
coef_list <- list(G_TotalM/20, B_TotalM/20, P_TotalM/20,
                   G_TotalS/20, B_TotalS/20, P_TotalS/20)
titles <- c("Gaussian_Mean", "Binary_Mean", "Poisson_Mean",
            "Gaussian_Sizer", "Binary_Sizer", "Poisson_Sizer")
tind <- seq(0, 10, length.out = 101)
xind <- seq(-10, 10, length.out = 100)
for (i in 1:6) {
  image.plot(
    tind, xind, t(coef_list[[i]]),
    main = paste("Estimated Coefficient Function -", titles[i]),
    xlab = "t", ylab = "x"
  )
}

MSE_GM_vec <- unlist(results[,1])
MSE_PM_vec <- unlist(results[,2])
MSE_BM_vec <- unlist(results[,3])
Bias_GM_vec <- unlist(results[,4])
Bias_PM_vec <- unlist(results[,5])
Bias_BM_vec <- unlist(results[,6])
MSE_GS_vec <- unlist(results[,7])
MSE_PS_vec <- unlist(results[,8])
MSE_BS_vec <- unlist(results[,9])
Bias_GS_vec <- unlist(results[,10])
Bias_PS_vec <- unlist(results[,11])
Bias_BS_vec <- unlist(results[,12])

par(mfrow = c(2,2))

# Model M - MSE
boxplot(MSE_GM_vec, MSE_PM_vec, MSE_BM_vec, 
        names = c("Gaussian", "Poisson", "Binary"),
        main = "Model M - MSE",
        xlab = NULL, ylab = "MSE",
        col = c("#1f77b4", "#2ca02c", "#ff7f0e"),
        border = "black",
        outline = TRUE)

# Model M - Bias
boxplot(Bias_GM_vec, Bias_PM_vec, Bias_BM_vec, 
        names = c("Gaussian", "Poisson", "Binary"),
        main = "Model M - Bias",
        xlab = NULL, ylab = "Bias",
        col = c("#1f77b4", "#2ca02c", "#ff7f0e"),
        border = "black",
        outline = TRUE)

# Model S - MSE
boxplot(MSE_GS_vec, MSE_PS_vec, MSE_BS_vec, 
        names = c("Gaussian", "Poisson", "Binary"),
        main = "Model S - MSE",
        xlab = NULL, ylab = "MSE",
        col = c("#1f77b4", "#2ca02c", "#ff7f0e"),
        border = "black",
        outline = TRUE)

# Model S - Bias
boxplot(Bias_GS_vec, Bias_PS_vec, Bias_BS_vec, 
        names = c("Gaussian", "Poisson", "Binary"),
        main = "Model S - Bias",
        xlab = NULL, ylab = "Bias",
        col = c("#1f77b4", "#2ca02c", "#ff7f0e"),
        border = "black",
        outline = TRUE)


#########################################
####           Mixed Method          ####
#########################################
#### Mixed Model ######
#### Simulations Results:
# 1.generate from Sizer and estimation from Mean
# 2.generate from Mean and estimation from Sizer
mixed <- function(N, J, tlen, func){
  tind <- seq(0, 10, length=tlen)
  X.arr <- array(dim=c(N,J,tlen))
  for (i in 1:N) {
    # subject-specific coefficients
    u_i1 <- rnorm(1, mean = 0, sd = 5)       # u_i1 ~ N(0, 25)
    u_i2 <- rnorm(1, mean = 0, sd = 0.2)     # u_i2 ~ N(0, 0.04)
    # repeated noisy observations (exchangeable)
    for (j in 1:J) {
      # latent true function X_i(t)
      X_true <- u_i1 + u_i2 * tind
      for (k in 1:10) {
        v_ik1 <- rnorm(1, mean = 0, sd = 1/k)     
        v_ik2 <- rnorm(1, mean = 0, sd = 1/k)    
        X_true <- X_true + v_ik1 * sin(2 * pi * k * tind / 10) +
          v_ik2 * cos(2 * pi * k * tind / 10)
      }
      X.arr[i, j, ] <- X_true + rnorm(tlen, mean = 0, sd = 0.1)
    }
  }
  g_fun <- g_fun_collection[[func]]
  
  #### Sizer
  trow <- matrix(tind, nrow = J, ncol = tlen, byrow = TRUE)
  tmat.wide <- matrix(rep(tind, N*J), ncol=tlen*J, byrow=TRUE)
  Xmat.wide <- t(apply(X.arr, 1, function(x) as.vector(t(x))))
  Lmat.wide <- matrix(1/(tlen*J), ncol=tlen*J, nrow=N)
  
  #### Mean_Error
  Xbar <- apply(X.arr, c(1,3), mean)
  trow_N <- matrix(tind, nrow = N, ncol = tlen, byrow = TRUE)
  tmat <- matrix(tind, nrow = N, ncol = tlen, byrow = TRUE)
  Lmat <- matrix(1/tlen, nrow = N, ncol = tlen)
  
  ### Gaussian ###
  Y_Sizer <- vapply(1:N, function(i) mean(g_fun(trow, X.arr[i, , ])), numeric(1))
  fit_fGAM_G_Sizer <- gam(Y_Sizer ~ ti(tmat, Xbar, by = Lmat,
                                bs = c('cr', 'cr'), k = c(10, 10), mc = c(F, T)), 
                          method = "REML")
  
  Y_Mean <- rowMeans(g_fun(trow_N, Xbar))
  fit_fGAM_G_Mean <- gam(Y_Mean ~ ti(tmat.wide, Xmat.wide, by = Lmat.wide,
                                 bs = c('cr', 'cr'), k = c(10, 10), mc = c(F, T)), 
                         method = "REML")
  
  ### Binary ###
  eta <- 10*vapply(1:N, function(i) mean(g_fun(trow, X.arr[i, , ])), numeric(1))
  p <- 1/(1 + exp(-eta))
  Y_Sizer <- rbinom(N, size = 1, prob = p)
  fit_fGAM_B_Sizer <- gam(Y_Sizer ~ ti(tmat, Xbar, by = Lmat,
                           bs = c('cr', 'cr'), k = c(10, 10), mc = c(F, T)),
                    family = binomial(), method = "REML")
  
  eta <- 10*rowMeans(g_fun(trow_N, Xbar))
  p <- 1/(1 + exp(-eta))
  Y_Mean <- rbinom(N, size = 1, prob = p)
  fit_fGAM_B_Mean <- gam(Y_Mean ~ ti(tmat.wide, Xmat.wide, by = Lmat.wide,
                           bs = c('cr', 'cr'), k = c(10, 10), mc = c(F, T)),
                    family = binomial(), method = "REML")
  
  ### Poisson ###
  eta <- 5*vapply(1:N, function(i) mean(g_fun(trow, X.arr[i, , ])), numeric(1))
  lambda <- exp(1 + eta)
  Y_Sizer <- rpois(N, lambda)
  fit_fGAM_P_Sizer <- gam(Y_Sizer ~ ti(tmat, Xbar, by = Lmat,
                           bs = c('cr', 'cr'), k = c(10, 10), mc = c(F, T)),
                    family = poisson(), method = "REML")
  
  eta <- 5*rowMeans(g_fun(trow_N, Xbar))
  lambda <- exp(1 + eta)
  Y_Mean <- rpois(N, lambda)
  fit_fGAM_P_Mean <- gam(Y_Mean ~ ti(tmat.wide, Xmat.wide, by = Lmat.wide,
                           bs = c('cr', 'cr'), k = c(10, 10), mc = c(F, T)),
                    family = poisson(), method = "REML")
  
  # prediction grid in [0, 1] × [0, 1]
  tind <- seq(0, 10, length.out = 101)
  xind <- seq(-10, 10, length.out = 100)
  
  # 1. coefficients for Sizer
  get_coefMat_Sizer <- function(fit_model, xind, tind) {
    coef <- predict(
      fit_model, type = "terms", se.fit = TRUE,
      newdata = data.frame(
        Xbar = rep(xind, each = length(tind)),
        tmat = rep(tind, length(xind)),
        Lmat = 1
      )
    )
    matrix(coef$fit, nrow = length(xind), ncol = length(tind), byrow = TRUE)
  }
  
  # 2. coefficient for Mean
  get_coefMat_Mean <- function(fit_model, xind, tind) {
    coef <- predict(
      fit_model, type = "terms", se.fit = TRUE,
      newdata = data.frame(
        Xmat.wide = rep(xind, each = length(tind)),
        tmat.wide = rep(tind,      length(xind)),
        Lmat.wide = 1
      )
    )
    matrix(coef$fit, nrow = length(xind), ncol = length(tind), byrow = TRUE)
  }
  
  coefMats <- list(
    G_Sizer = get_coefMat_Sizer(fit_fGAM_G_Sizer, xind, tind),
    G_Mean  = get_coefMat_Mean(fit_fGAM_G_Mean,  xind, tind),
    B_Sizer = get_coefMat_Sizer(fit_fGAM_B_Sizer, xind, tind),
    B_Mean  = get_coefMat_Mean(fit_fGAM_B_Mean,  xind, tind),
    P_Sizer = get_coefMat_Sizer(fit_fGAM_P_Sizer, xind, tind),
    P_Mean  = get_coefMat_Mean(fit_fGAM_P_Mean,  xind, tind)
  )
  return(coefMats)
}

test = mixed(5000, 10, 101, 'fun6')

# prediction grid
tind <- seq(0, 10, length.out = 101)
xind <- seq(-10, 10, length.out = 100)

par(mfrow = c(2,3))

coef_list <- list(test$G_Sizer, test$B_Sizer, test$P_Sizer, 
                  test$G_Mean, test$B_Mean, test$P_Mean)
titles <- c("Gaussian-Sizer_Mean", "Binary-Sizer_Mean","Poisson-Sizer_Mean",
            "Gaussian-Mean_Sizer", "Binary-Mean_Sizer","Poisson-Mean_Sizer")

for (i in 1:6) {
  image.plot(
    tind, xind, t(coef_list[[i]]),
    main = paste("Coefficient-", titles[i]),
    xlab = "t", ylab = "x"
  )
}


#########################################
####         Multiple Results        ####
#########################################
# parallel computation
n.cores <- parallel::detectCores() - 10
cl <- makeCluster(n.cores)
registerDoParallel(cl)
# initial
MSE_GSM <- MSE_PSM <- MSE_BSM <- Bias_GSM <- Bias_PSM <- Bias_BSM <- c()
MSE_GMS <- MSE_PMS <- MSE_BMS <- Bias_GMS <- Bias_PMS <- Bias_BMS <- c()
G_TotalSM <- P_TotalSM <- B_TotalSM <- matrix(0, nrow = 100, ncol = 101)
G_TotalMS <- P_TotalMS <- B_TotalMS <- matrix(0, nrow = 100, ncol = 101)
results <- foreach(i = 1:20, .combine = 'rbind', .packages = c("mgcv")) %dopar% {
  # estimated coefficients
  test <- mixed(5000, 10, 100, 'fun8')
  
  # true functions
  t <- seq(0, 10, len = 101)
  x <- seq(-10, 10, len = 100)
  z <- outer(t, x, function(t, x) sin(pi * t/10) * tanh(x))
  
  # MSE & Bias
  MSE_GSM <- mean((t(test$G_Sizer) - z)^2)
  MSE_PSM <- mean((t(test$P_Sizer)/5 - z)^2)
  MSE_BSM <- mean((t(test$B_Sizer)/10 - z)^2)
  Bias_GSM <- mean(abs(t(test$G_Sizer) - z))
  Bias_PSM <- mean(abs(t(test$P_Sizer)/5 - z))
  Bias_BSM <- mean(abs(t(test$B_Sizer)/10 - z))
  
  MSE_GMS <- mean((t(test$G_Mean) - z)^2)
  MSE_PMS <- mean((t(test$P_Mean)/5 - z)^2)
  MSE_BMS <- mean((t(test$B_Mean)/10 - z)^2)
  Bias_GMS <- mean(abs(t(test$G_Mean) - z))
  Bias_PMS <- mean(abs(t(test$P_Mean)/5 - z))
  Bias_BMS <- mean(abs(t(test$B_Mean)/10 - z))
  
  list(
    MSE_GSM = MSE_GSM, MSE_PSM = MSE_PSM, MSE_BSM = MSE_BSM,
    Bias_GSM = Bias_GSM, Bias_PSM = Bias_PSM, Bias_BSM = Bias_BSM,
    MSE_GMS = MSE_GMS, MSE_PMS = MSE_PMS, MSE_BMS = MSE_BMS,
    Bias_GMS = Bias_GMS, Bias_PMS = Bias_PMS, Bias_BMS = Bias_BMS,
    G_SM = test$G_Sizer, P_SM = test$P_Sizer/5, B_SM = test$B_Sizer/10,
    G_MS = test$G_Mean, P_MS = test$P_Mean/5, B_MS = test$B_Mean/10
  )
}
for (i in 1:20) {
  G_TotalSM <- G_TotalSM + results[i,13][[1]]
  P_TotalSM <- P_TotalSM + results[i,14][[1]]
  B_TotalSM <- B_TotalSM + results[i,15][[1]]
  G_TotalMS <- G_TotalMS + results[i,16][[1]]
  P_TotalMS <- P_TotalMS + results[i,17][[1]]
  B_TotalMS <- B_TotalMS + results[i,18][[1]]
}
stopCluster(cl)
#2D heatmap
par(mfrow = c(2,3))
coef_list <- list(G_TotalSM/20, B_TotalSM/20, P_TotalSM/20,
                  G_TotalMS/20, B_TotalMS/20, P_TotalMS/20)
titles <- c("Gaussian-Sizer_Mean", "Binary-Sizer_Mean","Poisson-Sizer_Mean",
            "Gaussian-Mean_Sizer", "Binary-Mean_Sizer","Poisson-Mean_Sizer")
tind <- seq(0, 10, length.out = 101)
xind <- seq(-10, 10, length.out = 100)
for (i in 1:6) {
  image.plot(
    tind, xind, t(coef_list[[i]]),
    main = paste("Estimated Coefficient Function -", titles[i]),
    xlab = "t", ylab = "x"
  )
}

MSE_GM_vec <- unlist(results[,1])
MSE_PM_vec <- unlist(results[,2])
MSE_BM_vec <- unlist(results[,3])
Bias_GM_vec <- unlist(results[,4])
Bias_PM_vec <- unlist(results[,5])
Bias_BM_vec <- unlist(results[,6])
MSE_GS_vec <- unlist(results[,7])
MSE_PS_vec <- unlist(results[,8])
MSE_BS_vec <- unlist(results[,9])
Bias_GS_vec <- unlist(results[,10])
Bias_PS_vec <- unlist(results[,11])
Bias_BS_vec <- unlist(results[,12])

par(mfrow = c(2,2))

# Model M - MSE
boxplot(MSE_GM_vec, MSE_PM_vec, MSE_BM_vec, 
        names = c("Gaussian", "Poisson", "Binary"),
        main = "Model SM - MSE",
        xlab = NULL, ylab = "MSE",
        col = c("#1f77b4", "#2ca02c", "#ff7f0e"),
        border = "black",
        outline = TRUE)

# Model M - Bias
boxplot(Bias_GM_vec, Bias_PM_vec, Bias_BM_vec, 
        names = c("Gaussian", "Poisson", "Binary"),
        main = "Model SM - Bias",
        xlab = NULL, ylab = "Bias",
        col = c("#1f77b4", "#2ca02c", "#ff7f0e"),
        border = "black",
        outline = TRUE)

# Model S - MSE
boxplot(MSE_GS_vec, MSE_PS_vec, MSE_BS_vec, 
        names = c("Gaussian", "Poisson", "Binary"),
        main = "Model MS - MSE",
        xlab = NULL, ylab = "MSE",
        col = c("#1f77b4", "#2ca02c", "#ff7f0e"),
        border = "black",
        outline = TRUE)

# Model S - Bias
boxplot(Bias_GS_vec, Bias_PS_vec, Bias_BS_vec, 
        names = c("Gaussian", "Poisson", "Binary"),
        main = "Model MS - Bias",
        xlab = NULL, ylab = "Bias",
        col = c("#1f77b4", "#2ca02c", "#ff7f0e"),
        border = "black",
        outline = TRUE)


#########################################
####        Predicted Results        ####
#########################################
###############################################
##  Functional Simulation + Prediction Test  ##
###############################################

library(mgcv)

## ----- 1. Function to simulate functional data -----
X_arr <- function(N, J, tlen) {
  tind <- seq(0, 10, length.out = tlen)
  X.arr <- array(dim = c(N, J, tlen))
  
  for (i in 1:N) {
    u_i1 <- rnorm(1, mean = 0, sd = 5)
    u_i2 <- rnorm(1, mean = 0, sd = 0.2)
    
    for (j in 1:J) {
      X_true <- u_i1 + u_i2 * tind
      for (k in 1:10) {
        v_ik1 <- rnorm(1, mean = 0, sd = 1 / k)
        v_ik2 <- rnorm(1, mean = 0, sd = 1 / k)
        X_true <- X_true +
          v_ik1 * sin(2 * pi * k * tind / 10) +
          v_ik2 * cos(2 * pi * k * tind / 10)
      }
      X.arr[i, j, ] <- X_true + rnorm(tlen, mean = 0, sd = 0.1)
    }
  }
  return(X.arr)
}

## ----- 2. Main simulation + prediction function -----
predict_fun1 <- function(N, N_new, J, tlen, func) {
  
  tind <- seq(0, 10, length.out = tlen)
  X.arr <- X_arr(N, J, tlen)
  g_fun <- g_fun_collection[[func]]
  
  # --- Setup ---
  trow <- matrix(tind, nrow = J, ncol = tlen, byrow = TRUE)
  tmat.wide <- matrix(rep(tind, N * J), ncol = tlen * J, byrow = TRUE)
  Xmat.wide <- t(apply(X.arr, 1, function(x) as.vector(t(x))))
  Lmat.wide <- matrix(1 / (tlen * J), ncol = tlen * J, nrow = N)
  
  Xbar <- apply(X.arr, c(1, 3), mean)
  trow_N <- matrix(tind, nrow = N, ncol = tlen, byrow = TRUE)
  tmat <- matrix(tind, nrow = N, ncol = tlen, byrow = TRUE)
  Lmat <- matrix(1 / tlen, nrow = N, ncol = tlen)
  
  ### ---- 2.1 Fit Sizer models ----
  Y <- vapply(1:N, function(i) mean(g_fun(trow, X.arr[i,,])), numeric(1))
  fit_fGAM_G <- gam(Y ~ ti(tmat.wide, Xmat.wide, by = Lmat.wide,
                           bs = c("cr", "cr"), k = c(10, 10), mc = c(F, T)), method = "REML")
  
  eta <- 10 * vapply(1:N, function(i) mean(g_fun(trow, X.arr[i,,])), numeric(1))
  p <- 1 / (1 + exp(-eta))
  Y <- rbinom(N, size = 1, prob = p)
  fit_fGAM_B <- gam(Y ~ ti(tmat.wide, Xmat.wide, by = Lmat.wide,
                           bs = c("cr", "cr"), k = c(10, 10), mc = c(F, T)),
                    family = binomial(), method = "REML")
  
  eta <- 5 * vapply(1:N, function(i) mean(g_fun(trow, X.arr[i,,])), numeric(1))
  lambda <- exp(1 + eta)
  Y <- rpois(N, lambda)
  fit_fGAM_P <- gam(Y ~ ti(tmat.wide, Xmat.wide, by = Lmat.wide,
                           bs = c("cr", "cr"), k = c(10, 10), mc = c(F, T)),
                    family = poisson(), method = "REML")
  
  ### ---- 2.2 Fit Mean models ----
  Y <- rowMeans(g_fun(trow_N, Xbar))
  fit_fGAM_G_M <- gam(Y ~ ti(tmat, Xbar, by = Lmat,
                             bs = c("cr", "cr"), k = c(10, 10), mc = c(F, T)),
                      method = "REML")
  
  eta <- 10 * rowMeans(g_fun(trow_N, Xbar))
  p <- 1 / (1 + exp(-eta))
  Y <- rbinom(N, size = 1, prob = p)
  fit_fGAM_B_M <- gam(Y ~ ti(tmat, Xbar, by = Lmat,
                             bs = c("cr", "cr"), k = c(10, 10), mc = c(F, T)),
                      family = binomial(), method = "REML")
  
  eta <- 5 * rowMeans(g_fun(trow_N, Xbar))
  lambda <- exp(1 + eta)
  Y <- rpois(N, lambda)
  fit_fGAM_P_M <- gam(Y ~ ti(tmat, Xbar, by = Lmat,
                             bs = c("cr", "cr"), k = c(10, 10), mc = c(F, T)),
                      family = poisson(), method = "REML")
  
  ## ----- 3. Generate new data -----
  X.arr.new <- X_arr(N_new, J, tlen)
  tmat.wide.new <- matrix(rep(tind, N_new * J), ncol = tlen * J, byrow = TRUE)
  Xmat.wide.new <- t(apply(X.arr.new, 1, function(x) as.vector(t(x))))
  Lmat.wide.new <- matrix(1 / (tlen * J), ncol = tlen * J, nrow = N_new)
  
  Xbar.new <- apply(X.arr.new, c(1, 3), mean)
  tmat.new <- matrix(tind, nrow = N_new, ncol = tlen, byrow = TRUE)
  Lmat.new <- matrix(1 / tlen, nrow = N_new, ncol = tlen)
  
  ## ----- 4. Predictions -----
  new_df_sizer <- data.frame(
    tmat.wide = I(tmat.wide.new),
    Xmat.wide = I(Xmat.wide.new),
    Lmat.wide = I(Lmat.wide.new)
  )
  
  # Gaussian
  Y_hat_G_subject <- predict(fit_fGAM_G, newdata = new_df_sizer, type = "response")
  #Y_hat_G_subject <- rowMeans(matrix(Y_hat_G_pts, nrow = N_new, ncol = J * tlen, byrow = TRUE))
  
  # Binary
  eta_hat_B_pts <- predict(fit_fGAM_B, newdata = new_df_sizer, type = "link")
  #eta_hat_B_subj <- rowMeans(matrix(eta_hat_B_pts, nrow = N_new, ncol = J * tlen, byrow = TRUE))
  p_hat_B <- plogis(eta_hat_B_pts)
  
  # Poisson
  eta_hat_P_pts <- predict(fit_fGAM_P, newdata = new_df_sizer, type = "link")
  #eta_hat_P_subj <- rowMeans(matrix(eta_hat_P_pts, nrow = N_new, ncol = J * tlen, byrow = TRUE))
  lambda_hat <- exp(eta_hat_P_pts)
  
  ## Mean models
  new_df_mean <- data.frame(
    tmat = I(tmat.new),
    Xbar = I(Xbar.new),
    Lmat = I(Lmat.new)
  )
  
  Y_hat_G_M <- predict(fit_fGAM_G_M, newdata = new_df_mean, type = "response")
  
  eta_hat_B_M <- predict(fit_fGAM_B_M, newdata = new_df_mean, type = "link")
  p_hat_B_M   <- plogis(eta_hat_B_M)
  
  eta_hat_P_M <- predict(fit_fGAM_P_M, newdata = new_df_mean, type = "link")
  lambda_hat_M <- exp(eta_hat_P_M)
  
  ## ----- 5. True values -----
  Y_G_true <- vapply(1:N_new, function(i) mean(g_fun(trow, X.arr.new[i,,])), 0.0)
  eta_B_true <- 10 * vapply(1:N_new, function(i) mean(g_fun(trow, X.arr.new[i,,])), 0.0)
  p_B_true <- plogis(eta_B_true)
  Y_B_true <- rbinom(N_new, 1, p_B_true)
  
  eta_P_true <- 5 * vapply(1:N_new, function(i) mean(g_fun(trow, X.arr.new[i,,])), 0.0)
  lambda_true <- exp(1 + eta_P_true)
  Y_P_true <- rpois(N_new, lambda_true)
  
  ## Mean model true
  Y_G_M_true <- rowMeans(g_fun(tmat.new, Xbar.new))
  eta_B_M_true <- 10 * rowMeans(g_fun(tmat.new, Xbar.new))
  p_B_M_true <- plogis(eta_B_M_true)
  Y_B_M_true <- rbinom(N_new, 1, p_B_M_true)
  
  eta_P_M_true <- 5 * rowMeans(g_fun(tmat.new, Xbar.new))
  lambda_true_M <- exp(1 + eta_P_M_true)
  Y_P_M_true <- rpois(N_new, lambda_true_M)
  
  ## ----- 6. Metrics -----
  MSE <- list(
    G_MSE = mean((Y_G_true - Y_hat_G_subject)^2),
    B_MSE_prob = mean((p_B_true - p_hat_B)^2),
    P_MSE = mean((log(lambda_true) - log(lambda_hat))^2),
    
    G_MSE_M = mean((Y_G_M_true - Y_hat_G_M)^2),
    B_MSE_M_prob = mean((p_B_M_true - p_hat_B_M)^2),
    P_MSE_M = mean((log(lambda_true_M) - log(lambda_hat_M))^2)
  )
  return(MSE)
}

test <- predict_fun1(5000, 1000, 10, 100, "fun9")
print(test)


