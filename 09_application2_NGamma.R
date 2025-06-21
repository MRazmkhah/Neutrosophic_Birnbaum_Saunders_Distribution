# ------------------------------------------------------------------------
# Neutrosophic Parameter Estimation and Inferences for Gamma Model (Application 2)
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Load Required Libraries
# ------------------------------------------------------------------------
library(GoFKernel)

# ------------------------------------------------------------------------
# Input Data
# ------------------------------------------------------------------------
xl <- c(304.12, 355.34, 310.93, 309.47, 309.12, 292.80, 327.49,
        280.33, 259.99, 238.19, 229.98, 226.77, 223.57, 233.12,
        216.37, 208.16, 206.30, 193.44, 177.31, 157.88, 153.18,
        143.93, 132.78, 127.87, 118.28, 116.75, 116.96, 114.21,
        106.86)

xu <- c(307.82, 355.34, 310.93, 309.47, 312.10, 292.80, 327.49,
        280.33, 259.99, 242.45, 229.98, 226.77, 223.57, 233.12,
        216.37, 208.16, 209.14, 193.44, 177.31, 157.88, 153.18,
        143.93, 132.78, 127.87, 118.28, 116.75, 116.96, 114.21,
        110.62)

L <- xl - 1e-9
U <- xu + 1e-9
n <- length(xl)

# ------------------------------------------------------------------------
# Parameter Estimation Functions
# ------------------------------------------------------------------------
loglikelihood <- function(x, alpha_hat, beta_hat) {
  sum(log(dgamma(x, shape=alpha_hat, scale = beta_hat, log = FALSE)))
}

ML <- function(x) {
  n <- length(x)
  DL_a <- function(aa){
    log(aa)-log(mean(x))-digamma(aa)-sum(log(x))/n
  }
  DL_a <- Vectorize(DL_a)
  res <- tryCatch({
    a0 <- uniroot(DL_a,lower=0.01, upper=10)$root
    b0 <- mean(x)/a0
    c(a0, b0)
  }, error = function(e) c(NA, NA))
  return(res)
}

get_bounds_with_likelihood <- function(LL, UU, estimator, n_samples = 100) {
  x_samples <- replicate(n_samples, runif(length(LL), min = LL, max = UU))
  x_samples <- t(x_samples)
  x_samples <- rbind(LL, UU, x_samples)

  results <- apply(x_samples, 1, function(x) {
    est <- estimator(x)
    if (any(is.na(est))) return(c(NA, NA, NA))
    alpha_hat <- est[1]
    beta_hat <- est[2]
    ll_val <- loglikelihood(x, alpha_hat, beta_hat)
    return(c(alpha_hat, beta_hat, ll_val))
  })

  results <- t(results)
  results <- results[complete.cases(results), ]
  if (nrow(results) == 0) return(rep(NA, 6))

  alpha_range <- range(results[, 1], na.rm = TRUE)
  beta_range <- range(results[, 2], na.rm = TRUE)
  loglik_range <- range(results[, 3], na.rm = TRUE)
  return(c(alpha_range, beta_range, loglik_range))
}

G <- 100
bounds <- get_bounds_with_likelihood(L, U, ML, n_samples = G)

ML_a_L <- bounds[1]
ML_a_U <- bounds[2]
ML_b_L <- bounds[3]
ML_b_U <- bounds[4]
loglike_min <- bounds[5]
loglike_max <- bounds[6]

# ------------------------------------------------------------------------
# AIC and BIC Calculations
# ------------------------------------------------------------------------
AIC <- function(y) { 2 * 2 - 2 * y }
BIC <- function(y) { 4 * log(n) - 2 * y }

AIC_min <- AIC(loglike_max)
AIC_max <- AIC(loglike_min)
BIC_min <- BIC(loglike_max)
BIC_max <- BIC(loglike_min)

# ------------------------------------------------------------------------
# Modified KS Test
# ------------------------------------------------------------------------
F_t <- function(t,a,b){
  pgamma(t, shape=a, scale = b, lower.tail = TRUE, log.p = FALSE)
}

Modified_KS <- function(t,a_hat,b_hat){
  n <- length(t)
  t_sort <- sort(t)
  V_hat <- F_t(t_sort,a_hat,b_hat)
  Y_hat <- qnorm(V_hat)
  Z <- (Y_hat - mean(Y_hat))/sqrt(sum((Y_hat-mean(Y_hat))^2)/(n-1))
  U_hat <- pnorm(Z)
  F_n <- ecdf(t_sort)
  ECDF <- F_n(t_sort)
  D <- NULL
  for (j in 1:n){
    if (j==1){
      D[j] <- abs(U_hat[j]-ECDF[j])
    } else {
      D[j] <- max(abs(ECDF[j]-U_hat[j]),abs(U_hat[j]-ECDF[j-1]))
    }
  }
  Max_D <- max(D)
  KS_star <- (sqrt(n)-0.01+0.85/sqrt(n))*Max_D
  return(c(Max_D , KS_star))
}

DDl <- Modified_KS(xl,ML_a_L,ML_b_L)
DDu <- Modified_KS(xu,ML_a_U,ML_b_U)

# ------------------------------------------------------------------------
# Monte Carlo Simulation for P-Value
# ------------------------------------------------------------------------
nn <- 10000

KS_Montecarlo <- function(data,a,b){
  KS_data <- Modified_KS(data,a,b)[2]
  F_t_q <-function(p) qgamma(p, shape=a, scale = b)
  RV <- function(n) {
    Q <- runif(n)
    lapply(Q, FUN=F_t_q)
  }

  KS_Monte <- numeric(nn)
  for (i in 1:nn){
    sample <- unlist(RV(n))
    KS_Monte[i] <- Modified_KS(sample,a,b)[2]
  }
  mean(KS_Monte > KS_data)
}

set.seed(12345)
p_value_l <- KS_Montecarlo(xl,ML_a_L,ML_b_L)
p_value_u <- KS_Montecarlo(xu,ML_a_U,ML_b_U)

set.seed(88888888)
paste0("[",round(ML_a_L,4),",",round(ML_a_U,4),
       "] & [",round(ML_b_L,4),",",round(ML_b_U,4),
       "] & [",round(loglike_min,4),",",round(loglike_max,4),
       "] & [",round(AIC_min,4),",",round(AIC_max,4), 
       "] & [",round(BIC_min,4),",",round(BIC_max,4),
       "] & [",round(DDl[2],4),",",round(DDu[2],4),
       "] & [",round(p_value_l,4),",",round(p_value_u,4),
       "]")
