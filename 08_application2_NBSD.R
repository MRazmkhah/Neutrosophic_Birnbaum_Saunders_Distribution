
# ------------------------------------------------------------------------
# Neutrosophic Parameter Estimation and Inference for Birnbaum-Saunders Model (Application 2)
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Load Required Libraries
# ------------------------------------------------------------------------
library(GoFKernel)

# ------------------------------------------------------------------------
# Input Data
# ------------------------------------------------------------------------
# Original data points
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

# Add small offsets to avoid optimization issues
L <- xl - 1e-9
U <- xu + 1e-9
n <- length(xl)  # Sample size

# ------------------------------------------------------------------------
# Parameter Estimation Functions
# ------------------------------------------------------------------------

# Log-likelihood function
loglikelihood <- function(x, alpha_hat, beta_hat) {
  -n * log(2) - n * log(alpha_hat) - n * log(beta_hat) - n * log(2 * pi) / 2 +
    sum(log(sqrt(beta_hat / x) + (beta_hat / x)^1.5)) -
    (1 / (2 * alpha_hat^2)) * sum(x / beta_hat + beta_hat / x - 2)
}

# Estimate alpha and beta using MLE
ML <- function(x) {
  n <- length(x)
  r <- n / sum(1 / x)
  s <- sum(x) / n
  DL_b <- function(bb) {
    K <- n / sum(1 / (bb + x))
    bb^2 - bb * (2 * r + K) + r * (s + K)
  }
  DL_b <- Vectorize(DL_b)
  res <- tryCatch({
    b0 <- uniroot(DL_b, lower = r, upper = s)$root
    a0 <- sqrt(s / b0 + b0 / r - 2)
    c(a0, b0)
  }, error = function(e) c(NA, NA))
  return(res)
}

# Extended function to estimate alpha, beta, and log-likelihood bounds
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

  alpha_range   <- range(results[, 1], na.rm = TRUE)
  beta_range    <- range(results[, 2], na.rm = TRUE)
  loglik_range  <- range(results[, 3], na.rm = TRUE)

  return(c(alpha_range, beta_range, loglik_range))
}

# Run estimation on real data
G <- 100
bounds <- get_bounds_with_likelihood(L, U, ML, n_samples = G)

# Store outputs
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
# Modified KS Test Function
# ------------------------------------------------------------------------
F_t <- function(t, a, b) {
  xx <- (1 / a) * (sqrt(t / b) - sqrt(b / t))
  pnorm(xx, mean = 0, sd = 1)
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
    }
    else{
      D[j] <- max(abs(ECDF[j]-U_hat[j]),abs(U_hat[j]-ECDF[j-1]))
    }
  }

  Max_D <- max(D)
  KS_star <- (sqrt(n)-0.01+0.85/sqrt(n))*Max_D

  c(Max_D , KS_star)
}

DDl <- Modified_KS(xl, ML_a_L, ML_b_L)
DDu <- Modified_KS(xu, ML_a_U, ML_b_U)

# ------------------------------------------------------------------------
# P-Value Calculation Using Monte Carlo Simulation
# ------------------------------------------------------------------------
KS_Montecarlo <- function(data, a, b, iterations = 5000) {
  KS_data <- Modified_KS(data, a, b)[2]

  F_tt <- function(t) {
    xx <- (1 / a) * (sqrt(t / b) - sqrt(b / t))
    pnorm(xx, mean = 0, sd = 1)
  }

  F_t_inv <- inverse(F_tt, lower = 0, upper = 1e5)

  RV <- function(n) {
    sapply(runif(n), F_t_inv)
  }

  KS_Monte <- replicate(iterations, {
    sample <- RV(length(data))
    Modified_KS(sample, a, b)[2]
  })

  mean(KS_Monte > KS_data)
}

set.seed(12345)
p_value_u <- KS_Montecarlo(xl, ML_a_L, ML_b_L)
p_value_l <- KS_Montecarlo(xu, ML_a_U, ML_b_U)

# ------------------------------------------------------------------------
# Final Results
# ------------------------------------------------------------------------
results <- paste0(
  "Estimated alpha = [", round(ML_a_L, 4), ", ", round(ML_a_U, 4), "] ",
  "& Estimated beta = [", round(ML_b_L, 4), ", ", round(ML_b_U, 4), "] ",
  "& Log-likelihood = [", round(loglike_min, 4), ", ", round(loglike_max, 4), "] ",
  "& AIC = [", round(AIC_min, 4), ", ", round(AIC_max, 4), "] ",
  "& BIC = [", round(BIC_min, 4), ", ", round(BIC_max, 4), "]",
  "& Modified_KS = [", round(DDl[2], 4), ",", round(DDu[2], 4), "]" ,
  "& p-value = [", round(p_value_l, 4), ",", round(p_value_u, 4), "]"
)

print(results)
