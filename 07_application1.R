# ------------------------------------------------------------------------
# Neutrosophic Parameter Estimation and Inferences for Birnbaum-Saunders Model (Application 1)
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Load Required Libraries
# ------------------------------------------------------------------------
library(GoFKernel)

# ------------------------------------------------------------------------
# Input Data
# ------------------------------------------------------------------------
# Indeterminacy noise
e <- 1  

# Real Birnbaum-Saunders data
real_data_BS <- c(70, 90, 96, 97, 99, 100, 103, 104, 104, 105, 107, 108, 108, 108, 109, 109, 112, 112, 113, 114, 114, 114,
                  116, 119, 120, 120, 120, 121, 121, 123, 124, 124, 124, 124, 124, 128, 128, 129, 129, 130, 130, 130, 131, 131,
                  131, 131, 131, 132, 132, 132, 133, 134, 134, 134, 134, 134, 136, 136, 137, 138, 138, 138, 139, 139, 141, 141,
                  142, 142, 142, 142, 142, 142, 144, 144, 145, 146, 148, 148, 149, 151, 151, 152, 155, 156, 157, 157, 157, 157,
                  158, 159, 162, 163, 163, 164, 166, 166, 168, 170, 174, 196, 212)

n <- length(real_data_BS)
L <- real_data_BS - e
U <- real_data_BS + e

# ------------------------------------------------------------------------
# Log-Likelihood Function
# ------------------------------------------------------------------------
loglikelihood <- function(x, alpha_hat, beta_hat) {
  -n * log(2) - n * log(alpha_hat) - n * log(beta_hat) - n * log(2 * pi) / 2 +
    sum(log(sqrt(beta_hat / x) + (beta_hat / x)^1.5)) -
    (1 / (2 * alpha_hat^2)) * sum(x / beta_hat + beta_hat / x - 2)
}

# ------------------------------------------------------------------------
# MLE Estimation Function
# ------------------------------------------------------------------------
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

# ------------------------------------------------------------------------
# Grid-Based Estimation Bounds
# ------------------------------------------------------------------------
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

# ------------------------------------------------------------------------
# Run Estimation
# ------------------------------------------------------------------------
G <- 100
bounds <- get_bounds_with_likelihood(L, U, ML, n_samples = G)

# Store results
ML_a_L <- bounds[1]
ML_a_U <- bounds[2]
ML_b_L <- bounds[3]
ML_b_U <- bounds[4]
loglike_min <- bounds[5]
loglike_max <- bounds[6]

# ------------------------------------------------------------------------
# AIC and BIC Calculations
# ------------------------------------------------------------------------
AIC <- function(y) {
  2 * 2 - 2 * y
}

BIC <- function(y) {
  4 * log(n) - 2 * y
}

AIC_min <- AIC(loglike_max)
AIC_max <- AIC(loglike_min)
BIC_min <- BIC(loglike_max)
BIC_max <- BIC(loglike_min)

# ------------------------------------------------------------------------
# Modified KS Test
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
    } else{
      D[j] <- max(abs(ECDF[j]-U_hat[j]),abs(U_hat[j]-ECDF[j-1]))
    }
  }

  Max_D <- max(D)
  KS_star <- (sqrt(n)-0.01+0.85/sqrt(n))*Max_D
  c(Max_D , KS_star)
}

DDl <- Modified_KS(L, ML_a_L, ML_b_L)
DDu <- Modified_KS(U, ML_a_U, ML_b_U)

# ------------------------------------------------------------------------
# Output Final Result
# ------------------------------------------------------------------------
set.seed(88888888)
result <- paste0(
  "Estimated alpha = [", round(ML_a_L, 4), ", ", round(ML_a_U, 4), "] ",
  "& Estimated beta = [", round(ML_b_L, 4), ", ", round(ML_b_U, 4), "] ",
  "& Log-likelihood = [", round(loglike_min, 4), ", ", round(loglike_max, 4), "] ",
  "& AIC = [", round(AIC_min, 4), ", ", round(AIC_max, 4), "] ",
  "& BIC = [", round(BIC_min, 4), ", ", round(BIC_max, 4), "]",
  "& Modified_KS = [", round(DDl[2], 4), ",", round(DDu[2], 4), "]" 
)
print(result)
