# ------------------------------------------------------------------------
# Clean Version: ECR Table for NBSD - Full-Interval vs. Endpoint Estimation
# ------------------------------------------------------------------------

library(GoFKernel)  # Required for inverse() function

# ------------------------------------------------------------------------
# User-Specified Parameter Configurations
# ------------------------------------------------------------------------

alpha_vals <- c(0.2)       # List of shape parameter values (α) to test
beta_vals <- c(3)          # List of scale parameter values (β) to test
epsilon_vals <- c(0.1)     # Indeterminacy levels (interval width)
n_vals <- c(100)           # Sample sizes
N <- 1000                  # Number of simulation replications for each setup

# ------------------------------------------------------------------------
# Estimator Functions for α and β (used in both methods)
# ------------------------------------------------------------------------

# Estimate beta using the method based on solving a non-linear score equation
estimate_beta <- function(x) {
  n <- length(x)
  r <- n / sum(1 / x)
  s <- mean(x)
  score <- function(beta) {
    K <- n / sum(1 / (beta + x))
    beta^2 - beta * (2 * r + K) + r * (s + K)
  }
  tryCatch(uniroot(score, lower = r, upper = s, extendInt = "yes")$root,
           error = function(e) NA)
}

# Estimate alpha by plugging in the estimated beta
estimate_alpha <- function(x) {
  n <- length(x)
  r <- n / sum(1 / x)
  s <- mean(x)
  score <- function(beta) {
    K <- n / sum(1 / (beta + x))
    beta^2 - beta * (2 * r + K) + r * (s + K)
  }
  b0 <- tryCatch(uniroot(score, lower = r, upper = s, extendInt = "yes")$root,
                 error = function(e) NA)
  if (is.na(b0)) return(NA)
  alpha_hat <- sqrt(s / b0 + b0 / r - 2)
  if (!is.finite(alpha_hat)) return(NA)
  return(alpha_hat)
}

# ------------------------------------------------------------------------
# Grid-Based Estimation Function (Used in Full-Interval Method)
# ------------------------------------------------------------------------

# This function estimates parameter bounds by sampling uniformly from intervals
get_bounds <- function(LL, UU, estimator, n_samples = 1000) {
  x_samples <- replicate(n_samples, runif(length(LL), min = LL, max = UU))
  x_samples <- t(x_samples)
  x_samples <- rbind(LL, UU, x_samples)  # Ensure edges included
  values <- apply(x_samples, 1, estimator)
  values <- values[!is.na(values)]
  return(c(min(values), max(values)))
}

# ------------------------------------------------------------------------
# Main Simulation Function for One Configuration
# ------------------------------------------------------------------------

simulate_case <- function(eps, alpha, beta, n) {
  # True CDF and its inverse for Birnbaum-Saunders with given α and β
  F_t <- function(t) {
    z <- (1 / alpha) * (sqrt(t / beta) - sqrt(beta / t))
    pnorm(z)
  }
  F_t_inv <- inverse(F_t, lower = 0.00001, upper = 1000)
  
  # Storage vectors for parameter bounds
  A_L <- A_U <- B_L <- B_U <- A_L_o <- A_U_o <- B_L_o <- B_U_o <- numeric()
  N_a <- N_b <- N_a_o <- N_b_o <- 0  # Counters for ECRs
  
  for (i in 1:N) {
    # Step 1: Generate exact values
    Q <- runif(n)
    xl <- unlist(lapply(Q, F_t_inv))
    
    # Step 2: Add indeterminacy noise to form intervals
    noise <- runif(n, min = -eps, max = eps)
    xu <- pmax(xl + noise, 1e-7)
    L <- pmin(xl, xu)
    U <- pmax(xl, xu)
    
    # Step 3a: Full-interval estimation using grid sampling
    BB <- get_bounds(L, U, estimate_beta)
    B_L[i] <- BB[1]; B_U[i] <- BB[2]
    AA <- get_bounds(L, U, estimate_alpha)
    A_L[i] <- AA[1]; A_U[i] <- AA[2]
    
    # Step 3b: Endpoint-only estimation using just L and U
    B_L_o[i] <- min(estimate_beta(L), estimate_beta(U))
    B_U_o[i] <- max(estimate_beta(L), estimate_beta(U))
    A_L_o[i] <- min(estimate_alpha(L), estimate_alpha(U))
    A_U_o[i] <- max(estimate_alpha(L), estimate_alpha(U))
    
    # Step 4: Count coverage for ECR calculations
    if (B_L[i] <= beta & beta <= B_U[i]) N_b <- N_b + 1
    if (B_L_o[i] <= beta & beta <= B_U_o[i]) N_b_o <- N_b_o + 1
    if (A_L[i] <= alpha & alpha <= A_U[i]) N_a <- N_a + 1
    if (A_L_o[i] <= alpha & alpha <= A_U_o[i]) N_a_o <- N_a_o + 1
  }
  
  # ------------------------------------------------------------------------
  # Return final output as one row of a data frame
  # ------------------------------------------------------------------------
  out <- data.frame(
    epsilon = eps, alpha = alpha, beta = beta, n = n,
    alpha_A = sprintf("[%.3f, %.3f]", mean(A_L), mean(A_U)),
    beta_A  = sprintf("[%.3f, %.3f]", mean(B_L), mean(B_U)),
    ECR_A_alpha = round(100 * N_a / N, 1),
    ECR_A_beta  = round(100 * N_b / N, 1),
    alpha_B = sprintf("[%.3f, %.3f]", mean(A_L_o), mean(A_U_o)),
    beta_B  = sprintf("[%.3f, %.3f]", mean(B_L_o), mean(B_U_o)),
    ECR_B_alpha = round(100 * N_a_o / N, 1),
    ECR_B_beta  = round(100 * N_b_o / N, 1)
  )
  return(out)
}

# ------------------------------------------------------------------------
# Run Simulation for All Configurations and Combine Results
# ------------------------------------------------------------------------

results <- do.call(rbind, lapply(epsilon_vals, function(eps) {
  do.call(rbind, lapply(alpha_vals, function(a) {
    do.call(rbind, lapply(beta_vals, function(b) {
      do.call(rbind, lapply(n_vals, function(n) {
        simulate_case(eps, a, b, n)
      }))
    }))
  }))
}))

# ------------------------------------------------------------------------
# Print Final Table
# ------------------------------------------------------------------------

print(results)
