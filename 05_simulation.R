
# ------------------------------------------------------------------------
# Simulation Study for Neutrosophic Birnbaum-Saunders Model (Grid-based)
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Load Required Libraries
# ------------------------------------------------------------------------
library(GoFKernel)
library(parallel)

# ------------------------------------------------------------------------
# Simulation Parameters
# ------------------------------------------------------------------------
N <- 1000  # Number of simulation repetitions
G <- 100   # Grid points per repetition
a <- 0.5
b <- 1.5
e <- 0.1

# ------------------------------------------------------------------------
# Simulation Function
# ------------------------------------------------------------------------
f2 <- function(n) {
  A_L <- A_U <- B_L <- B_U <- numeric()

  F_t <- function(t) {
    xx <- (1 / a) * (sqrt(t / b) - sqrt(b / t))
    pnorm(xx)
  }

  F_t.inv <- inverse(F_t, lower = 0.001, upper = 500)

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

  # Get bounds using grid
  get_bounds <- function(LL, UU, estimator, n_samples = G) {
    x_samples <- replicate(n_samples, runif(length(LL), min = LL, max = UU))
    x_samples <- t(x_samples)
    x_samples <- rbind(LL, UU, x_samples)

    values <- t(apply(x_samples, 1, estimator))
    values <- values[complete.cases(values), ]

    if (nrow(values) == 0) return(rep(NA, 4))

    alpha_range <- range(values[, 1], na.rm = TRUE)
    beta_range  <- range(values[, 2], na.rm = TRUE)
    return(c(alpha_range, beta_range))
  }

  for (i in 1:N) {
    RV <- function(n) {
      Q <- runif(n)
      sapply(Q, F_t.inv)
    }

    xl <- unlist(RV(n))
    NN <- runif(n, min = -e, max = e)
    xu <- xl + NN

    MM <- cbind(xl, xu)
    L <- apply(MM, 1, min)
    U <- apply(MM, 1, max)

    bounds <- get_bounds(L, U, ML, n_samples = G)
    if (any(is.na(bounds))) next

    A_L[i] <- bounds[1]
    A_U[i] <- bounds[2]
    B_L[i] <- bounds[3]
    B_U[i] <- bounds[4]
  }

  # Compute statistics
  Ave_Bias_a_L <- mean(A_L) - a
  Ave_Bias_a_U <- mean(A_U) - a
  MSE_a_L <- sqrt(mean((A_L - a)^2))
  MSE_a_U <- sqrt(mean((A_U - a)^2))

  Ave_Bias_b_L <- mean(B_L) - b
  Ave_Bias_b_U <- mean(B_U) - b
  MSE_b_L <- sqrt(mean((B_L - b)^2))
  MSE_b_U <- sqrt(mean((B_U - b)^2))

  c <- round(c(mean(A_L), mean(A_U),
               mean(B_L), mean(B_U),
               Ave_Bias_a_L, Ave_Bias_a_U,
               Ave_Bias_b_L, Ave_Bias_b_U,
               MSE_a_L, MSE_a_U,
               MSE_b_L, MSE_b_U), 4)
  return(c)
}

f3 <- Vectorize(f2)

# ------------------------------------------------------------------------
# Run Simulation for Multiple Sample Sizes
# ------------------------------------------------------------------------
set.seed(88888888)
for (n in c(30, 50, 100, 200)) {
  start_time <- Sys.time()
  d <- f3(n)
  cat(paste0("[", d[1], ",", d[2],
             "] & [", d[3], ",", d[4],
             "] & [", d[5], ",", d[6],
             "] & [", d[7], ",", d[8],
             "] & [", d[9], ",", d[10],
             "] & [", d[11], ",", d[12], "]
"))
  end_time <- Sys.time()
  print(end_time - start_time)
}
