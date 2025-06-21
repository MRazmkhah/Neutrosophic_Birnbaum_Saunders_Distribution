# ------------------------------------------------------------------------
# Statistical Properties of the Neutrosophic Birnbaum–Saunders Distribution
# ------------------------------------------------------------------------

# ------------------------------------------------------------------------
# Functions for Statistical Properties
# ------------------------------------------------------------------------

EX <- function(x) {
  a <- x[1]
  b <- x[2]
  b / 2 * (a^2 + 2)
}

VarX <- function(x) {
  a <- x[1]
  b <- x[2]
  (b^2 / 4) * (5 * a^4 + 4 * a^2)
}

CV <- function(x) {
  a <- x[1]
  b <- x[2]
  sqrt(5 * a^4 + 4 * a^2) / (a^2 + 2)
}

sk <- function(x) {
  a <- x[1]
  b <- x[2]
  (44 * a^3 + 24 * a) / (5 * a^2 + 4)^1.5
}

ku <- function(x) {
  a <- x[1]
  b <- x[2]
  3 + (558 * a^4 + 240 * a^2) / (5 * a^2 + 4)^2
}

# ------------------------------------------------------------------------
# Parameter Boundaries (Lower and Upper bounds for a and b)
# ------------------------------------------------------------------------

A <- matrix(c(
  0.1, 0.35, 0.5, 0.75, 1, 1.5, 2, 3,
  0.1, 0.35, 0.5, 0.75, 1, 1.5, 2, 3,
  0.1, 0.35, 0.5, 0.75, 1, 1.5, 2, 3
), ncol = 2, byrow = TRUE)

B <- matrix(c(
  0.5, 1, 0.5, 1, 0.5, 1, 0.5, 1,
  1, 2, 1, 2, 1, 2, 1, 2,
  2, 3, 2, 3, 2, 3, 2, 3
), ncol = 2, byrow = TRUE)

# ------------------------------------------------------------------------
# Compute Bounds for Mean, Variance, CV, Skewness, and Kurtosis
# ------------------------------------------------------------------------

D <- matrix(NA, nrow = 12, ncol = 1)

for (j in 1:nrow(A)) {
  a_L <- A[j, 1]
  a_U <- A[j, 2]
  b_L <- B[j, 1]
  b_U <- B[j, 2]

  a0 <- (a_L + a_U) / 2
  b0 <- (b_L + b_U) / 2

  # Compute lower and upper bounds using optimization
  mu_L <- round(optim(par = c(a0, b0), fn = EX, lower = c(a_L, b_L), upper = c(a_U, b_U), 
                      method = "L-BFGS-B")$value, 3)
  mu_U <- round(optim(par = c(a0, b0), fn = EX, lower = c(a_L, b_L), upper = c(a_U, b_U), 
                      method = "L-BFGS-B", control = list(fnscale = -1))$value, 3)

  var_L <- round(optim(par = c(a0, b0), fn = VarX, lower = c(a_L, b_L), upper = c(a_U, b_U), 
                       method = "L-BFGS-B")$value, 3)
  var_U <- round(optim(par = c(a0, b0), fn = VarX, lower = c(a_L, b_L), upper = c(a_U, b_U), 
                       method = "L-BFGS-B", control = list(fnscale = -1))$value, 3)

  CV_L <- round(optim(par = c(a0, b0), fn = CV, lower = c(a_L, b_L), upper = c(a_U, b_U), 
                      method = "L-BFGS-B")$value, 3)
  CV_U <- round(optim(par = c(a0, b0), fn = CV, lower = c(a_L, b_L), upper = c(a_U, b_U), 
                      method = "L-BFGS-B", control = list(fnscale = -1))$value, 3)

  skew_L <- round(optim(par = c(a0, b0), fn = sk, lower = c(a_L, b_L), upper = c(a_U, b_U), 
                        method = "L-BFGS-B")$value, 3)
  skew_U <- round(optim(par = c(a0, b0), fn = sk, lower = c(a_L, b_L), upper = c(a_U, b_U), 
                        method = "L-BFGS-B", control = list(fnscale = -1))$value, 3)

  kur_L <- round(optim(par = c(a0, b0), fn = ku, lower = c(a_L, b_L), upper = c(a_U, b_U), 
                       method = "L-BFGS-B")$value, 3)
  kur_U <- round(optim(par = c(a0, b0), fn = ku, lower = c(a_L, b_L), upper = c(a_U, b_U), 
                       method = "L-BFGS-B", control = list(fnscale = -1))$value, 3)

  D[j, ] <- paste0("[", b_L, ", ", b_U, "] & [", a_L, ", ", a_U, "] & ",
                   "[", mu_L, ", ", mu_U, "] & [", var_L, ", ", var_U, "] & ",
                   "[", CV_L, ", ", CV_U, "] & [", skew_L, ", ", skew_U, "] & ",
                   "[", kur_L, ", ", kur_U, "]")
}

# ------------------------------------------------------------------------
# Display Results
# ------------------------------------------------------------------------
print(D)