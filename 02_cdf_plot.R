# ------------------------------------------------------------------------
# CDF Plot for the Neutrosophic Birnbaum-Saunders Distribution (NBSD)
# ------------------------------------------------------------------------

# -----------------------------
# CDF Function for NBSD
# -----------------------------
F_N <- function(t, a, b) {
  x <- (1 / a) * (sqrt(t / b) - sqrt(b / t))
  cdf <- pnorm(x, mean = 0, sd = 1)
  return(cdf)
}

# -----------------------------
# Define Input Parameters
# -----------------------------
t <- seq(0, 3, 0.01)     # Grid of t values (BS variable)
a_L <- 0.1               # Lower bound of alpha
a_U <- 0.35              # Upper bound of alpha
b_L <- 1.0               # Lower bound of beta
b_U <- 1.0               # Upper bound of beta

# -----------------------------
# Open Plot Window
# -----------------------------
windows(3, 3)            # For Windows OS (use quartz() or x11() if needed)

# -----------------------------
# Initial CDF Plot at Lower Bound
# -----------------------------
y1 <- F_N(t, a_L, b_L)

# -----------------------------
# Determine Max Y-Axis Limit for Plotting
# -----------------------------
ff1 <- function(t) { F_N(t, a_L, b_L) }
ff2 <- function(t) { F_N(t, a_U, b_U) }
ans1 <- optimize(ff1, interval = c(0, 3), maximum = TRUE)
ans2 <- optimize(ff2, interval = c(0, 3), maximum = TRUE)
L <- max(ans1$objective, ans2$objective)  # Maximum CDF value

# -----------------------------
# Create Initial Plot
# -----------------------------
plot(t, y1, type = "l", 
     ylab = "CDF", xlab = expression(paste(t[N])), 
     ylim = c(0, L))
axis(side = 1, lwd = 3)
axis(side = 2, lwd = 3)

# -----------------------------
# Fill Neutrosophic Area Between Curves
# -----------------------------
for (b in seq(b_L, b_U, 0.05)) {
  for (a in seq(a_L, a_U, 0.05)) {
    y2 <- F_N(t, a, b)
    polygon(c(t, rev(t)), c(y2, rev(y1)), col = "lightblue", border = NA)
    y1 <- y2  # Update y1 for next iteration
  }
}

# -----------------------------
# Draw Outer Boundaries
# -----------------------------
y3 <- F_N(t, a_U, b_U)
polygon(c(t, rev(t)), c(y3, rev(y3)), lwd = 1, col = "lightblue")

y4 <- F_N(t, a_L, b_L)
polygon(c(t, rev(t)), c(y4, rev(y4)), lwd = 1, col = "lightblue")
