# ------------------------------------------------------------------------
# PDF Plot for the Neutrosophic Birnbaum-Saunders Distribution (NBSD)
# ------------------------------------------------------------------------

# -----------------------------
# PDF Function for NBSD
# -----------------------------
f_N <- function(t, a, b) { 
  pdf <- (1 / sqrt(2 * pi)) * 
    exp(-1 / (2 * (a^2)) * (t / b + b / t - 2)) * 
    (t^(-3/2) * (t + b)) / (2 * a * sqrt(b))
  return(pdf)
}

# -----------------------------
# Define Input Parameters
# -----------------------------
t <- seq(0, 3, 0.01)    # Grid of t values (BS variable)
a_L <- 0.25             # Lower bound of alpha
a_U <- 0.25             # Upper bound of alpha
b_L <- 0.5              # Lower bound of beta
b_U <- 1.0              # Upper bound of beta

# -----------------------------
# Create Plot Window
# -----------------------------
windows(3, 3)           # (Use quartz() on Mac, x11() on Linux if needed)

# -----------------------------
# Plot Initial PDF Curve (Lower Bound)
# -----------------------------
y1 <- f_N(t, a_L, b_L)
plot(t, y1, type = "l", 
     ylab = "PDF", xlab = expression(paste(t[N])), 
     lwd = 1, ylim = c(0, 5))
axis(side = 1, lwd = 3)
axis(side = 2, lwd = 3)

# -----------------------------
# Fill Neutrosophic Uncertainty Area
# -----------------------------
for (b in seq(b_L, b_U, 0.01)) {
  for (a in seq(a_L, a_U, 0.01)) {
    y2 <- f_N(t, a, b)
    polygon(c(t, rev(t)), c(y2, rev(y1)), col = "#6BD7AF", border = NA)
    y1 <- y2  # Update for next layer
  }
}

# -----------------------------
# Plot Outer Boundaries (Final Curves)
# -----------------------------
y3 <- f_N(t, a_U, b_U)
polygon(c(t, rev(t)), c(y3, rev(y1)), lwd = 1, col = "#6BD7AF")

y4 <- f_N(t, a_L, b_L)
polygon(c(t, rev(t)), c(y4, rev(y4)), lwd = 1, col = "#6BD7AF")
