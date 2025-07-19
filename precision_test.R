# Simple test for numerical precision issue
alpha <- 0.04
N <- 20

# This is what setBound("LOND", alpha, N) does:
bound <- rep(alpha/N, N)

cat("Alpha:", alpha, "\n")
cat("N:", N, "\n")
cat("alpha/N:", alpha/N, "\n")
cat("Sum of bound:", sum(bound), "\n")
cat("Difference:", sum(bound) - alpha, "\n")
cat("Machine epsilon:", .Machine$double.eps, "\n")
cat("Current tolerance:", .Machine$double.eps * length(bound), "\n")
cat("Suggested tolerance (100x):", 100 * .Machine$double.eps * length(bound), "\n")
cat("Is current tolerance sufficient?", abs(sum(bound) - alpha) <= .Machine$double.eps * length(bound), "\n")
cat("Is 100x tolerance sufficient?", abs(sum(bound) - alpha) <= 100 * .Machine$double.eps * length(bound), "\n")

# More robust tolerance test
tolerance_multipliers <- c(1, 10, 100, 1000)
for (mult in tolerance_multipliers) {
  tol <- mult * .Machine$double.eps * length(bound)
  is_ok <- abs(sum(bound) - alpha) <= tol
  cat("Multiplier", mult, "- Tolerance:", tol, "- OK:", is_ok, "\n")
}
