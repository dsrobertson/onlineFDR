# Check if our tolerance fixes are causing mathematical issues

# Test the specific case from vignette
alpha <- 0.04
N <- 20

# Test different tolerance multipliers
multipliers <- c(1, 10, 100, 1000)

cat("Testing tolerance effects:\n\n")

for (mult in multipliers) {
  # Create a bound that exceeds alpha slightly (simulate worst case)
  bound <- rep(alpha/N, N)
  # Add small error to simulate platform differences
  bound[1] <- bound[1] + 1e-15
  
  sum_bound <- sum(bound)
  diff <- sum_bound - alpha
  tolerance <- mult * .Machine$double.eps * N
  
  cat("Multiplier:", mult, "\n")
  cat("  Sum of bound:", sprintf("%.20f", sum_bound), "\n")
  cat("  Alpha:       ", sprintf("%.20f", alpha), "\n")
  cat("  Difference:  ", sprintf("%.2e", diff), "\n")
  cat("  Tolerance:   ", sprintf("%.2e", tolerance), "\n")
  cat("  Passes check:", diff <= tolerance, "\n")
  cat("  Mathematical validity:", sum_bound <= alpha + 0.001, "\n") # reasonable bound
  cat("\n")
}

# Test if 100x tolerance might be allowing invalid cases
cat("=== Checking if 100x tolerance is too permissive ===\n")

# Simulate a case where betai genuinely exceeds alpha
bad_bound <- rep(alpha/N * 1.01, N)  # 1% over
sum_bad <- sum(bad_bound)
diff_bad <- sum_bad - alpha
tolerance_100x <- 100 * .Machine$double.eps * N

cat("Intentionally bad bound (1% over alpha):\n")
cat("Sum:", sprintf("%.6f", sum_bad), "\n")
cat("Alpha:", sprintf("%.6f", alpha), "\n")
cat("Difference:", sprintf("%.6f", diff_bad), "\n")
cat("100x tolerance:", sprintf("%.2e", tolerance_100x), "\n")
cat("Would bad case pass 100x tolerance?", diff_bad <= tolerance_100x, "\n")
