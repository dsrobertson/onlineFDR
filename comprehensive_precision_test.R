# Test to simulate potential macOS precision issues
# Create the exact scenario from the vignette

# Test different alpha values that might have precision issues on macOS
test_alphas <- c(0.04, 0.05, 0.01, 0.03, 0.1)
test_Ns <- c(10, 20, 100)

cat("Testing setBound LOND precision across different alpha and N values:\n\n")

for (alpha in test_alphas) {
  for (N in test_Ns) {
    # Simulate setBound("LOND", alpha, N)
    bound <- rep(alpha/N, N)
    sum_bound <- sum(bound)
    diff <- sum_bound - alpha
    
    # Test with different tolerance levels
    basic_tol <- .Machine$double.eps * N
    our_tol <- 100 * .Machine$double.eps * N
    
    cat(sprintf("Alpha: %.3f, N: %d\n", alpha, N))
    cat(sprintf("  Sum: %.20f\n", sum_bound))
    cat(sprintf("  Difference: %.2e\n", diff))
    cat(sprintf("  Basic tolerance: %.2e - Pass: %s\n", basic_tol, abs(diff) <= basic_tol))
    cat(sprintf("  Our tolerance:   %.2e - Pass: %s\n", our_tol, abs(diff) <= our_tol))
    cat("\n")
  }
}

# Test with the exact values from the failing vignette
cat("=== EXACT VIGNETTE TEST ===\n")
alpha <- 0.04
N <- 20
bound <- rep(alpha/N, N)
diff <- sum(bound) - alpha

cat("Vignette case (alpha=0.04, N=20):\n")
cat("Difference:", sprintf("%.2e", diff), "\n")
cat("100x tolerance:", sprintf("%.2e", 100 * .Machine$double.eps * N), "\n")
cat("Would pass validation:", abs(diff) <= 100 * .Machine$double.eps * N, "\n")

# Simulate a potential macOS difference (artificially adding small error)
mac_error <- 1e-15  # Simulate potential macOS precision difference
bound_mac <- bound + mac_error / N
diff_mac <- sum(bound_mac) - alpha

cat("\nSimulated macOS precision difference:\n")
cat("Difference with error:", sprintf("%.2e", diff_mac), "\n")
cat("Would pass old validation:", abs(diff_mac) <= .Machine$double.eps * N, "\n")
cat("Would pass new validation:", abs(diff_mac) <= 100 * .Machine$double.eps * N, "\n")
