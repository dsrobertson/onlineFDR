# Test the fixed setBound function
source("R/setBound.R")

cat("Testing fixed setBound function:\n\n")

# Test the problematic case
alpha <- 0.04
N <- 20

# Test original vs fixed approach
old_bound <- rep(alpha/N, N)
new_bound <- setBound("LOND", alpha = alpha, N)

cat("Original approach:\n")
cat("  Sum:", sprintf("%.20f", sum(old_bound)), "\n")
cat("  Alpha:", sprintf("%.20f", alpha), "\n")
cat("  Difference:", sprintf("%.2e", sum(old_bound) - alpha), "\n")
cat("  Exactly equal:", sum(old_bound) == alpha, "\n\n")

cat("Fixed approach:\n")
cat("  Sum:", sprintf("%.20f", sum(new_bound)), "\n")
cat("  Alpha:", sprintf("%.20f", alpha), "\n")
cat("  Difference:", sprintf("%.2e", sum(new_bound) - alpha), "\n")
cat("  Exactly equal:", sum(new_bound) == alpha, "\n\n")

# Test different cases
test_cases <- list(
  list(alpha = 0.05, N = 10),
  list(alpha = 0.01, N = 100),
  list(alpha = 0.1, N = 5)
)

cat("Testing various alpha/N combinations:\n")
for (i in seq_along(test_cases)) {
  tc <- test_cases[[i]]
  bound <- setBound("LOND", alpha = tc$alpha, tc$N)
  
  cat("Case", i, "- alpha:", tc$alpha, ", N:", tc$N, "\n")
  cat("  Sum exactly equals alpha:", sum(bound) == tc$alpha, "\n")
  cat("  Difference:", sprintf("%.2e", sum(bound) - tc$alpha), "\n\n")
}
