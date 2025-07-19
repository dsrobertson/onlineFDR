# Test script to verify the tolerance fix works
suppressPackageStartupMessages({
  library(devtools)
  load_all()
})

# Create test data similar to the vignette
sample.df <- data.frame(
  id = c('A15432', 'B90969', 'C18705', 'B49731', 'E99902',
         'C38292', 'A30619', 'D46627', 'E29198', 'A41418',
         'D51456', 'C88669', 'E03673', 'A63155', 'B66033'),
  date = as.Date(c(rep("2014-12-01",3),
                   rep("2015-09-21",5),
                   rep("2016-05-19",2),
                   "2016-11-12",
                   rep("2017-03-27",4))),
  pval = c(2.90e-14, 0.06743, 0.01514, 0.08174, 0.00171,
           3.61e-05, 0.79149, 0.27201, 0.28295, 7.59e-08,
           0.69274, 0.30443, 0.000487, 0.72342, 0.54757))

# Test the problematic case
alpha <- 0.04
N <- 20
bound <- setBound("LOND", alpha = alpha, N)

cat("Testing LOND with alpha =", alpha, "and N =", N, "\n")
cat("Sum of bound:", sum(bound), "\n")
cat("Alpha:", alpha, "\n")
cat("Difference:", sum(bound) - alpha, "\n")
cat("Old tolerance:", .Machine$double.eps * length(bound), "\n")
cat("New tolerance (100x):", 100 * .Machine$double.eps * length(bound), "\n")
cat("Is difference < new tolerance?", abs(sum(bound) - alpha) < 100 * .Machine$double.eps * length(bound), "\n\n")

# Test if LOND runs without error
set.seed(1)
tryCatch({
  result <- LOND(sample.df, alpha = alpha, betai = bound)
  cat("SUCCESS: LOND ran without error!\n")
  cat("Number of rejections:", sum(result$R), "\n")
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
})

cat("\n=== Additional tests ===\n")

# Test edge cases
test_cases <- list(
  list(alpha = 0.05, N = 10),
  list(alpha = 0.01, N = 100),
  list(alpha = 0.1, N = 5)
)

for (i in seq_along(test_cases)) {
  tc <- test_cases[[i]]
  cat("Test case", i, ": alpha =", tc$alpha, ", N =", tc$N, "\n")
  
  bound_test <- setBound("LOND", alpha = tc$alpha, tc$N)
  diff <- abs(sum(bound_test) - tc$alpha)
  tolerance <- 100 * .Machine$double.eps * length(bound_test)
  
  cat("  Difference:", diff, "\n")
  cat("  Tolerance:", tolerance, "\n")
  cat("  OK:", diff < tolerance, "\n\n")
}
