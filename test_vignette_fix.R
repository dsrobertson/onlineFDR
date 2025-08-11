# Test the exact vignette scenario with our fix
suppressPackageStartupMessages({
  library(devtools)
  load_all()
})

# Test data from vignette
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

# The exact failing line from vignette
bound <- setBound("LOND", alpha = 0.04, 20)

cat("setBound result:\n")
cat("Sum of bound:", sum(bound), "\n")
cat("Alpha:", 0.04, "\n")
cat("Exactly equal:", sum(bound) == 0.04, "\n\n")

# Test LOND function
set.seed(1)
tryCatch({
  result <- LOND(sample.df, alpha = 0.04, betai = bound)
  cat("SUCCESS: LOND function executed without error!\n")
  cat("Number of discoveries:", sum(result$R), "\n")
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
})
