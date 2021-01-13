test.df1 <- data.frame(
  id = 'A',
  pval = 1e-04,
  batch = 1
)

test.df2 <- data.frame(
  id = c('A', 'B', 'C'),
  pval = c(1e-04, 0.1, 1e-04),
  batch = c(1,1,2)
)

test_that("Errors for edge cases", {
  expect_error(BatchPRDS(matrix(NA, nrow=2, ncol=2)),
               "d must be a dataframe")
  
  expect_error(BatchPRDS(rep(0.01, 4)),
               "d must be a dataframe")
  
  expect_error(BatchPRDS(test.df1, alpha=2),
               "alpha must be between 0 and 1.")
  
  expect_error(BatchPRDS(test.df1, gammai=-1),
               "All elements of gammai must be non-negative.")
  
  expect_error(BatchPRDS(test.df1, gammai=2),
               "The sum of the elements of gammai must not be greater than 1.")
})

test_that("Correct rejections for sample dataframes", {
  expect_identical(BatchPRDS(test.df1)$R, 1)
  expect_identical(BatchPRDS(test.df2)$R, c(1,0,1))
})