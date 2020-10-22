test.df1 <- data.frame(
  id = 'A',
  pval = 1e-04,
  decision.times = 1
)

test.df2 <- data.frame(
  id = c('A','B'),
  pval = c(1e-04,1e-05),
  decision.times = c(1,NA)
)

test.df3 <- data.frame(
  id = 'A',
  pval = 1e-04,
  decision.times = 1.2
)

test.df4 <- data.frame(
  id = 'A',
  pval = 1e-04,
  decision.times = -1
)

test.df5 <- data.frame(
  id = c('A','B'),
  pval = c(1e-04,1e-05),
  lags = c(1,NA)
)

test.df6 <- data.frame(
  id = 'A',
  pval = 1e-04,
  lags = 1.2
)

test.df7 <- data.frame(
  id = 'A',
  pval = 1e-04,
  lags = -1
)

test.df8 <- data.frame(
  id = 'A',
  pval = 1e-04
)

test_that("Check error messages", {
    expect_error(checkStarVersion(test.df1, N=1, version='ver'),
                 "version must be 'async', 'dep' or 'batch'.")
  
    expect_error(checkStarVersion(test.df2, N=1, version='async'),
                 "Please provide a decision time for all p-values.")
    
    expect_error(checkStarVersion(test.df3, N=1, version='async'),
                 "All decision times should be integers.")
    
    expect_error(checkStarVersion(test.df4, N=1, version='async'),
                 "All decision times should be positive integers.")
    
    expect_error(checkStarVersion(test.df5, N=1, version='dep'),
                 "Please provide lags for all p-values.")
    
    expect_error(checkStarVersion(test.df6, N=1, version='dep'),
                 "All lags should be integers.")
    
    expect_error(checkStarVersion(test.df7, N=1, version='dep'),
                 "All lags should be positive integers.")
    
    expect_error(checkStarVersion(test.df8, N=1, version='batch', batch.sizes = 2),
                 "The sum of the batch sizes must equal the number of p-values observed.")
    
    expect_error(checkStarVersion(test.df8, N=3, version='batch', batch.sizes = c(1.2, 1.8)),
                 "All batch sizes should be integers.")
    
    expect_error(checkStarVersion(test.df8, N=4, version='batch', batch.sizes = c(-1, 5)),
                 "All batch sizes should be positive integers.")
})
