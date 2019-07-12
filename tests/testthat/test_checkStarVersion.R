test_that("Check error messages", {
    expect_error(checkStarVersion(N=1, version='ver'),
                 "version must be 'async', 'dep' or 'batch'.")
  
    expect_error(checkStarVersion(N=1, version='async', decision.times=c(1,1)),
                 "Please provide a decision time for all p-values.")
    
    expect_error(checkStarVersion(N=1, version='async', decision.times=1.2),
                 "All decision times should be integers.")
    
    expect_error(checkStarVersion(N=1, version='dep', lags=c(1,1)),
                 "Please provide lags for all p-values.")
    
    expect_error(checkStarVersion(N=1, version='dep', lags=1.2),
                 "All lags should be integers.")
    
    expect_error(checkStarVersion(N=1, version='batch', batch.sizes=2),
                 "The sum of the batch sizes must equal the number of p-values observed.")
    
    expect_error(checkStarVersion(N=3, version='batch', batch.sizes=c(1.2, 1.8)),
                 "All batch sizes should be integers.")
    
    expect_error(checkStarVersion(N=4, version='batch', batch.sizes=c(-1, 5)),
                 "All batch sizes should be positive integers.")
})
