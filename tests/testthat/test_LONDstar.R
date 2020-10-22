test.pval <- c(1e-07, 3e-04, 0.1, 5e-04)
test.decision.times <- seq_len(4)
test.decision.times2 <- rep(4,4)
test.lags <- rep(0,4)
test.lags2 <- rep(4,4)
test.batch.sizes <- rep(1,4)
test.batch.sizes <- 4
test.df <- data.frame(id = seq_len(4), pval = test.pval, decision.times = test.decision.times)
test.df2 <- data.frame(id = seq_len(4), pval = test.pval, decision.times = test.decision.times2)
test.df3 <-  data.frame(id = seq_len(4), pval = test.pval, lags = test.lags)
test.df4 <-  data.frame(id = seq_len(4), pval = test.pval, lags = test.lags2)

test_that("Errors for edge cases", {
    
    expect_error(LONDstar(test.df, alpha = -0.1, version='async'),
                 "alpha must be between 0 and 1.")
    
    expect_error(LONDstar(test.df, betai = -1, version='async'),
                 "All elements of betai must be non-negative.")
    
    expect_error(LONDstar(test.df, betai = 2, version='async'),
                 "The sum of the elements of betai must be <= alpha.")

})

test_that("Correct rejections for version async", {
    
    expect_identical(LONDstar(test.df, version='async')$R,
                     c(1,1,0,1))
    
    expect_identical(LONDstar(test.df2, version='async')$R,
                     c(1,1,0,0))
})

test_that("Correct rejections for version dep", {
    
    expect_identical(LONDstar(test.df3, version='dep')$R,
                     c(1,1,0,1))
    
    expect_identical(LONDstar(test.df4, version='dep')$R,
                     c(1,1,0,0))
})


test_that("Correct rejections for version batch", {
    
    expect_identical(LONDstar(0.1, version='batch', batch.sizes=1)$R, 0)
    
    expect_identical(LONDstar(test.pval, version='batch',
                              batch.sizes = rep(1,4))$R,
                     c(1,1,0,1))
    
    expect_identical(LONDstar(test.pval, version='batch',
                              batch.sizes = 4)$R,
                     c(1,1,0,0))
})


test_that("Check that LOND is a special case of the LONDstar
          algorithms", {
    expect_equal(LOND(test.df, original=FALSE)$alphai,
                 LONDstar(test.df, version='async')$alphai)
              
    expect_equal(LOND(test.df, original=FALSE)$alphai,
                 LONDstar(test.df3, version='dep')$alphai)
              
    expect_equal(LOND(test.df, original=FALSE)$alphai,
                 LONDstar(test.pval, version='batch',
                          batch.sizes = rep(1,4))$alphai)              
})

test_that("LONDstar inputs are correct when given vector input", {
    expect_error(LONDstar(test.pval, version="async"),
                 "d needs to be a dataframe with a column of decision.times")
    expect_error(LONDstar(test.pval, version="dep"),
                 "d needs to be a dataframe with a column of lags")
})
