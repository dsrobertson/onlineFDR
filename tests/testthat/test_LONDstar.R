test.pval <- c(1e-07, 3e-04, 0.1, 5e-04)
test.df <- data.frame(id = seq_len(4), pval = test.pval)


test_that("Errors for edge cases", {
    
    expect_error(LONDstar(0.1, version='async', decision.times=1, alpha = -0.1),
                 "alpha must be between 0 and 1.")
    
    expect_error(LONDstar(0.1, version='async', decision.times=1, betai = -1),
                 "All elements of betai must be non-negative.")
    
    expect_error(LONDstar(0.1, version='async', decision.times=1, betai=2),
                 "The sum of the elements of betai must be <= alpha.")

})

test_that("Correct rejections for version async", {
    
    expect_identical(LONDstar(0.1, version='async', decision.times=1)$R, 0)
    
    expect_identical(LONDstar(test.pval, version='async',
                              decision.times = seq_len(4))$R,
                     c(1,1,0,1))
    
    expect_identical(LONDstar(test.pval, version='async',
                              decision.times = rep(4,4))$R,
                     c(1,1,0,0))
})

test_that("Correct rejections for version dep", {
    
    expect_identical(LONDstar(0.1, version='dep', lags=1)$R, 0)
    
    expect_identical(LONDstar(test.pval, version='dep',
                              lags = rep(0,4))$R,
                     c(1,1,0,1))
    
    expect_identical(LONDstar(test.pval, version='dep',
                              lags = rep(4,4))$R,
                     c(1,1,0,0))
})


test_that("Correct rejections for version batch", {
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
                 LONDstar(test.pval, version='async',
                          decision.times = seq_len(4))$alphai)
              
    expect_equal(LOND(test.df, original=FALSE)$alphai,
                 LONDstar(test.pval, version='dep',
                          lags = rep(0,4))$alphai)
              
    expect_equal(LOND(test.df, original=FALSE)$alphai,
                 LONDstar(test.pval, version='batch',
                          batch.sizes = rep(1,4))$alphai)              
})
