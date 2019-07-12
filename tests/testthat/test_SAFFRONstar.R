test.pval <- c(1e-07, 3e-04, 0.1, 6e-04)
test.df <- data.frame(id = seq_len(4), pval = test.pval)


test_that("Errors for edge cases", {
    
    expect_error(SAFFRONstar(0.1, version='async', decision.times=1, alpha = -0.1),
                 "alpha must be between 0 and 1.")
    
    expect_error(SAFFRONstar(0.1, version='async', decision.times=1, lambda=2),
                 "lambda must be between 0 and 1.")
    
    
    expect_error(SAFFRONstar(0.1, version='async', decision.times=1, gammai = -1),
                 "All elements of gammai must be non-negative.")
    
    expect_error(SAFFRONstar(0.1, version='async', decision.times=1, gammai=2),
                 "The sum of the elements of gammai must not be greater than 1.")
    
    expect_error(SAFFRONstar(0.1, version='async', decision.times=1, w0 = -0.01),
                 "w0 must be non-negative.")

})

test_that("Correct rejections for version async", {
    
    expect_identical(SAFFRONstar(0.1, version='async', decision.times=1)$R, 0)
    
    expect_identical(SAFFRONstar(test.pval, version='async',
                              decision.times = seq_len(4))$R,
                     c(1,1,0,1))
    
    expect_identical(SAFFRONstar(test.pval, version='async',
                              decision.times = rep(4,4))$R,
                     c(1,1,0,0))
})

test_that("Correct rejections for version dep", {
    
    expect_identical(SAFFRONstar(0.1, version='dep', lags=1)$R, 0)
    
    expect_identical(SAFFRONstar(test.pval, version='dep',
                              lags = rep(0,4))$R,
                     c(1,1,0,1))
    
    expect_identical(SAFFRONstar(test.pval, version='dep',
                              lags = rep(4,4))$R,
                     c(1,1,0,0))
})


test_that("Correct rejections for version batch", {
    
    expect_identical(SAFFRONstar(0.1, version='batch', batch.sizes=1)$R, 0)
    
    expect_identical(SAFFRONstar(c(0.1,0.1), version='batch',
                                 batch.sizes=c(1,1))$R, c(0,0))
    
    expect_identical(SAFFRONstar(test.pval, version='batch',
                              batch.sizes = rep(1,4))$R,
                     c(1,1,0,1))
    
    expect_identical(SAFFRONstar(test.pval, version='batch',
                              batch.sizes = 4)$R,
                     c(1,1,0,0))
})


test_that("Check that SAFFRON is a special case of the SAFFRONstar
          algorithms", {
              expect_equal(SAFFRON(test.df)$alphai,
                           SAFFRONstar(test.pval, version='async',
                                    decision.times = seq_len(4))$alphai)
              
              expect_equal(SAFFRON(test.df)$alphai,
                           SAFFRONstar(test.pval, version='dep',
                                    lags = rep(0,4))$alphai)
              
              expect_equal(SAFFRON(test.df)$alphai,
                           SAFFRONstar(test.pval, version='batch',
                                    batch.sizes = rep(1,4))$alphai)              
})
