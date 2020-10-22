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
    
    expect_error(SAFFRONstar(test.df, version='async', alpha = -0.1),
                 "alpha must be between 0 and 1.")
    
    expect_error(SAFFRONstar(test.df, version='async', lambda=2),
                 "lambda must be between 0 and 1.")
    
    expect_error(SAFFRONstar(test.df, version='async', gammai = -1),
                 "All elements of gammai must be non-negative.")
    
    expect_error(SAFFRONstar(test.df, version='async', gammai=2),
                 "The sum of the elements of gammai must not be greater than 1.")
    
    expect_error(SAFFRONstar(test.df, version='async', w0 = -0.01),
                 "w0 must be non-negative.")
    
    # expect_error(SAFFRONstar(0.1, version='async', decision.times=1, w0 = 1),
    #              "w0 must be less than (1-lambda)*alpha")

})

test_that("Correct rejections for version async", {
    
    expect_identical(SAFFRONstar(test.df, version='async')$R,
                     c(1,1,0,1))
    
    expect_identical(SAFFRONstar(test.df2, version='async')$R,
                     c(1,1,0,1))
})

test_that("Correct rejections for version dep", {
    
    expect_identical(SAFFRONstar(test.df3, version='dep')$R,
                     c(1,1,0,1))
    
    expect_identical(SAFFRONstar(test.df4, version='dep')$R,
                     c(1,1,0,1))
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
                     c(1,1,0,1))
})


test_that("Check that SAFFRON is a special case of the SAFFRONstar
          algorithms", {
              expect_equal(SAFFRON(test.df)$alphai,
                           SAFFRONstar(test.df, version='async')$alphai)
              
              expect_equal(SAFFRON(test.df)$alphai,
                           SAFFRONstar(test.df3, version='dep')$alphai)
              
              expect_equal(SAFFRON(test.df)$alphai,
                           SAFFRONstar(test.pval, version='batch',
                                    batch.sizes = rep(1,4))$alphai)              
})

test_that("SAFFRONstar inputs are correct when given vector input", {
    expect_error(SAFFRONstar(test.pval, version="async"),
                 "d needs to be a dataframe with a column of decision.times")
    expect_error(SAFFRONstar(test.pval, version="dep"),
                 "d needs to be a dataframe with a column of lags")
})
