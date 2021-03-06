test.pval <- c(1e-07, 3e-04, 0.1, 6e-04)
test.df <- data.frame(id = seq_len(4), pval = test.pval, lags = rep(3,4))
test.df2 <- data.frame(id = seq_len(4), pval = test.pval, lags = rep(0,4))

test_that("Errors for edge cases", {

  expect_error(ADDIS_spending(test.df, alpha = -0.1),
               "alpha must be between 0 and 1.")
  
  expect_error(ADDIS_spending(test.df, lambda = -0.1),
               "lambda must be between 0 and 1.")
  
  expect_error(ADDIS_spending(test.df, tau = -0.1),
               "tau must be between 0 and 1.")
  
  expect_error(ADDIS_spending(test.df, gammai = -1),
               "All elements of gammai must be non-negative.")
  
  expect_error(ADDIS_spending(test.df, gammai=2),
               "The sum of the elements of gammai must not be greater than 1.")
  
  expect_error(ADDIS_spending(test.df, lambda=0.5, tau=0.25),
              "lambda must be less than tau.")
})

test_that("Correct rejections", {
    
    expect_identical(ADDIS_spending(c(0.1,0.1))$R, c(0,0))
  
    expect_identical(ADDIS_spending(test.pval)$R, c(1,1,0,1))
    
    expect_identical(ADDIS_spending(test.df, dep=TRUE)$R,
                     c(1,1,0,0))
})

test_that("Check that ADDIS-spending with dep=FALSE is a special case of dep=TRUE", {
              expect_equal(ADDIS_spending(test.pval)$alphai,
                           ADDIS_spending(test.df2, dep=TRUE)$alphai)
})
