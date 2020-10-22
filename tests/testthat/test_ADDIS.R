test.pval <- c(1e-07, 3e-04, 0.1, 6e-04)
test.df <- data.frame(id = seq_len(4), pval = test.pval)
test.df2 <- data.frame(id = seq_len(4), pval = test.pval, decision.times = c(2,3,4,NA))
test.df3 <- data.frame(id = seq_len(4), pval = test.pval, decision.times = seq_len(4)+3)
test.df4 <- data.frame(id = seq_len(4), pval = test.pval, decision.times = seq_len(4))

test_that("Errors for edge cases", {

  expect_error(ADDIS(test.df, alpha = -0.1),
               "alpha must be between 0 and 1.")
  
  expect_error(ADDIS(test.df, lambda = -0.1),
               "lambda must be between 0 and 1.")
  
  expect_error(ADDIS(test.df, tau = -0.1),
               "tau must be between 0 and 1.")
  
  expect_error(ADDIS(test.df, gammai = -1),
               "All elements of gammai must be non-negative.")
  
  expect_error(ADDIS(test.df, gammai=2),
               "The sum of the elements of gammai must not be greater than 1.")
  
  expect_error(ADDIS(test.df, w0 = -0.01),
               "w0 must be non-negative.")
  
  # expect_error(ADDIS(0.1, w0 = 1),
  #              "w0 must be less than tau*lambda*alpha")
  
  expect_error(ADDIS(test.df2, async=TRUE),
              "Please provide a decision time for each p-value.")
})


test_that("Correct rejections", {
  
    expect_identical(ADDIS(test.df)$R, c(1,1,0,1))
    
    expect_identical(ADDIS(test.df3,
                           async=TRUE)$R,
                           c(1,1,0,0))
})

test_that("Check that ADDIS with async=FALSE is a special case of async=TRUE", {
              expect_equal(ADDIS(test.pval)$alphai,
                           ADDIS(test.df4, async=TRUE)$alphai)
})

test_that("Check that ADDIS gives same results as SAFFRON and SAFFRONstar
          algorithms with discard=TRUE", {
              expect_equal(ADDIS(test.pval)$alphai,
                           SAFFRON(test.df, discard=TRUE)$alphai)
              
              expect_equal(ADDIS(test.df3, async=TRUE)$alphai,
                           SAFFRONstar(test.df3, version='async', discard=TRUE)$alphai)
})

test_that("ADDIS inputs are correct with async=TRUE", {
  expect_error(ADDIS(test.df, async=TRUE),
               "d needs to have a column of decision.times")
  expect_error(ADDIS(c(0.1, 0.1), async=TRUE),
               "d needs to be a dataframe with a column of decision.times")
})
