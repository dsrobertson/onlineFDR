test.pval <- c(1e-07, 3e-04, 0.1, 6e-04)
test.df <- data.frame(id = seq_len(4), pval = test.pval)

test_that("Errors for edge cases", {

  expect_error(ADDIS(0.1, alpha = -0.1),
               "alpha must be between 0 and 1.")
  
  expect_error(ADDIS(0.1, lambda = -0.1),
               "lambda must be between 0 and 1.")
  
  expect_error(ADDIS(0.1, tau = -0.1),
               "tau must be between 0 and 1.")
  
  expect_error(ADDIS(0.1, gammai = -1),
               "All elements of gammai must be non-negative.")
  
  expect_error(ADDIS(0.1, gammai=2),
               "The sum of the elements of gammai must not be greater than 1.")
  
  expect_error(ADDIS(0.1, w0 = -0.01),
               "w0 must be non-negative.")
  
  # expect_error(ADDIS(0.1, w0 = 1),
  #              "w0 must be less than tau*lambda*alpha")
  
  expect_error(ADDIS(c(0.1,0.1), async=TRUE, decision.times=1),
              "Please provide a decision time for each p-value.")
})


test_that("Correct rejections", {
    
    expect_identical(ADDIS(0.1)$R, 0)
    expect_identical(ADDIS(0.1, async=TRUE, decision.times=1)$R, 0)
    
    expect_identical(ADDIS(c(0.1,0.1))$R, c(0,0))
  
    expect_identical(ADDIS(test.pval)$R, c(1,1,0,1))
    
    expect_identical(ADDIS(test.pval,
                           async=TRUE, decision.times=seq_len(4)+3)$R,
                     c(1,1,0,0))
})

test_that("Check that ADDIS with async=FALSE is a special case of async=TRUE", {
              expect_equal(ADDIS(test.pval)$alphai,
                           ADDIS(test.pval, async=TRUE,
                                 decision.times=seq_len(4))$alphai)
})

test_that("Check that ADDIS gives same results as SAFFRON and SAFFRONstar
          algorithms with discard=TRUE", {
              expect_equal(ADDIS(test.pval)$alphai,
                           SAFFRON(test.df, discard=TRUE)$alphai)
              
              expect_equal(ADDIS(test.pval, async=TRUE,
                                 decision.times=seq_len(4)+2)$alphai,
                           SAFFRONstar(test.pval, version='async', discard=TRUE,
                                       decision.times=seq_len(4)+2)$alphai)
})
