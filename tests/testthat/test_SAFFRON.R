test.df1 <- data.frame(
    id = 'A',
    date = as.Date("2014-12-01"),
    pval = 1e-04
)

test.df2 <- data.frame(
    id = c('A', 'B', 'C'),
    date = as.Date(rep("2014-12-01",3)),
    pval = c(1e-07, 0.1, 1e-07)
)

test2 <- SAFFRON(test.df2, random=FALSE)$R


test_that("Errors for edge cases", {
   
    expect_error(SAFFRON(matrix(NA, nrow=2, ncol=2)),
                 "d must either be a dataframe or a vector of p-values.")
    
    expect_error(SAFFRON(0.1, alpha = -0.1),
                 "alpha must be between 0 and 1.")
    
    expect_error(SAFFRON(0.1, lambda=2),
                 "lambda must be between 0 and 1.")
    
    expect_error(SAFFRON(0.1, gammai = -1),
                 "All elements of gammai must be non-negative.")
    
    expect_error(SAFFRON(0.1, gammai=2),
    "The sum of the elements of gammai must not be greater than 1.")
    
    expect_error(SAFFRON(0.1, w0 = -0.01),
                 "w0 must be non-negative.")
})

test_that("Correct rejections for sample dataframes", {
    expect_identical(SAFFRON(test.df1)$R, 1)
    expect_identical(test2, c(1,0,1))
    
    expect_identical(SAFFRON(c(0.1,0.1))$R, c(0,0))
})
