test.df1 <- data.frame(
    id = 'A',
    date = as.Date("2014-12-01"),
    pval = 1e-04
)

test.df2 <- data.frame(
    id = c('A', 'B', 'C'),
    date = as.Date(rep("2014-12-01",3)),
    pval = c(1e-07, 0.1, 6e-04)
)

test2 <- LOND(test.df2, random=FALSE)$R
test2dep <- LOND(test.df2, dep=TRUE, random=FALSE)$R


test_that("Errors for edge cases", {
    expect_error(LOND(matrix(NA, nrow=2, ncol=2)),
                 "d must either be a dataframe or a vector of p-values.")
    
    expect_error(LOND(0.1, alpha = -0.1),
                 "alpha must be between 0 and 1.")
    
    expect_error(LOND(0.1, betai = -1),
                 "All elements of betai must be non-negative.")
    
    expect_error(LOND(0.1, betai=2),
    "The sum of the elements of betai must not be greater than alpha.")
})

test_that("LOND gives same results when dep = TRUE for N = 1", {
    expect_equal(LOND(test.df1), LOND(test.df1, dep = TRUE))
})

test_that("Correct rejections for sample data", {
    expect_identical(LOND(test.df1)$R, 1)
    expect_identical(test2, c(1,0,1))
    expect_identical(test2dep, c(1,0,0))
})
