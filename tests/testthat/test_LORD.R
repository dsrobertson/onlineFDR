test.df1 <- data.frame(
    id = 'A',
    date = as.Date("2014-12-01"),
    pval = 1e-04
)

test.df2 <- data.frame(
    id = c('A', 'B', 'C', 'D'),
    date = as.Date(rep("2014-12-01",4)),
    pval = c(1e-07, 0.1, 0.00025, 0.07)
)

test2plus <- LORD(test.df2, version='++', random=FALSE)$R
test2discard <- LORD(test.df2, version='discard', random=FALSE)$R
test2v3 <- LORD(test.df2, version=3, random=FALSE)$R
test2dep <- LORD(test.df2, version='dep', random=FALSE)$R


test_that("Errors for edge cases", {
    expect_error(LORD(matrix(NA, nrow=2, ncol=2)),
                 "d must either be a dataframe or a vector of p-values.")
    
    expect_error(LORD(0.1, alpha = -0.1),
                 "alpha must be between 0 and 1.")
    
    expect_error(LORD(0.1, gammai = -1),
                 "All elements of gammai must be non-negative.")
    
    expect_error(LORD(0.1, gammai=2),
                 "The sum of the elements of gammai must be <= 1.")
    
    expect_error(LORD(0.1, version='dep', gammai = -1),
                 "All elements of gammai must be non-negative.")
    
    expect_error(LORD(0.1, version='dep', gammai=2),
                 "The sum of the elements of gammai must be <= alpha/b0.")
    
    expect_error(LORD(0.1, w0 = -0.01),
                 "w0 must be non-negative.")
    
    expect_error(LORD(0.1, version=3, b0 = -0.01),
                 "b0 must be positive.")
    
    expect_error(LORD(0.1, version=3, alpha=0.05, b0=0.02, w0=0.04),
                 "The sum of w0 and b0 must not be greater than alpha.")
    
    expect_error(LORD(0.1, alpha=0.05, w0=0.06),
                 "w0 must not be greater than alpha.")
})

test_that("LORD gives same results when N = 1", {
    expect_equal(LORD(test.df1), LORD(test.df1, version=3))
})

test_that("Correct rejections for sample data", {
    expect_identical(test2plus, c(1,0,1,0))
    expect_identical(test2discard, c(1,0,1,0))
    expect_identical(test2v3, c(1,0,1,0))
    expect_identical(test2dep, c(1,0,1,0))
    
    expect_identical(LORD(0.1, version='discard')$R, 0)
    expect_identical(LORD(c(0.1,0.1), version='discard')$R, c(0,0))
    
    expect_identical(LORD(0.1, version='dep')$R, 0)
})
