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

test2 <- supLORD(test.df2, eps = 0.15, r = 2, eta = 0.15, rho = 2, 
                 random=FALSE)$R

test_that("Errors for edge cases", {
    expect_error(supLORD(matrix(NA, nrow=2, ncol=2)),
                 "d must either be a dataframe or a vector of p-values.")

    expect_error(supLORD(0.1, gammai = -1),
                 "All elements of gammai must be non-negative.")
    
    expect_error(supLORD(0.1, gammai=2),
                 "The sum of the elements of gammai must be <= 1.")
    
    expect_error(supLORD(0.1, delta = -0.01),
                 "delta must be between 0 and 1.")
    
    expect_error(supLORD(0.1, eps = -0.01),
                 "eps must be between 0 and 1.")
    
    expect_error(supLORD(0.1, eps = 0.15, r = 2.1),
                 "r must be a positive integer")
    
    expect_error(supLORD(0.1, eps = 0.15, r = 2, eta = -0.01),
                 "eta must be a positive real number.")
    
    expect_error(supLORD(0.1, eps = 0.15, r = 2, eta = 0.15, rho=2.1),
                 "rho must be a positive integer.")
    
})

test_that("Correct rejections for sample data", {
    expect_identical(test2, c(1,0,0,0))
})
