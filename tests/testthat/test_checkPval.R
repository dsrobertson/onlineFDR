test_that("Check error messages", {
    expect_warning(checkPval(c(0.1,NA)), "Missing p-values were ignored")
    
    expect_error(checkPval(c(0.1,'A')),
                 "The vector of p-values contain at least one non-numeric element.")
    
    expect_error(checkPval(c(0.1,2)), "All p-values must be between 0 and 1.")
})
