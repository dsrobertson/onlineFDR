test.df0 <- data.frame()
test.df1 <- data.frame(id = 'A')

test.df2 <- data.frame(
    id = 'A',
    date = as.Date("2014-12-01")
)

test.df3 <- data.frame(
    id = 'A',
    date = "date",
    pval = 0.01
)

test_that("Check error messages", {
    
    expect_error(checkdf(test.df0),
                 "The dataframe d is missing a column 'id' of identifiers.")
    
    expect_error(checkdf(test.df1),
                 "The dataframe d is missing a column 'pval' of p-values.")
    
    expect_error(checkdf(test.df2),
                 "The dataframe d is missing a column 'pval' of p-values.")
    
    expect_error(checkdf(test.df3, date.format="%Y-%m-%d"),
                 "One or more dates are not in the correct format.")
})
