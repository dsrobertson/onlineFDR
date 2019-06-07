test.df0 <- data.frame()
test.df1 <- data.frame(id = 'A')

test.df2 <- data.frame(
    id = 'A',
    date = as.Date("2014-12-01")
)

test.df3 <- data.frame(
    id = 'A',
    date = as.Date("2014-12-01"),
    pval = 1e-04
)

test.df4 <- data.frame(
    id = c('A', 'B', 'C'),
    date = as.Date(rep("2014-12-01",3)),
    pval = c(1e-07, 0.1, 1e-07)
)

test4 <- SAFFRON(test.df4, random=FALSE)$R


test_that("Errors for edge cases", {
    expect_error(SAFFRON(1e-04), "d must be a dataframe.")
    expect_error(SAFFRON(test.df0), "The dataframe d is missing a column 'id' of identifiers.")
    expect_error(SAFFRON(test.df1), "The dataframe d is missing a column 'pval' of p-values.")
    expect_error(SAFFRON(test.df2), "The dataframe d is missing a column 'pval' of p-values.")
})

test_that("Correct rejections for sample dataframes", {
    expect_identical(SAFFRON(test.df3)$R, 1)
    expect_identical(test4, c(1,0,1))
})
