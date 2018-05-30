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
    pval = c(1e-07, 0.00055, 0.1)
)


test_that("Errors for edge cases", {
    expect_error(LORD(1e-04), "d must be a dataframe.")
    expect_error(LORD(test.df0), "The dataframe d is missing a column 'id' of identifiers.")
    expect_error(LORD(test.df1), "The dataframe d is missing a column 'pval' of p-values.")
    expect_error(LORD(test.df2), "The dataframe d is missing a column 'pval' of p-values.")
})

test_that("LORD gives same results for all three versions when N = 1", {
    expect_equal(LORD(test.df3), LORD(test.df3, version=2))
    expect_equal(LORD(test.df3, version=2), LORD(test.df3, version=1))
})

test_that("Correct rejections for sample dataframes", {
    expect_identical(LORD(test.df3)$R, 1)
    expect_identical(LORD(test.df4, seed=1)$R, c(1,0,1))
    expect_identical(LORD(test.df4, version=2, seed=1)$R, c(1,0,1))
    expect_identical(LORD(test.df4, version=1, seed=1)$R, c(1,0,0))
})
