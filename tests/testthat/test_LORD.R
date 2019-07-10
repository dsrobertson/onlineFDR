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
    pval = c(1e-07, 0.1, 0.00025)
)

test4plus <- LORD(test.df4, version='++', random=FALSE)$R
test4discard <- LORD(test.df4, version='discard', random=FALSE)$R
test4v3 <- LORD(test.df4, version=3, random=FALSE)$R
test4dep <- LORD(test.df4, version='dep', random=FALSE)$R


test_that("Errors for edge cases", {
    expect_error(LORD(test.df0), "The dataframe d is missing a column 'id' of identifiers.")
    expect_error(LORD(test.df1), "The dataframe d is missing a column 'pval' of p-values.")
    expect_error(LORD(test.df2), "The dataframe d is missing a column 'pval' of p-values.")
})

test_that("LORD gives same results when N = 1", {
    expect_equal(LORD(test.df3), LORD(test.df3, version=3))
})

test_that("Correct rejections for sample dataframes", {
    expect_identical(test4plus, c(1,0,1))
    expect_identical(test4discard, c(1,0,1))
    expect_identical(test4v3, c(1,0,1))
    expect_identical(test4dep, c(1,0,1))
    
})
