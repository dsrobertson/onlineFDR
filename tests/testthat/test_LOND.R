test.df <- data.frame(
  id = c('A', 'B', 'C'),
  date = as.Date(rep("2014-12-01",3)),
  pval = c(1e-07, 6e-04, 0.1)
)

test_that('Rejections', {
  expect_equal(LOND(test.df)$R, c(1,0,1))
  expect_equal(LOND(test.df, dep=TRUE)$R, c(1,0,0))
})
