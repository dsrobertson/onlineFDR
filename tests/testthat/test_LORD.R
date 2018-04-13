test.df <- data.frame(
  id = c('A', 'B', 'C'),
  date = as.Date(rep("2014-12-01",3)),
  pval = c(1e-07, 0.00055, 0.1)
)

test_that('Rejections', {
  expect_equal(LORD(test.df)$R, c(1,0,1))
  expect_equal(LORD(test.df, version=2)$R, c(1,0,1))
  expect_equal(LORD(test.df, version=1)$R, c(1,0,0))
})
