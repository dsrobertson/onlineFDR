test.df <- data.frame(
  id = c('A', 'B', 'C'),
  date = as.Date(rep("2014-12-01",3)),
  pval = c(1e-07, 1e-07, 0.1)
)

test_that('Rejections', {
  expect_equal(LORDdep(test.df)$R, c(1,0,1))
})
