test.pval <- c(1e-07, 3e-04, 0.1, 5e-04)
test.df <- data.frame(id = seq_len(4), pval = test.pval)
test.df2 <- data.frame(id = seq_len(4), pval = test.pval, decision.times = seq_len(4))

test_that("Dataframe does not contain certain columns", {
  expect_error(LONDstar(test.df, version = 'async'),
               "Your dataset does not contain the column 'decision.times'.")
  
  expect_error(LONDstar(test.df2, version = 'dep'),
               "Your dataset does not contain the column 'lags'.")
})
