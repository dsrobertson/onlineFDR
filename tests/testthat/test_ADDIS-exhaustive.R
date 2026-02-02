test_that("ADDIS_exhaustive returns expected structure for numeric vector", {
  set.seed(1)
  p <- runif(10)
  res <- ADDIS_exhaustive(p, alpha = 0.05)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("pval", "alphai", "R") %in% names(res)))
  expect_equal(nrow(res), length(p))
  expect_true(all(res$R %in% c(0, 1)))
})

test_that("ADDIS_exhaustive returns expected structure for data.frame", {
  set.seed(2)
  p <- runif(8)
  d <- data.frame(id = letters[1:8], pval = p)
  res <- ADDIS_exhaustive(d, alpha = 0.1)
  expect_s3_class(res, "data.frame")
  expect_true(all(c("pval", "alphai", "R", "id") %in% names(res)))
  expect_equal(nrow(res), nrow(d))
  expect_equal(res$id, d$id)
  expect_true(all(res$R %in% c(0, 1)))
})

test_that("ADDIS_exhaustive validates inputs", {
  expect_error(ADDIS_exhaustive(c(-0.1, 0.2)))
  expect_error(ADDIS_exhaustive(c(0.1, 1.2)))
  expect_error(ADDIS_exhaustive(c(0.1, 0.2), alpha = -0.1))
  expect_error(ADDIS_exhaustive(c(0.1, 0.2), alpha = 1.5))
})
