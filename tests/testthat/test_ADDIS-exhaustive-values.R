test_that("ADDIS_exhaustive levels match the E-ADDIS-Spending formula", {
  # Fischer et al. (2024): alpha_i = alpha * gamma_t * (tau - lambda) / (1 - alpha^(i)),
  # where t and alpha^(i) only change when lambda < p_i <= tau
  p <- c(0.3, 0.001, 0.4, 0.6)
  alpha <- 0.05; tau <- 0.5; lambda <- 0.25
  g <- 0.4374901658/(seq_len(5)^1.6)
  a1 <- alpha                       # alpha^(1)
  a2 <- alpha - alpha * g[1]        # p1 = 0.3 is spent
  a4 <- a2 - alpha * g[2]           # p3 = 0.4 is spent
  expected <- c(alpha * g[1] * (tau - lambda) / (1 - a1),
                alpha * g[2] * (tau - lambda) / (1 - a2),
                alpha * g[2] * (tau - lambda) / (1 - a2),
                alpha * g[3] * (tau - lambda) / (1 - a4))
  res <- ADDIS_exhaustive(p, alpha = alpha, tau = tau, lambda = lambda)
  expect_equal(res$alphai, expected)
  expect_identical(res$R, c(0, 1, 0, 0))
})

test_that("ADDIS_exhaustive rejects invalid tuning parameters", {
  p <- c(0.1, 0.2, 0.3)
  expect_error(ADDIS_exhaustive(p, alpha = 1),
               "alpha must be between 0 and 1.")
  expect_error(ADDIS_exhaustive(p, tau = 0.5, lambda = 0.5),
               "lambda must be between 0 and tau.")
  expect_error(ADDIS_exhaustive(p, alpha = 0.5, tau = 1, lambda = 0.01),
               "lambda must be at least alpha * tau.", fixed = TRUE)
  expect_error(ADDIS_exhaustive(p, gamma = rep(0.5, 4)),
               "The sum of the elements of gamma must not be greater than 1.")
})

test_that("ADDIS_exhaustive never rejects a p-value above lambda", {
  # At the boundary lambda = alpha * tau with all weight on gamma[1], the
  # first level equals lambda exactly
  set.seed(3)
  p <- runif(200)
  res <- ADDIS_exhaustive(p, alpha = 0.2, tau = 0.8, lambda = 0.16,
                          gamma = c(1, rep(0, 200)))
  expect_equal(res$alphai[1], 0.16)
  expect_true(all(res$alphai <= 0.16 + 1e-12))
  expect_true(all(res$R[p > 0.16] == 0))
})
