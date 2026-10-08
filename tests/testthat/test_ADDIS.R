test.pval <- c(1e-07, 3e-04, 0.1, 6e-04)
test.df <- data.frame(id = seq_len(4), pval = test.pval)
test.df2 <- data.frame(id = seq_len(4), pval = test.pval, decision.times = c(2,3,4,NA))
test.df3 <- data.frame(id = seq_len(4), pval = test.pval, decision.times = seq_len(4)+3)
test.df4 <- data.frame(id = seq_len(4), pval = test.pval, decision.times = seq_len(4))

test_that("Errors for edge cases", {

  expect_error(ADDIS(test.df, alpha = -0.1),
               "alpha must be between 0 and 1.")
  
  expect_error(ADDIS(test.df, lambda = -0.1),
               "lambda must be between 0 and tau.")
  
  expect_error(ADDIS(test.df, tau = 0.5, lambda = 0.5),
               "lambda must be between 0 and tau.")

  expect_error(ADDIS(test.df, tau = -0.1),
               "tau must be between 0 and 1.")
  
  expect_error(ADDIS(test.df, gammai = -1),
               "All elements of gammai must be non-negative.")
  
  expect_error(ADDIS(test.df, gammai=2),
               "The sum of the elements of gammai must not be greater than 1.")
  
  expect_error(ADDIS(test.df, w0 = -0.01),
               "w0 must be non-negative.")
  
  expect_error(ADDIS(0.1, w0 = 1),
                "w0 must be less than alpha.")
  
  expect_error(ADDIS(test.df2, async=TRUE),
              "Please provide a decision time for each p-value.")
})


test_that("Correct rejections", {
  
    expect_identical(ADDIS(test.df)$R, c(1,1,0,1))
    
    expect_identical(ADDIS(test.df3,
                           async=TRUE)$R,
                           c(1,1,0,0))
})

test_that("Check that ADDIS with async=FALSE is a special case of async=TRUE", {
              expect_equal(ADDIS(test.pval, async=FALSE)$alphai,
                           ADDIS(test.df4, async=TRUE)$alphai)
})

test_that("ADDIS inputs are correct with async=TRUE", {
  expect_error(ADDIS(test.df, async=TRUE),
               "d needs to have a column of decision.times")
  expect_error(ADDIS(c(0.1, 0.1), async=TRUE),
               "d needs to be a dataframe with a column of decision.times")
})

test_that("Async ADDIS treats p-values up to lambda as candidates, as sync ADDIS does", {
  # p-values in (tau*lambda, lambda] = (0.125, 0.25] separate the old and new
  # candidate thresholds in addis_async_faster
  p <- c(1e-07, 0.2, 0.15, 6e-04, 0.01)
  d <- data.frame(id = seq_along(p), pval = p, decision.times = seq_along(p))
  expect_equal(ADDIS(p, async = FALSE)$alphai,
               ADDIS(d, async = TRUE)$alphai)
})

test_that("Async ADDIS follows Tian and Ramdas (2019), Algorithm 3", {
  # Reference implementation written from the paper: kappa_j, kappa_j^* and C_j^+
  # are defined by decision times, and S^t counts tests still running as selected
  ref <- function(p, E, alpha, g, w0, lambda, tau) {
    a <- numeric(length(p))
    for (t in seq_along(p)) {
      prev <- seq_len(t - 1)
      fin <- prev[E[prev] < t]
      S <- sum(p[fin] <= tau) + sum(E[prev] >= t)
      kap <- sort(E[fin][p[fin] <= a[fin]])
      val <- w0 * g[S - sum(p[fin] <= lambda) + 1]
      for (j in seq_along(kap)) {
        ks <- sum(p[prev] <= tau & E[prev] <= kap[j])
        Cj <- sum(p[prev] <= lambda & E[prev] > kap[j] & E[prev] < t)
        val <- val + (if (j == 1) alpha - w0 else alpha) * g[S - ks - Cj + 1]
      }
      a[t] <- min(lambda, (tau - lambda) * val)
    }
    a
  }
  set.seed(1)
  n <- 60
  p <- pnorm(-(rnorm(n) + 3 * (runif(n) < 0.3)))
  E <- seq_len(n) + sample(0:5, n, replace = TRUE)
  g <- 0.4374901658 / seq_len(n + 1)^1.6
  d <- data.frame(id = seq_len(n), pval = p, decision.times = E)
  out <- ADDIS(d, alpha = 0.1, gammai = g, w0 = 0.05, lambda = 0.25, tau = 0.5, async = TRUE)
  expect_equal(out$alphai, ref(p, E, 0.1, g, 0.05, 0.25, 0.5))
  expect_gt(sum(out$R), 0)
})
