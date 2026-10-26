test.pval <- c(1e-07, 3e-04, 0.1, 5e-04)
test.decision.times <- seq_len(4)
test.decision.times2 <- rep(4,4)
test.lags <- rep(0,4)
test.lags2 <- rep(4,4)
test.batch.sizes <- rep(1,4)
test.batch.sizes <- 4
test.df <- data.frame(id = seq_len(4), pval = test.pval, decision.times = test.decision.times)
test.df2 <- data.frame(id = seq_len(4), pval = test.pval, decision.times = test.decision.times2)
test.df3 <-  data.frame(id = seq_len(4), pval = test.pval, lags = test.lags)
test.df4 <-  data.frame(id = seq_len(4), pval = test.pval, lags = test.lags2)


test_that("Errors for edge cases", {
    
    expect_error(SAFFRONstar(test.df, version='async', alpha = -0.1),
                 "alpha must be between 0 and 1.")
    
    expect_error(SAFFRONstar(test.df, version='async', lambda=2),
                 "lambda must be between 0 and 1.")
    
    expect_error(SAFFRONstar(test.df, version='async', gammai = -1),
                 "All elements of gammai must be non-negative.")
    
    expect_error(SAFFRONstar(test.df, version='async', gammai=2),
                 "The sum of the elements of gammai must not be greater than 1.")
    
    expect_error(SAFFRONstar(test.df, version='async', w0 = -0.01),
                 "w0 must be non-negative.")
    
    # expect_error(SAFFRONstar(0.1, version='async', decision.times=1, w0 = 1),
    #              "w0 must be less than (1-lambda)*alpha")

})

test_that("Correct rejections for version async", {
    
    expect_identical(SAFFRONstar(test.df, version='async')$R,
                     c(1,1,0,1))
    
    expect_identical(SAFFRONstar(test.df2, version='async')$R,
                     c(1,1,0,1))
})

test_that("Correct rejections for version dep", {
    
    expect_identical(SAFFRONstar(test.df3, version='dep')$R,
                     c(1,1,0,1))
    
    expect_identical(SAFFRONstar(test.df4, version='dep')$R,
                     c(1,1,0,1))
})


test_that("Correct rejections for version batch", {
    
    expect_identical(SAFFRONstar(0.1, version='batch', batch.sizes=1)$R, 0)
    
    expect_identical(SAFFRONstar(c(0.1,0.1), version='batch',
                                 batch.sizes=c(1,1))$R, c(0,0))
    
    expect_identical(SAFFRONstar(test.pval, version='batch',
                              batch.sizes = rep(1,4))$R,
                     c(1,1,0,1))
    
    expect_identical(SAFFRONstar(test.pval, version='batch',
                              batch.sizes = 4)$R,
                     c(1,1,0,1))
})


test_that("Check that SAFFRON is a special case of the SAFFRONstar
          algorithms", {
              expect_equal(SAFFRON(test.df)$alphai,
                           SAFFRONstar(test.df, version='async')$alphai)
              
              expect_equal(SAFFRON(test.df)$alphai,
                           SAFFRONstar(test.df3, version='dep')$alphai)
              
              expect_equal(SAFFRON(test.df)$alphai,
                           SAFFRONstar(test.pval, version='batch',
                                    batch.sizes = rep(1,4))$alphai)              
})

test_that("Async SAFFRONstar follows Zrnic et al. (2021), Algorithm 3", {
  # Reference implementation of Algorithm 3 with asynchronous conflict sets:
  # test t uses the outcomes of tests i < t with E[i] < t, r_k is the k-th
  # smallest decision time among those rejections, and C_k^+ counts the
  # candidates among tests that started after r_k
  ref <- function(p, E, alpha, g, w0, lambda) {
    a <- numeric(length(p))
    for (t in seq_along(p)) {
      prev <- seq_len(t - 1)
      known <- prev[E[prev] < t]
      cand <- known[p[known] <= lambda]
      r <- sort(E[known][p[known] <= a[known]])
      val <- w0 * g[t - length(cand)]
      for (k in seq_along(r)) {
        val <- val + (if (k == 1) alpha - w0 else alpha) *
          g[t - r[k] - sum(cand > r[k])]
      }
      a[t] <- min(lambda, (1 - lambda) * val)
    }
    a
  }
  for (maxdelay in c(10, 50)) {
    set.seed(maxdelay)
    n <- 80
    p <- pnorm(-(rnorm(n) + 3 * (runif(n) < 0.3)))
    E <- seq_len(n) + sample(0:maxdelay, n, replace = TRUE)
    g <- 0.4374901658 / seq_len(n + 1)^1.6
    d <- data.frame(id = seq_len(n), pval = p, decision.times = E)
    out <- SAFFRONstar(d, alpha = 0.1, version = "async", gammai = g,
                       w0 = 0.04, lambda = 0.4)
    expect_equal(out$alphai, ref(p, E, 0.1, g, 0.04, 0.4))
    expect_gt(sum(out$R), 0)
  }
})

test_that("SAFFRONstar inputs are correct when given vector input", {
    expect_error(SAFFRONstar(test.pval, version="async"),
                 "d needs to be a dataframe with a column of decision.times")
    expect_error(SAFFRONstar(test.pval, version="dep"),
                 "d needs to be a dataframe with a column of lags")
})
