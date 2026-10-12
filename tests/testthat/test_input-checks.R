# Input checks and fixes from the second code review (version 2.21.2)

p6 <- c(1e-6, 0.2, 1e-5, 0.6, 1e-4, 0.03)

test_that("empty or all-missing input gives an error rather than a crash", {
    for (f in list(LOND, LORD, SAFFRON, ADDIS, Alpha_investing, ADDIS_spending,
                   online_fallback, Alpha_spending, StoreyBH)) {
        expect_error(f(numeric(0)), "d must contain at least one p-value.")
    }
    expect_error(suppressWarnings(LORD(c(NA_real_, NA_real_))),
                 "d must contain at least one p-value.")
    expect_error(LORDstar(data.frame(pval = numeric(0), decision.times = numeric(0)),
                          version = "async"),
                 "d must contain at least one p-value.")
})

test_that("missing p-values remove only their own rows", {
    d <- data.frame(id = c("a", NA, "c", "d"), pval = c(1e-6, 0.2, NA, 0.6))
    expect_warning(out <- LORD(d), "Missing p-values were omitted.")
    expect_equal(nrow(out), 3)
    expect_equal(out$pval, c(1e-6, 0.2, 0.6))

    d_async <- data.frame(pval = c(1e-6, NA, 0.6), decision.times = 2:4)
    expect_error(LORDstar(d_async, version = "async"),
                 "Missing p-values cannot be omitted")
})

test_that("a data frame without a pval column gives the specific error", {
    expect_error(LORD(data.frame(p = p6)),
                 "The dataframe d is missing a column 'pval' of p-values.")
})

test_that("LOND, Alpha_investing and online_fallback return the identifiers", {
    d <- data.frame(id = paste0("t", 1:6),
                    date = as.Date(rep(c("2020-01-01", "2020-01-02"), each = 3)),
                    pval = p6)
    for (f in list(LOND, Alpha_investing, online_fallback)) {
        set.seed(3)
        out <- f(d)
        expect_setequal(out$id, d$id)
        expect_equal(out$pval, d$pval[match(out$id, d$id)])
    }
})

test_that("decision times must be finite, and very large ones mean 'not yet known'", {
    p <- rep(1e-6, 5)
    expect_error(LORDstar(data.frame(pval = p, decision.times = c(Inf, 3:6)),
                          version = "async"),
                 "All decision times should be finite.")
    late <- data.frame(pval = p, decision.times = c(1000, 3:6))
    huge <- data.frame(pval = p, decision.times = c(3e9, 3:6))
    expect_equal(LORDstar(huge, version = "async")$alphai,
                 LORDstar(late, version = "async")$alphai)
    expect_equal(SAFFRONstar(huge, version = "async")$alphai,
                 SAFFRONstar(late, version = "async")$alphai)
    expect_equal(LONDstar(huge, version = "async")$alphai,
                 LONDstar(late, version = "async")$alphai)
    expect_equal(ADDIS(huge, async = TRUE)$alphai,
                 ADDIS(late, async = TRUE)$alphai)
})

test_that("ADDIS with async = TRUE checks the decision times", {
    expect_error(ADDIS(data.frame(pval = p6, decision.times = c(1.5, 3:7)), async = TRUE),
                 "All decision times should be integers.")
    expect_error(ADDIS(data.frame(pval = p6, decision.times = c(-3, 3:7)), async = TRUE),
                 "All decision times should be positive integers.")
})

test_that("ADDIS_spending with dep = TRUE needs a column of lags", {
    expect_error(ADDIS_spending(p6, dep = TRUE),
                 "d needs to be a dataframe with a column of lags")
    expect_error(ADDIS_spending(data.frame(id = 1:6, pval = p6), dep = TRUE),
                 "d needs to be a dataframe with a column of lags")
})

test_that("a too-short gammai gives an error rather than NA levels", {
    expect_error(Alpha_spending(p6, gammai = rep(0.1, 3)),
                 "gammai must have length at least equal to the number of p-values.")
    db <- data.frame(id = 1:6, batch = rep(1:3, each = 2), pval = p6)
    for (f in list(BatchBH, BatchPRDS, BatchStBH)) {
        expect_error(f(db, gammai = c(0.5, 0.3)),
                     "gammai must have length at least equal to the number of batches.")
    }
})

test_that("LORD dep default with w0 > b0 does not depend on N and uses the full sum", {
    set.seed(1)
    p <- pnorm(-rnorm(1000, mean = 3 * rbinom(1000, 1, 0.3)))
    a_short <- LORD(p[1:6], version = "dep", w0 = 0.04, b0 = 0.01)$alphai
    a_long <- LORD(p, version = "dep", w0 = 0.04, b0 = 0.01)$alphai
    expect_equal(a_short, a_long[1:6])
    gamma1 <- 0.05/(5.068667247 * 0.04 + 2.109742802 * 0.01)/log(2)^3
    expect_equal(a_long[1], gamma1 * 0.04)
})

test_that("the LORD dep constants bound the infinite sums from above", {
    M <- 1e6
    i <- seq_len(M)
    lm <- log(pmax(i, 2))
    A <- sum(1/(i * lm^3)) + 1/(2 * log(M)^2) - 1/(2 * M * log(M)^3)
    B <- sum(log(i)/(i * lm^3)) + 1/log(M) - 1/(2 * M * log(M)^2)
    expect_true(A <= 5.068667247 && A > 5.068667247 - 1e-8)
    expect_true(B <= 2.109742802 && B > 2.109742802 - 1e-8)
    # sum of (1 + log i)/(i log(max(i,2))^3) times 0.1393066 is at most 1
    expect_lte(0.1393066 * (A + B), 1)
    expect_gt(0.1393066 * (A + B), 1 - 1e-6)
})

test_that("LOND with dep = TRUE accepts a betai longer than N", {
    expect_no_warning(LOND(p6, betai = rep(0.001, 7), dep = TRUE))
})

test_that("corrected error messages", {
    expect_error(LORD(p6, version = 1), "superseded by LORD\\+\\+")
    expect_error(SAFFRONstar(p6, version = "batch", batch.sizes = c(3, 3), w0 = 1),
                 "w0 must not be greater than alpha.")
    expect_error(setBound("x", N = 10), "Alpha_investing, Alpha_spending, ADDIS_spending")
})
