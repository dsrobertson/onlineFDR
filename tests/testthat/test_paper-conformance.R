## Checks that test levels follow the source papers (fixes in version 2.21.2)

test_that("Default gamma from Javanmard and Montanari (2018), equation 31, sums to 1", {
    C <- 0.07908197
    M <- 1e6
    m <- seq_len(M)
    s <- sum(C * log(pmax(m, 2))/(m * exp(sqrt(log(m)))))
    # tail: integral of log(x)/(x exp(sqrt(log x))) from M to infinity
    u <- sqrt(log(M))
    tail <- C * 2 * exp(-u) * (u^3 + 3*u^2 + 6*u + 6)
    expect_equal(s + tail, 1, tolerance = 1e-6)

    expect_equal(LORD(c(0.5, 0.5))$alphai[1], C * log(2) * 0.005)
    expect_equal(LOND(c(0.5, 0.5))$alphai[1], C * log(2) * 0.05)
})

test_that("LORD 3 follows Javanmard and Montanari (2018), equations 9, 15 and 18", {
    set.seed(3)
    n <- 60
    w0 <- 0.005
    b0 <- 0.045
    g <- 0.07908197 * log(pmax(seq_len(n), 2))/(seq_len(n) * exp(sqrt(log(seq_len(n)))))
    for (s in seq_len(20)) {
        p <- ifelse(runif(n) < 0.3, rbeta(n, 0.1, 10), runif(n))
        a <- numeric(n)
        R <- logical(n)
        W <- w0       # W(i)
        Wlast <- w0   # W at the last discovery time
        last <- 0     # last discovery time
        for (i in seq_len(n)) {
            a[i] <- g[i - last] * Wlast
            R[i] <- p[i] <= a[i]
            W <- W - a[i] + R[i] * b0
            if (R[i]) {
                last <- i
                Wlast <- W
            }
        }
        out <- LORD(p, version = 3, w0 = w0, b0 = b0, gammai = g)
        expect_equal(out$alphai, a)
        expect_equal(as.logical(out$R), R)
    }
})

test_that("LORD 3 and 'dep' need w0 < alpha when b0 is not given", {
    expect_error(LORD(c(0.01, 0.5), version = 3, w0 = 0.2),
                 "w0 must be less than alpha.")
    expect_error(LORD(c(0.01, 0.5), version = "dep", w0 = 0.05),
                 "w0 must be less than alpha.")
})

test_that("LORD discard needs w0 <= tau.discard * alpha (Tian and Ramdas, 2019)", {
    expect_error(LORD(c(0.01, 0.5), version = "discard", tau.discard = 0.05, w0 = 0.005),
                 "w0 must not be greater than tau.discard \\* alpha.")
    expect_error(LORD(c(0.01, 0.5), version = "discard", tau.discard = 0.5, w0 = 0.025), NA)

    # default w0 = min(alpha/10, tau.discard * alpha/2)
    C <- 0.07908197
    expect_equal(LORD(c(0.5, 0.5), version = "discard", tau.discard = 0.05)$alphai[1],
                 C * log(2) * 0.05 * 0.05/2)
    expect_equal(LORD(c(0.5, 0.5), version = "discard")$alphai[1],
                 C * log(2) * 0.005)
})

test_that("LORD dep checks the condition of Javanmard and Montanari (2018), example 3.7", {
    # sum(xi) = 0.9 <= alpha/b0 = 1.11, but sum(xi * (1 + log(i))) = 1.44
    xi <- c(0.3, 0.3, 0.3)
    expect_error(LORD(c(0.01, 0.5, 0.5), version = "dep", gammai = xi),
                 "gammai\\*\\(1 \\+ log\\(i\\)\\) must be <= alpha/b0")
    expect_error(suppressWarnings(LORDdep(c(0.01, 0.5, 0.5), xi = xi)),
                 "xi\\*\\(1 \\+ log\\(i\\)\\) must not be greater than alpha/b0")

    # the default sequence satisfies the condition
    i <- seq_len(1e6)
    expect_lt(sum(0.139307/(i * log(pmax(i, 2))^3) * (1 + log(i))), 1)

    for (N in c(1, 10, 1000)) {
        xi <- setBound("LORDdep", N = N)
        expect_equal(sum(xi * (1 + log(seq_len(N)))), 1)
        expect_error(LORD(runif(N), version = "dep", gammai = xi), NA)
    }
})

test_that("BatchStBH and StoreyBH do not reject above the Storey-BH threshold when alpha_t >= 1", {
    # batch 2 has one p-value: pi0 = (1 + 1)/(1 * 0.5) = 4, threshold alpha_2/4
    d <- data.frame(id = seq_len(41), batch = c(rep(1, 40), 2),
                    pval = c(rep(1e-10, 40), 0.9))
    out <- BatchStBH(d, alpha = 0.05)
    expect_gt(out$alphai[41], 1)
    expect_equal(out$R[41], 0)

    # pi0 = (4 + 1)/(0.5 * 4) = 2.5, so adjusted p-values are at least 2.25
    expect_equal(StoreyBH(rep(0.9, 4), alpha = 1)$R, rep(0, 4))
})

test_that("SAFFRONstar batch caps the first batch at lambda", {
    out <- SAFFRONstar(rep(0.5, 6), alpha = 0.5, version = "batch",
                       batch.sizes = c(3, 3), lambda = 0.01)
    expect_true(all(out$alphai <= 0.01))
})

test_that("ADDIS_spending with dep = TRUE follows Tian and Ramdas (2021), equation 23", {
    ref <- function(p, L, g, alpha, lambda, tau) {
        n <- length(p)
        a <- numeric(n)
        for (i in seq_len(n)) {
            j <- seq_len(max(0, i - L[i] - 1))    # j < i - L_i
            t <- 1 + min(L[i], i - 1) + sum(p[j] <= tau) - sum(p[j] <= lambda)
            a[i] <- alpha * (tau - lambda) * g[t]
        }
        a
    }
    set.seed(23)
    n <- 50
    g <- 0.4374901658/(seq_len(n)^(1.6))
    for (s in seq_len(20)) {
        L <- numeric(n)
        for (i in seq_len(n)[-1]) L[i] <- min(L[i - 1] + 1, sample(0:4, 1))
        p <- ifelse(runif(n) < 0.3, rbeta(n, 0.1, 10), runif(n))
        d <- data.frame(id = seq_len(n), pval = p, lags = L)
        out <- ADDIS_spending(d, dep = TRUE)
        expect_equal(out$alphai, ref(p, L, g, 0.05, 0.25, 0.5))
    }
})

test_that("LORDstar async works for a single p-value", {
    d <- data.frame(id = 1, pval = 0.001, decision.times = 1)
    for (s in seq_len(20)) {
        out <- LORDstar(d, version = "async")
        expect_equal(out$alphai, 0.07908197 * log(2) * 0.005)
    }
})

test_that("A gammai (or betai) shorter than the number of p-values gives an error", {
    p <- c(0.01, 0.2, 0.5, 0.7)
    msg <- "must have length at least equal to the number of p-values."
    g <- c(0.5, 0.3)
    expect_error(ADDIS(p, gammai = g), msg)
    expect_error(ADDIS_spending(p, gammai = g), msg)
    expect_error(Alpha_investing(p, gammai = g), msg)
    expect_error(LORD(p, gammai = g), msg)
    expect_error(LORD(p, version = 3, gammai = g), msg)
    expect_error(LORD(p, version = "dep", gammai = c(0.1, 0.1)), msg)
    expect_error(SAFFRON(p, gammai = g), msg)
    expect_error(online_fallback(p, gammai = g), msg)

    d <- data.frame(id = seq_len(4), pval = p, decision.times = seq_len(4))
    expect_error(LORDstar(d, version = "async", gammai = g), msg)
    expect_error(SAFFRONstar(d, version = "async", gammai = g), msg)
    expect_error(LONDstar(d, version = "async", betai = c(0.01, 0.01)), msg)
})

test_that("setBound output passes the sum checks (rounding error is allowed for)", {
    for (N in c(2, 10, 37, 100, 263)) {
        p <- runif(N)
        d <- data.frame(id = seq_len(N), pval = p, decision.times = seq_len(N))
        expect_error(LOND(p, betai = setBound("LOND", N = N)), NA)
        expect_error(LORD(p, gammai = setBound("LORD", N = N)), NA)
        expect_error(SAFFRON(p, gammai = setBound("SAFFRON", N = N)), NA)
        expect_error(ADDIS(p, gammai = setBound("ADDIS", N = N)), NA)
        expect_error(Alpha_investing(p, gammai = setBound("Alpha_investing", N = N)), NA)
        expect_error(ADDIS_spending(p, gammai = setBound("ADDIS_spending", N = N)), NA)
        expect_error(online_fallback(p, gammai = setBound("online_fallback", N = N)), NA)
        expect_error(LONDstar(d, version = "async", betai = setBound("LONDstar", N = N)), NA)
        expect_error(LORDstar(d, version = "async", gammai = setBound("LORDstar", N = N)), NA)
        expect_error(SAFFRONstar(d, version = "async", gammai = setBound("SAFFRONstar", N = N)), NA)
    }
})

test_that("Batch versions of the *star functions keep genuine zero levels", {
    p <- c(0.001, 0.5, 0.5, 0.001, 0.5, 0.5)
    g <- c(0.6, 0.4, 0, 0, 0, 0)
    out <- LORDstar(p, version = "batch", batch.sizes = c(3, 3), gammai = g)
    expect_equal(nrow(out), 6)
    expect_equal(out$alphai[3], 0)
    out <- SAFFRONstar(p, version = "batch", batch.sizes = c(3, 3), gammai = g)
    expect_equal(nrow(out), 6)
    expect_equal(out$alphai[3], 0)
    out <- LONDstar(p, version = "batch", batch.sizes = c(3, 3), betai = 0.05 * g)
    expect_equal(nrow(out), 6)
    expect_equal(out$alphai[3], 0)
})

test_that("gammai must be non-increasing where the papers require it", {
    p <- c(0.01, 0.2, 0.5, 0.7)
    g <- c(0.1, 0.3, 0.2, 0.1)
    msg <- "The elements of gammai must be non-increasing."
    expect_error(LORD(p, gammai = g), msg)
    expect_error(SAFFRON(p, gammai = g), msg)
    expect_error(ADDIS(p, gammai = g), msg)
    expect_error(Alpha_investing(p, gammai = g), msg)
    d <- data.frame(id = seq_len(4), pval = p, decision.times = seq_len(4))
    expect_error(LORDstar(d, version = "async", gammai = g), msg)
    expect_error(SAFFRONstar(d, version = "async", gammai = g), msg)
    # not required for LORD dep, ADDIS-spending or Alpha-spending
    expect_error(ADDIS_spending(p, gammai = g), NA)
    expect_error(Alpha_spending(p, gammai = g), NA)
})

test_that("lambda = 1 gives an error for SAFFRON, SAFFRONstar, BatchStBH and StoreyBH", {
    p <- c(0.01, 0.2, 0.5, 0.7)
    msg <- "lambda must be between 0 and 1."
    expect_error(SAFFRON(p, lambda = 1), msg)
    d <- data.frame(id = seq_len(4), pval = p, decision.times = seq_len(4))
    expect_error(SAFFRONstar(d, version = "async", lambda = 1), msg)
    expect_error(BatchStBH(data.frame(id = 1:4, batch = c(1, 1, 2, 2), pval = p), lambda = 1), msg)
    expect_error(StoreyBH(p, lambda = 1), msg)
})

test_that("StoreyBH accepts a data frame", {
    p <- c(2.90e-08, 0.06743, 0.01514, 0.08174, 0.00171, 3.60e-05, 0.79149)
    out <- StoreyBH(data.frame(id = seq_along(p), pval = p))
    expect_equal(out$R, StoreyBH(p)$R)
    expect_equal(out$id, seq_along(p))
})

test_that("bonfInfinite and LORDdep sort data frames by date", {
    d <- data.frame(id = c("a", "b", "c"),
                    date = as.Date(c("2020-01-03", "2020-01-01", "2020-01-02")),
                    pval = c(0.5, 1e-6, 0.2))
    out <- suppressWarnings(bonfInfinite(d, random = FALSE))
    expect_equal(out$id, c("b", "c", "a"))
    expect_equal(out$R, c(1, 0, 0))
    out <- suppressWarnings(LORDdep(d, random = FALSE))
    expect_equal(out$id, c("b", "c", "a"))
    expect_equal(out$R, c(1, 0, 0))
})

test_that("ADDIS_spending levels are at most lambda (Tian and Ramdas 2021)", {
    # alpha * (tau - lambda) * gamma_1 = 0.05 * 0.49 * 0.437 = 0.0107 > lambda
    out <- ADDIS_spending(c(0.0105, 0.5, 0.5), lambda = 0.01, tau = 0.5)
    expect_true(all(out$alphai <= 0.01))
    expect_equal(out$R[1], 0)
    d <- data.frame(id = 1:3, pval = c(0.0105, 0.5, 0.5), lags = c(0, 1, 1))
    out <- ADDIS_spending(d, lambda = 0.01, tau = 0.5, dep = TRUE)
    expect_true(all(out$alphai <= 0.01))
})

test_that("supLORD follows Xu and Ramdas (2021), equations 10, 16, 17, 20 and 21", {
    ref <- function(p, delta, eps, r, eta, rho, g) {
        L <- log(1/delta)
        f <- function(a) log(1 + L/a) - L/(a + L) - L/(eps*r)
        a <- exp(uniroot(function(la) f(exp(la)), c(-50, 20), tol = 1e-13)$root)
        logbar <- L/(a*log(1 + L/a))
        beta0 <- (eps*r/logbar - a)/r
        b <- eps/logbar
        gbar <- function(i, W) {
            e <- eta*W/beta0
            if (e <= 1) return(g[i])
            if (i > rho) return(0)
            # g[i]^e/sum(g[seq_len(rho)]^e), scaled to avoid underflow
            lg <- log(g[seq_len(rho)])
            exp(e*(log(g[i]) - max(lg)))/sum(exp(e*(lg - max(lg))))
        }
        n <- length(p)
        alpha <- R <- numeric(n)
        W <- beta0
        t <- Wt <- bt <- numeric(0)
        for (k in seq_len(n)) {
            alpha[k] <- gbar(k, beta0)*beta0
            for (j in seq_along(t)) alpha[k] <- alpha[k] + gbar(k - t[j], Wt[j])*bt[j]
            R[k] <- p[k] <= alpha[k]
            bk <- if (length(t) + 1 <= r - 1) beta0 else b
            W <- W + bk*R[k] - alpha[k]
            if (R[k]) {
                t <- c(t, k)
                Wt <- c(Wt, W)
                bt <- c(bt, bk)
            }
        }
        list(alpha = alpha, R = R)
    }
    set.seed(17)
    n <- 150
    g <- 0.07908197 * log(pmax(seq_len(n), 2))/(seq_len(n) * exp(sqrt(log(seq_len(n)))))
    for (s in seq_len(20)) {
        p <- ifelse(runif(n) < 0.4, rbeta(n, 0.05, 20), runif(n))
        eta <- sample(c(0.05, 1, 3, 20), 1)
        rho <- sample(c(5, 30, n), 1)
        r <- sample(c(1, 5, 30), 1)
        out <- supLORD(p, delta = 0.05, eps = 0.15, r = r, eta = eta, rho = rho, gammai = g)
        rf <- ref(p, 0.05, 0.15, r, eta, rho, g)
        expect_equal(out$alphai, rf$alpha, tolerance = 1e-8)
        expect_equal(out$R, as.numeric(rf$R))
    }
})

test_that("supLORD spends the first rejection's boost straight away", {
    # with eta = 20, every rejection's boost is spent within rho = 2 steps
    p <- rep(1e-10, 20)
    out <- supLORD(p, delta = 0.05, eps = 0.15, r = 30, eta = 20, rho = 2)
    expect_equal(sum(out$R), 20)
})

test_that("supLORD handles large eta, rho > N and short gammai", {
    set.seed(5)
    p <- ifelse(runif(200) < 0.5, rbeta(200, 0.05, 50), runif(200))
    out <- supLORD(p, delta = 0.05, eps = 0.15, r = 5, eta = 500, rho = 10)
    expect_true(all(is.finite(out$alphai)) && all(out$alphai >= 0))
    expect_gt(sum(out$R), 0)
    
    out <- supLORD(runif(10), delta = 0.05, eps = 0.15, r = 3, eta = 2, rho = 50)
    expect_true(all(is.finite(out$alphai)))
    
    expect_error(supLORD(runif(10), eps = 0.15, r = 3, eta = 2, rho = 50, 
                         gammai = rep(0.05, 10)),
                 "rho must not be greater than the length of gammai.")
    expect_error(supLORD(runif(10), eps = 0.15, r = 3, eta = 2, rho = 5, 
                         gammai = rep(0.05, 8)),
                 "gammai must have length at least equal to the number of p-values.")
    
    # small eps * r: the offset a is still the positive root, so levels are positive
    out <- supLORD(runif(20), delta = 0.05, eps = 0.01, r = 1, eta = 0.5, rho = 5)
    expect_true(all(out$alphai > 0))
})
