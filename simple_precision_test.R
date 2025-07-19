# Simple precision test
alpha <- 0.04
N <- 20

# Manual calculation of what setBound does for LOND
bound <- rep(alpha/N, N)

cat("Sum of bound:", format(sum(bound), digits=20), "\n")
cat("Alpha:       ", format(alpha, digits=20), "\n")
cat("Difference:  ", format(sum(bound) - alpha, digits=20), "\n")

# Test tolerances
eps <- .Machine$double.eps
tol1 <- eps * N
tol100 <- 100 * eps * N

cat("Machine epsilon:", format(eps, digits=20), "\n")
cat("Basic tolerance (eps * N):", format(tol1, digits=20), "\n")
cat("100x tolerance:", format(tol100, digits=20), "\n")

diff <- abs(sum(bound) - alpha)
cat("Absolute difference:", format(diff, digits=20), "\n")
cat("Passes basic tolerance:", diff <= tol1, "\n")
cat("Passes 100x tolerance:", diff <= tol100, "\n")
