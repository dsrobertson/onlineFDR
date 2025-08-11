# Diagnostic script for numerical precision issues on macOS (or any platform)
# Add this to your repo and run it in GitHub Actions to see the output

alpha <- 0.04
N <- 20

cat('alpha:', alpha, '\n')
cat('N:', N, '\n')

betai <- rep(alpha/N, N)
sum_betai <- sum(betai)
diff <- sum_betai - alpha
eps <- .Machine$double.eps

cat('sum(rep(alpha/N, N)):', format(sum_betai, digits=20), '\n')
cat('Difference (sum - alpha):', format(diff, digits=20), '\n')
cat('Machine epsilon:', format(eps, digits=20), '\n')
cat('Tolerance (eps * N):', format(eps * N, digits=20), '\n')
cat('Is difference < tolerance?', abs(diff) < eps * N, '\n')

# Also print the actual betai vector for inspection
cat('betai:', paste(format(betai, digits=20), collapse=', '), '\n')
