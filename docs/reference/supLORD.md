# supLORD: Online control of the false discovery exceedance (FDX) and the FDR at stopping times

Implements the supLORD procedure, which controls both FDX and FDR,
including the FDR at stopping times, as presented by Xu and Ramdas
(2021).

## Usage

``` r
supLORD(
  d,
  delta = 0.05,
  eps,
  r,
  eta,
  rho,
  gammai,
  random = TRUE,
  display_progress = FALSE,
  date.format = "%Y-%m-%d"
)
```

## Arguments

- d:

  Either a vector of p-values, or a dataframe with three columns: an
  identifier (\`id'), date (\`date') and p-value (\`pval'). If no column
  of dates is provided, then the p-values are treated as being ordered
  in sequence, arriving one at a time.

- delta:

  The probability at which the FDP exceeds eps (at any time step after
  making r rejections). Must be between 0 and 1, defaults to 0.05.

- eps:

  The upper bound on the FDP. Must be between 0 and 1.

- r:

  The threshold of rejections after which the error control begins to
  apply. Must be a positive integer.

- eta:

  Controls the pace at which wealth is spent as a function of the
  algorithm's current wealth. Must be a positive real number.

- rho:

  Controls the length of time before the spending sequence exhausts the
  wealth earned from a rejection. Must be a positive integer.

- gammai:

  Optional vector of \\\gamma_i\\. A default is provided as proposed by
  Javanmard and Montanari (2018).

- random:

  Logical. If `TRUE` (the default), then the order of the p-values in
  each batch (i.e. those that have exactly the same date) is randomised.

- display_progress:

  Logical. If `TRUE` prints out a progress bar for the algorithm
  runtime.

- date.format:

  Optional string giving the format that is used for dates.

## Value

- d.out:

  A dataframe with the original data `d` (which will be reordered if
  there are batches and `random = TRUE`), the supLORD-adjusted
  significance thresholds \\\alpha_i\\ and the indicator function of
  discoveries `R`. Hypothesis \\i\\ is rejected if the \\i\\-th p-value
  is less than or equal to \\\alpha_i\\, in which case `R[i] = 1`
  (otherwise `R[i] = 0`).

## Details

The function takes as its input either a vector of p-values or a
dataframe with three columns: an identifier (\`id'), date (\`date') and
p-value (\`pval'). The case where p-values arrive in batches corresponds
to multiple instances of the same date. If no column of dates is
provided, then the p-values are treated as being ordered in sequence,
arriving one at a time..

The supLORD procedure provably controls the FDX for p-values that are
conditionally superuniform under the null. supLORD also controls the
supFDR and hence the FDR (even at stopping times). Given an overall
significance level \\\alpha\\, we choose a sequence of non-negative
non-increasing numbers \\\gamma_i\\ that sum to 1.

supLORD requires the user to specify r, a threshold of rejections after
which the error control begins to apply, eps, the upper bound on the
false discovery proportion (FDP), and delta, the probability at which
the FDP exceeds eps at any time step after making r rejections. As well,
the user should specify the variables eta, which controls the pace at
which wealth is spent (as a function of the algorithm's current wealth),
and rho, which controls the length of time before the spending sequence
exhausts the wealth earned from a rejection.

Further details of the supLORD procedure can be found in Xu and Ramdas
(2021).

## References

Xu, Z. and Ramdas, A. (2021). Dynamic Algorithms for Online Multiple
Testing. *Annual Conference on Mathematical and Scientific Machine
Learning*, PMLR, 145:955-986.

## Examples

``` r
set.seed(1)
N <- 1000
B <- rbinom(N, 1, 0.5)
Z <- rnorm(N, mean = 3*B)
pval <- pnorm(-Z)

out <- supLORD(pval, eps=0.15, r=30, eta=0.05, rho=30, random=FALSE)
head(out)
#>              d       alphai R
#> 1 4.691912e-01 0.0019301522 0
#> 2 6.167166e-01 0.0004197471 0
#> 3 3.462711e-02 0.0003575072 0
#> 4 1.300690e-03 0.0002973164 0
#> 5 1.606961e-01 0.0002520627 0
#> 6 2.174486e-06 0.0002180533 1
sum(out$R)
#> [1] 390
```
