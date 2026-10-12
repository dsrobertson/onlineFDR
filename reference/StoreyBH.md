# StoreyBH: Offline FDR control using the St-BH procedure

Implements the Storey-BH algorithm for offline FDR control, as presented
by Storey (2002), with the null proportion estimate of Storey, Taylor
and Siegmund (2004).

## Usage

``` r
StoreyBH(d, alpha = 0.05, lambda = 0.5)
```

## Arguments

- d:

  Either a vector of p-values, or a dataframe with the column: p-value
  (\`pval').

- alpha:

  Overall significance level of the FDR procedure, the default is 0.05.

- lambda:

  Threshold for Storey-BH, must be in (0, 1). Defaults to 0.5.

## Value

- out:

  A dataframe with the original data `d` and the indicator function of
  discoveries `R`. Hypothesis \\i\\ is rejected if the \\i\\-th p-value
  is less than or equal to \\(k/n)\alpha / \hat{\pi}\_0\\, where \\k\\
  is the largest \\r\\ such that the \\r\\-th smallest p-value is at
  most \\(r/n)\alpha / \hat{\pi}\_0\\, \\n\\ is the total number of
  hypotheses and \\\hat{\pi}\_0 = (1 + n\_\lambda) / (n(1 - \lambda))\\
  is Storey's estimate of the proportion of true null hypotheses, with
  \\n\_\lambda\\ the number of p-values greater than \\\lambda\\. If
  hypothesis \\i\\ is rejected, `R[i] = 1` (otherwise `R[i] = 0`).

## Details

The function takes as its input either a vector of p-values, or a
dataframe with a column of p-values (\`pval').

## References

Storey, J.D. (2002). A direct approach to false discovery rates.
*Journal of the Royal Statistical Society (Series B)*, 64(3):479-498.

Storey, J.D., Taylor, J.E. and Siegmund, D. (2004). Strong control,
conservative point estimation and simultaneous conservative consistency
of false discovery rates: a unified approach. *Journal of the Royal
Statistical Society (Series B)*, 66(1):187-205.

## Examples

``` r

pvals <- c(2.90e-08, 0.06743, 0.01514, 0.08174, 0.00171,
        3.60e-05, 0.79149, 0.27201, 0.28295, 7.59e-08,
        0.69274, 0.30443, 0.00136, 0.72342, 0.54757)

StoreyBH(pvals)
#>          pval R
#> 1  2.9000e-08 1
#> 2  6.7430e-02 0
#> 3  1.5140e-02 1
#> 4  8.1740e-02 0
#> 5  1.7100e-03 1
#> 6  3.6000e-05 1
#> 7  7.9149e-01 0
#> 8  2.7201e-01 0
#> 9  2.8295e-01 0
#> 10 7.5900e-08 1
#> 11 6.9274e-01 0
#> 12 3.0443e-01 0
#> 13 1.3600e-03 1
#> 14 7.2342e-01 0
#> 15 5.4757e-01 0
```
