# StoreyBH: Offline FDR control using the St-BH procedure

Implements the Storey-BH algorithm for offline FDR control, as presented
by Storey (2002).

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

  Threshold for Storey-BH, must be between 0 and 1. Defaults to 0.5.

## Value

- ordered_d:

  A dataframe with the original data `d` and the indicator function of
  discoveries `R`. Hypothesis \\i\\ is rejected if the \\i\\-th p-value
  is less than or equal to \\(r/n)\alpha\\, where \\r\\ is the rank of
  the \\i\\-th p-value within an ordered set and \\n\\ is the total
  number of hypotheses. If hypothesis \\i\\ is rejected, `R[i] = 1`
  (otherwise `R[i] = 0`).

## Details

The function takes as its input either a vector of p-values, or a
dataframe with a column of p-values (\`pval').

## References

Storey, J.D. (2002). A direct approach to false discovery rates. *J. R.
Statist. Soc. B*: 64, Part 3, 479-498.

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
