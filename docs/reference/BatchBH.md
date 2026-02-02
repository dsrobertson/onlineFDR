# BatchBH: Online batch FDR control using the BH procedure

Implements the BatchBH algorithm for online FDR control, as presented by
Zrnic et al. (2020).

## Usage

``` r
BatchBH(d, alpha = 0.05, gammai, display_progress = FALSE)
```

## Arguments

- d:

  A dataframe with three columns: identifiers (\`id'), batch numbers
  (\`batch') and p-values (\`pval').

- alpha:

  Overall significance level of the FDR procedure, the default is 0.05.

- gammai:

  Optional vector of \\\gamma_i\\. A default is provided with
  \\\gamma_j\\ proportional to \\1/j^(1.6)\\.

- display_progress:

  Logical. If `TRUE` prints out a progress bar for the algorithm
  runtime.

## Value

- out:

  A dataframe with the original data `d` and the indicator function of
  discoveries `R`. Hypothesis \\i\\ is rejected if the \\i\\-th p-value
  within the \\t\\-th batch is less than or equal to \\(r/n)\alpha_t\\,
  where \\r\\ is the rank of the \\i\\-th p-value within an ordered set
  and \\n\\ is the total number of hypotheses within the \\t\\-th batch.
  If hypothesis \\i\\ is rejected, `R[i] = 1` (otherwise `R[i] = 0`).

## Details

The function takes as its input a dataframe with three columns:
identifiers (\`id'), batch numbers (\`batch') and p-values (\`pval').

The BatchBH algorithm controls the FDR when the p-values in a batch are
independent, and independent across batches. Given an overall
significance level \\\alpha\\, we choose a sequence of non-negative
numbers \\\gamma_i\\ such that they sum to 1. The algorithm runs the
Benjamini-Hochberg procedure on each batch, where the values of the
adjusted significance thresholds \\\alpha\_{t+1}\\ depend on the number
of previous discoveries.

Further details of the BatchBH algorithm can be found in Zrnic et al.
(2020).

## References

Zrnic, T., Jiang D., Ramdas A. and Jordan M. (2020). The Power of
Batching in Multiple Hypothesis Testing. *International Conference on
Artificial Intelligence and Statistics*, 3806-3815.

## Examples

``` r
sample.df <- data.frame(
id = c('A15432', 'B90969', 'C18705', 'B49731', 'E99902',
    'C38292', 'A30619', 'D46627', 'E29198', 'A41418',
    'D51456', 'C88669', 'E03673', 'A63155', 'B66033'),
pval = c(2.90e-08, 0.06743, 0.01514, 0.08174, 0.00171,
        3.60e-05, 0.79149, 0.27201, 0.28295, 7.59e-08,
        0.69274, 0.30443, 0.00136, 0.72342, 0.54757),
batch = c(rep(1,5), rep(2,6), rep(3,4)))

BatchBH(sample.df)
#>        id       pval batch R      alphai
#> 1  A15432 2.9000e-08     1 1 0.021874508
#> 2  B90969 6.7430e-02     1 0 0.021874508
#> 3  C18705 1.5140e-02     1 0 0.021874508
#> 4  B49731 8.1740e-02     1 0 0.021874508
#> 5  E99902 1.7100e-03     1 1 0.021874508
#> 6  C38292 3.6000e-05     2 1 0.009621196
#> 7  A30619 7.9149e-01     2 0 0.009621196
#> 8  D46627 2.7201e-01     2 0 0.009621196
#> 9  E29198 2.8295e-01     2 0 0.009621196
#> 10 A41418 7.5900e-08     2 1 0.009621196
#> 11 D51456 6.9274e-01     2 0 0.009621196
#> 12 C88669 3.0443e-01     3 0 0.025012888
#> 13 E03673 1.3600e-03     3 1 0.025012888
#> 14 A63155 7.2342e-01     3 0 0.025012888
#> 15 B66033 5.4757e-01     3 0 0.025012888
```
