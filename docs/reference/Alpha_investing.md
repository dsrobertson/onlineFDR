# Alpha-investing for online FDR control

Implements a variant of the Alpha-investing algorithm of Foster and
Stine (2008) that guarantees FDR control, as proposed by Ramdas et al.
(2018). This procedure uses SAFFRON's update rule with the constant
\\\lambda\\ replaced by a sequence \\\lambda_i = \alpha_i\\. This is
also equivalent to using the ADDIS algorithm with \\\tau = 1\\ and
\\\lambda_i = \alpha_i\\.

## Usage

``` r
Alpha_investing(
  d,
  alpha = 0.05,
  gammai,
  w0,
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
  in sequence.

- alpha:

  Overall significance level of the FDR procedure, the default is 0.05.

- gammai:

  Optional vector of \\\gamma_i\\. A default is provided with
  \\\gamma_j\\ proportional to \\1/j^(1.6)\\.

- w0:

  Initial \`wealth' of the procedure, defaults to \\\alpha/2\\. Must be
  between 0 and \\\alpha\\.

- random:

  Logical. If `TRUE` (the default), then the order of the p-values in
  each batch (i.e. those that have exactly the same date) is randomised.

- display_progress:

  Logical. If `TRUE` prints out a progress bar for the algorithm
  runtime.

- date.format:

  Optional string giving the format that is used for dates.

## Value

- out:

  A dataframe with the original data `d` (which will be reordered if
  there are batches and `random = TRUE`), the LORD-adjusted significance
  thresholds \\\alpha_i\\ and the indicator function of discoveries `R`.
  Hypothesis \\i\\ is rejected if the \\i\\-th p-value is less than or
  equal to \\\alpha_i\\, in which case `R[i] = 1` (otherwise
  `R[i] = 0`).

## Details

The function takes as its input either a vector of p-values or a
dataframe with three columns: an identifier (\`id'), date (\`date') and
p-value (\`pval'). The case where p-values arrive in batches corresponds
to multiple instances of the same date. If no column of dates is
provided, then the p-values are treated as being ordered in sequence.

The Alpha-investing procedure provably controls FDR for independent
p-values. Given an overall significance level \\\alpha\\, we choose a
sequence of non-negative non-increasing numbers \\\gamma_i\\ that sum
to 1. Alpha-investing depends on a constant \\w_0\\, which satisfies \\0
\le w_0 \le \alpha\\ and represents the initial \`wealth' of the
procedure.

Further details of the Alpha-investing procedure and its modification
can be found in Foster and Stine (2008) and Ramdas et al. (2018).

## References

Foster, D. and Stine R. (2008). \\\alpha\\-investing: a procedure for
sequential control of expected false discoveries. *Journal of the Royal
Statistical Society (Series B)*, 29(4):429-444.

Ramdas, A., Zrnic, T., Wainwright M.J. and Jordan, M.I. (2018). SAFFRON:
an adaptive algorithm for online control of the false discovery rate.
*Proceedings of the 35th International Conference in Machine Learning*,
80:4286-4294.

## See also

[`SAFFRON`](https://dsrobertson.github.io/onlineFDR/reference/SAFFRON.md)
uses the update rule of Alpha-investing but with constant \\\lambda\\.

## Examples

``` r
sample.df <- data.frame(
id = c('A15432', 'B90969', 'C18705', 'B49731', 'E99902',
    'C38292', 'A30619', 'D46627', 'E29198', 'A41418',
    'D51456', 'C88669', 'E03673', 'A63155', 'B66033'),
date = as.Date(c(rep('2014-12-01',3),
               rep('2015-09-21',5),
                rep('2016-05-19',2),
                '2016-11-12',
               rep('2017-03-27',4))),
pval = c(2.90e-08, 0.06743, 0.01514, 0.08174, 0.00171,
        3.60e-05, 0.79149, 0.27201, 0.28295, 7.59e-08,
        0.69274, 0.30443, 0.00136, 0.72342, 0.54757))

Alpha_investing(sample.df, random=FALSE)
#>          pval      alphai R
#> 1  2.9000e-08 0.010818925 1
#> 2  6.7430e-02 0.021406257 0
#> 3  1.5140e-02 0.007164201 0
#> 4  8.1740e-02 0.003757589 0
#> 5  1.7100e-03 0.002374706 1
#> 6  3.6000e-05 0.023680499 1
#> 7  7.9149e-01 0.044095287 0
#> 8  2.7201e-01 0.015842430 0
#> 9  2.8295e-01 0.008711190 0
#> 10 7.5900e-08 0.005700295 1
#> 11 6.9274e-01 0.026865786 0
#> 12 3.0443e-01 0.011205456 0
#> 13 1.3600e-03 0.006863124 1
#> 14 7.2342e-01 0.027979664 0
#> 15 5.4757e-01 0.011945806 0

set.seed(1); Alpha_investing(sample.df)
#>          pval      alphai R
#> 1  2.9000e-08 0.010818925 1
#> 2  6.7430e-02 0.021406257 0
#> 3  1.5140e-02 0.007164201 0
#> 4  8.1740e-02 0.003757589 0
#> 5  1.7100e-03 0.002374706 1
#> 6  2.7201e-01 0.023680499 0
#> 7  3.6000e-05 0.008803369 1
#> 8  7.9149e-01 0.029838354 0
#> 9  7.5900e-08 0.012084065 1
#> 10 2.8295e-01 0.032981505 0
#> 11 6.9274e-01 0.014137539 0
#> 12 7.2342e-01 0.008529625 0
#> 13 3.0443e-01 0.005905508 0
#> 14 5.4757e-01 0.004412045 0
#> 15 1.3600e-03 0.003461425 1

set.seed(1); Alpha_investing(sample.df, alpha=0.1, w0=0.025)
#>          pval      alphai R
#> 1  2.9000e-08 0.010818925 1
#> 2  6.7430e-02 0.041915265 0
#> 3  1.5140e-02 0.014226480 0
#> 4  8.1740e-02 0.007487045 0
#> 5  1.7100e-03 0.004738159 1
#> 6  2.7201e-01 0.046265410 0
#> 7  3.6000e-05 0.017453092 1
#> 8  7.9149e-01 0.057947646 0
#> 9  7.5900e-08 0.023879569 1
#> 10 2.8295e-01 0.063856912 0
#> 11 6.9274e-01 0.027880910 0
#> 12 7.2342e-01 0.016914972 0
#> 13 3.0443e-01 0.011741675 0
#> 14 5.4757e-01 0.008785330 0
#> 15 1.3600e-03 0.006898971 1
```
