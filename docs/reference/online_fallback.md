# Online fallback procedure for FWER control

Implements the online fallback procedure of Tian and Ramdas (2021),
which guarantees strong FWER control under arbitrary dependence of the
p-values.

## Usage

``` r
online_fallback(
  d,
  alpha = 0.05,
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

- alpha:

  Overall significance level of the FDR procedure, the default is 0.05.

- gammai:

  Optional vector of \\\gamma_i\\. A default is provided as proposed by
  Javanmard and Montanari (2018), equation 31.

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
provided, then the p-values are treated as being ordered in sequence,
arriving one at a time. Given an overall significance level \\\alpha\\,
we choose a sequence of non-negative non-increasing numbers \\\gamma_i\\
that sum to 1.

The online fallback procedure provides a uniformly more powerful method
than Alpha-spending, by saving the significance level of a previous
rejection. More specifically, the procedure tests hypothesis \\H_i\\ at
level \$\$\alpha_i = \alpha \gamma_i + R\_{i-1} \alpha\_{i-1}\$\$ where
\\R_i = 1\\p_i \leq \alpha_i\\\\ denotes a rejected hypothesis.

Further details of the online fallback procedure can be found in Tian
and Ramdas (2021).

## References

Tian, J. and Ramdas, A. (2021). Online control of the familywise error
rate. *Statistical Methods for Medical Research*, 30(4):976–993.

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

online_fallback(sample.df, random=FALSE)
#>          pval       alphai R
#> 1  2.9000e-08 0.0026758385 1
#> 2  6.7430e-02 0.0032577488 0
#> 3  1.5140e-02 0.0004956249 0
#> 4  8.1740e-02 0.0004121803 0
#> 5  1.7100e-03 0.0003494435 0
#> 6  3.6000e-05 0.0003022950 1
#> 7  7.9149e-01 0.0005682672 0
#> 8  2.7201e-01 0.0002372613 0
#> 9  2.8295e-01 0.0002140474 0
#> 10 7.5900e-08 0.0001949126 1
#> 11 6.9274e-01 0.0003737922 0
#> 12 3.0443e-01 0.0001652568 0
#> 13 1.3600e-03 0.0001535420 0
#> 14 7.2342e-01 0.0001433627 0
#> 15 5.4757e-01 0.0001344368 0

set.seed(1); online_fallback(sample.df)
#>          pval       alphai R
#> 1  2.9000e-08 0.0026758385 1
#> 2  6.7430e-02 0.0032577488 0
#> 3  1.5140e-02 0.0004956249 0
#> 4  8.1740e-02 0.0004121803 0
#> 5  1.7100e-03 0.0003494435 0
#> 6  2.7201e-01 0.0003022950 0
#> 7  3.6000e-05 0.0002659722 1
#> 8  7.9149e-01 0.0005032335 0
#> 9  7.5900e-08 0.0002140474 1
#> 10 2.8295e-01 0.0004089600 0
#> 11 6.9274e-01 0.0001788796 0
#> 12 7.2342e-01 0.0001652568 0
#> 13 3.0443e-01 0.0001535420 0
#> 14 5.4757e-01 0.0001433627 0
#> 15 1.3600e-03 0.0001344368 0

set.seed(1); online_fallback(sample.df, alpha=0.1)
#>          pval       alphai R
#> 1  2.9000e-08 0.0053516771 1
#> 2  6.7430e-02 0.0065154977 0
#> 3  1.5140e-02 0.0009912499 0
#> 4  8.1740e-02 0.0008243606 0
#> 5  1.7100e-03 0.0006988870 0
#> 6  2.7201e-01 0.0006045900 0
#> 7  3.6000e-05 0.0005319444 1
#> 8  7.9149e-01 0.0010064670 0
#> 9  7.5900e-08 0.0004280949 1
#> 10 2.8295e-01 0.0008179201 0
#> 11 6.9274e-01 0.0003577593 0
#> 12 7.2342e-01 0.0003305137 0
#> 13 3.0443e-01 0.0003070841 0
#> 14 5.4757e-01 0.0002867254 0
#> 15 1.3600e-03 0.0002688736 0
```
