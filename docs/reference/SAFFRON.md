# SAFFRON: Adaptive online FDR control

Implements the SAFFRON procedure for online FDR control, where SAFFRON
stands for Serial estimate of the Alpha Fraction that is Futilely
Rationed On true Null hypotheses, as presented by Ramdas et al. (2018).
The algorithm is based on an estimate of the proportion of true null
hypotheses. More precisely, SAFFRON sets the adjusted test levels based
on an estimate of the amount of alpha-wealth that is allocated to
testing the true null hypotheses.

## Usage

``` r
SAFFRON(
  d,
  alpha = 0.05,
  gammai,
  w0,
  lambda = 0.5,
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

  Optional vector of \\\gamma_i\\. A default is provided with
  \\\gamma_j\\ proportional to \\1/j^(1.6)\\.

- w0:

  Initial \`wealth' of the procedure, defaults to \\\alpha/2\\. Must be
  between 0 and \\\alpha\\.

- lambda:

  Optional threshold for a \`candidate' hypothesis, must be between 0
  and 1. Defaults to 0.5.

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
arriving one at a time.

SAFFRON procedure provably controls FDR for independent p-values. Given
an overall significance level \\\alpha\\, we choose a sequence of
non-negative non-increasing numbers \\\gamma_i\\ that sum to 1.

SAFFRON depends on constants \\w_0\\ and \\\lambda\\, where \\w_0\\
satisfies \\0 \le w_0 \le \alpha\\ and represents the initial \`wealth'
of the procedure, and \\0 \< \lambda \< 1\\ represents the threshold for
a \`candidate' hypothesis. A \`candidate' refers to p-values smaller
than \\\lambda\\, since SAFFRON will never reject a p-value larger than
\\\lambda\\.

Note that FDR control also holds for the SAFFRON procedure if only the
p-values corresponding to true nulls are mutually independent, and
independent from the non-null p-values.

The SAFFRON procedure can lose power in the presence of conservative
nulls, which can be compensated for by adaptively \`discarding' these
p-values. This option is called by setting `discard=TRUE`, which is the
same algorithm as ADDIS.

Further details of the SAFFRON procedure can be found in Ramdas et al.
(2018).

## References

Ramdas, A., Zrnic, T., Wainwright M.J. and Jordan, M.I. (2018). SAFFRON:
an adaptive algorithm for online control of the false discovery rate.
*Proceedings of the 35th International Conference in Machine Learning*,
80:4286-4294.

## See also

[`SAFFRONstar`](https://dsrobertson.github.io/onlineFDR/reference/SAFFRONstar.md)
presents versions of SAFFRON for *asynchronous* testing, i.e. where each
hypothesis test can itself be a sequential process and the tests can
overlap in time.

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

SAFFRON(sample.df, random=FALSE)
#>          pval      alphai R     id
#> 1  2.9000e-08 0.005468627 1 A15432
#> 2  6.7430e-02 0.010937254 0 B90969
#> 3  1.5140e-02 0.010937254 0 C18705
#> 4  8.1740e-02 0.010937254 0 B49731
#> 5  1.7100e-03 0.010937254 1 E99902
#> 6  3.6000e-05 0.021874508 1 C38292
#> 7  7.9149e-01 0.032811762 0 A30619
#> 8  2.7201e-01 0.010823845 0 D46627
#> 9  2.8295e-01 0.010823845 0 E29198
#> 10 7.5900e-08 0.010823845 1 A41418
#> 11 6.9274e-01 0.021761099 0 D51456
#> 12 3.0443e-01 0.009265591 0 C88669
#> 13 1.3600e-03 0.009265591 1 E03673
#> 14 7.2342e-01 0.020202846 0 A63155
#> 15 5.4757e-01 0.009064367 0 B66033

set.seed(1); SAFFRON(sample.df)
#>          pval      alphai R     id
#> 1  2.9000e-08 0.005468627 1 A15432
#> 2  6.7430e-02 0.010937254 0 B90969
#> 3  1.5140e-02 0.010937254 0 C18705
#> 4  8.1740e-02 0.010937254 0 B49731
#> 5  1.7100e-03 0.010937254 1 E99902
#> 6  2.7201e-01 0.021874508 0 D46627
#> 7  3.6000e-05 0.021874508 1 C38292
#> 8  7.9149e-01 0.032811762 0 A30619
#> 9  7.5900e-08 0.010823845 1 A41418
#> 10 2.8295e-01 0.021761099 0 E29198
#> 11 6.9274e-01 0.021761099 0 D51456
#> 12 7.2342e-01 0.009265591 0 A63155
#> 13 3.0443e-01 0.005456418 0 C88669
#> 14 5.4757e-01 0.005456418 0 B66033
#> 15 1.3600e-03 0.003688669 1 E03673

set.seed(1); SAFFRON(sample.df, alpha=0.1, w0=0.025)
#>          pval      alphai R     id
#> 1  2.9000e-08 0.005468627 1 A15432
#> 2  6.7430e-02 0.021874508 0 B90969
#> 3  1.5140e-02 0.021874508 1 C18705
#> 4  8.1740e-02 0.043749017 0 B49731
#> 5  1.7100e-03 0.043749017 1 E99902
#> 6  2.7201e-01 0.065623525 0 D46627
#> 7  3.6000e-05 0.065623525 1 C38292
#> 8  7.9149e-01 0.087498033 0 A30619
#> 9  7.5900e-08 0.028863587 1 A41418
#> 10 2.8295e-01 0.050738095 0 E29198
#> 11 6.9274e-01 0.050738095 0 D51456
#> 12 7.2342e-01 0.022302945 0 A63155
#> 13 3.0443e-01 0.013293195 0 C88669
#> 14 5.4757e-01 0.013293195 0 B66033
#> 15 1.3600e-03 0.009042997 1 E03673

```
