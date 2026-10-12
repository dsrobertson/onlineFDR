# LOND: Online FDR control based on number of discoveries

Implements the LOND algorithm for online FDR control, where LOND stands
for (significance) Levels based On Number of Discoveries, as presented
by Javanmard and Montanari (2015).

## Usage

``` r
LOND(
  d,
  alpha = 0.05,
  betai,
  dep = FALSE,
  random = TRUE,
  display_progress = FALSE,
  date.format = "%Y-%m-%d",
  original = TRUE
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

- betai:

  Optional vector of \\\beta_i\\. A default is provided as proposed by
  Javanmard and Montanari (2018), equation 31.

- dep:

  Logical. If `TRUE`, runs the modified LOND algorithm which guarantees
  FDR control for *dependent* p-values. Defaults to `FALSE`.

- random:

  Logical. If `TRUE` (the default), then the order of the p-values in
  each batch (i.e. those that have exactly the same date) is randomised.

- display_progress:

  Logical. If `TRUE` prints out a progress bar for the algorithm
  runtime.

- date.format:

  Optional string giving the format that is used for dates.

- original:

  Logical. If `TRUE`, runs the original LOND algorithm of Javanmard and
  Montanari (2015), otherwise runs the modified algorithm with
  \\\alpha_i = \max(D(i-1), 1)\beta_i\\ (Javanmard and Montanari, 2015,
  Remark 2.3), which also controls the FDR under positive dependence
  (Zrnic et al., 2021). Defaults to `TRUE`.

## Value

- out:

  A dataframe with the p-values `pval` (reordered if there are batches
  and `random = TRUE`), the LOND-adjusted significance thresholds
  \\\alpha_i\\, the indicator function of discoveries `R` and, if `d` is
  a dataframe with a column `id`, the identifiers `id`. Hypothesis \\i\\
  is rejected if the \\i\\-th p-value is less than or equal to
  \\\alpha_i\\, in which case `R[i] = 1` (otherwise `R[i] = 0`).

## Details

The function takes as its input either a vector of p-values, or a
dataframe with three columns: an identifier (\`id'), date (\`date') and
p-value (\`pval'). The case where p-values arrive in batches corresponds
to multiple instances of the same date. If no column of dates is
provided, then the p-values are treated as being ordered in sequence,
arriving one at a time.

The LOND algorithm controls the FDR for independent p-values (see below
for the modification for dependent p-values). Given an overall
significance level \\\alpha\\, we choose a sequence of non-negative
numbers \\\beta_i\\ such that they sum to \\\alpha\\. The values of the
adjusted significance thresholds \\\alpha_i\\ are chosen as follows:
\$\$\alpha_i = (D(i-1) + 1)\beta_i\$\$ where \\D(n)\\ denotes the number
of discoveries in the first \\n\\ hypotheses.

A slightly modified version of LOND with thresholds \\\alpha_i =
\max(D(i-1), 1)\beta_i\\ (Javanmard and Montanari, 2015, Remark 2.3)
provably controls the FDR under positive dependence (PRDS condition),
see Zrnic et al. (2021).

For arbitrarily dependent p-values, LOND controls the FDR if it is
modified with \\\beta_i / H(i)\\ in place of \\\beta_i\\, where \\H(i)\\
is the i-th harmonic number.

Further details of the LOND algorithm can be found in Javanmard and
Montanari (2015).

## References

Javanmard, A. and Montanari, A. (2015). On Online Control of False
Discovery Rate. *arXiv preprint*, <https://arxiv.org/abs/1502.06197>.

Javanmard, A. and Montanari, A. (2018). Online Rules for Control of
False Discovery Rate and False Discovery Exceedance. *Annals of
Statistics*, 46(2):526-554.

Zrnic, T., Ramdas, A. and Jordan, M.I. (2021). Asynchronous Online
Testing of Multiple Hypotheses. *Journal of Machine Learning Research*,
22(33):1-39.

## See also

[`LONDstar`](https://dsrobertson.github.io/onlineFDR/reference/LONDstar.md)
presents versions of LOND for *asynchronous* testing, i.e. where each
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

set.seed(1); LOND(sample.df)
#>          pval       alphai R     id
#> 1  2.9000e-08 0.0027407722 1 A15432
#> 2  6.7430e-02 0.0011920626 0 B90969
#> 3  1.5140e-02 0.0010153042 0 C18705
#> 4  8.1740e-02 0.0008443651 0 B49731
#> 5  1.7100e-03 0.0007158466 0 E99902
#> 6  2.7201e-01 0.0006192614 0 D46627
#> 7  3.6000e-05 0.0005448529 1 C38292
#> 8  7.9149e-01 0.0007290564 0 A30619
#> 9  7.5900e-08 0.0006577249 1 A41418
#> 10 2.8295e-01 0.0007985699 0 E29198
#> 11 6.9274e-01 0.0007328818 0 D51456
#> 12 7.2342e-01 0.0006770682 0 A63155
#> 13 3.0443e-01 0.0006290720 0 C88669
#> 14 5.4757e-01 0.0005873666 0 B66033
#> 15 1.3600e-03 0.0005507966 0 E03673

LOND(sample.df, random=FALSE)
#>          pval       alphai R     id
#> 1  2.9000e-08 0.0027407722 1 A15432
#> 2  6.7430e-02 0.0011920626 0 B90969
#> 3  1.5140e-02 0.0010153042 0 C18705
#> 4  8.1740e-02 0.0008443651 0 B49731
#> 5  1.7100e-03 0.0007158466 0 E99902
#> 6  3.6000e-05 0.0006192614 1 C38292
#> 7  7.9149e-01 0.0008172794 0 A30619
#> 8  2.7201e-01 0.0007290564 0 D46627
#> 9  2.8295e-01 0.0006577249 0 E29198
#> 10 7.5900e-08 0.0005989274 1 A41418
#> 11 6.9274e-01 0.0007328818 0 D51456
#> 12 3.0443e-01 0.0006770682 0 C88669
#> 13 1.3600e-03 0.0006290720 0 E03673
#> 14 7.2342e-01 0.0005873666 0 A63155
#> 15 5.4757e-01 0.0005507966 0 B66033

set.seed(1); LOND(sample.df, alpha=0.1)
#>          pval      alphai R     id
#> 1  2.9000e-08 0.005481544 1 A15432
#> 2  6.7430e-02 0.002384125 0 B90969
#> 3  1.5140e-02 0.002030608 0 C18705
#> 4  8.1740e-02 0.001688730 0 B49731
#> 5  1.7100e-03 0.001431693 0 E99902
#> 6  2.7201e-01 0.001238523 0 D46627
#> 7  3.6000e-05 0.001089706 1 C38292
#> 8  7.9149e-01 0.001458113 0 A30619
#> 9  7.5900e-08 0.001315450 1 A41418
#> 10 2.8295e-01 0.001597140 0 E29198
#> 11 6.9274e-01 0.001465764 0 D51456
#> 12 7.2342e-01 0.001354136 0 A63155
#> 13 3.0443e-01 0.001258144 0 C88669
#> 14 5.4757e-01 0.001174733 0 B66033
#> 15 1.3600e-03 0.001101593 0 E03673
```
