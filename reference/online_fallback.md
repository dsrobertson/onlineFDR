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

  Overall significance level of the FWER procedure, the default is 0.05.

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

  A dataframe with the p-values `pval` (reordered if there are batches
  and `random = TRUE`), the adjusted significance thresholds
  \\\alpha_i\\, the indicator function of discoveries `R` and, if `d` is
  a dataframe with a column `id`, the identifiers `id`. Hypothesis \\i\\
  is rejected if the \\i\\-th p-value is less than or equal to
  \\\alpha_i\\, in which case `R[i] = 1` (otherwise `R[i] = 0`).

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
level \$\$\alpha_i = \alpha \gamma_i + R\_{i-1} \alpha\_{i-1}\$\$ for
\\i \ge 2\\, with \\\alpha_1 = \alpha\gamma_1\\, where \\R_i = 1\\p_i
\leq \alpha_i\\\\ denotes a rejected hypothesis.

Further details of the online fallback procedure can be found in Tian
and Ramdas (2021).

## References

Javanmard, A. and Montanari, A. (2018). Online Rules for Control of
False Discovery Rate and False Discovery Exceedance. *Annals of
Statistics*, 46(2):526-554.

Tian, J. and Ramdas, A. (2021). Online control of the familywise error
rate. *Statistical Methods in Medical Research*, 30(4):976-993.

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
#>          pval       alphai R     id
#> 1  2.9000e-08 0.0027407722 1 A15432
#> 2  6.7430e-02 0.0033368035 0 B90969
#> 3  1.5140e-02 0.0005076521 0 C18705
#> 4  8.1740e-02 0.0004221825 0 B49731
#> 5  1.7100e-03 0.0003579233 0 E99902
#> 6  3.6000e-05 0.0003096307 1 C38292
#> 7  7.9149e-01 0.0005820572 0 A30619
#> 8  2.7201e-01 0.0002430188 0 D46627
#> 9  2.8295e-01 0.0002192416 0 E29198
#> 10 7.5900e-08 0.0001996425 1 A41418
#> 11 6.9274e-01 0.0003828629 0 D51456
#> 12 3.0443e-01 0.0001692671 0 C88669
#> 13 1.3600e-03 0.0001572680 0 E03673
#> 14 7.2342e-01 0.0001468417 0 A63155
#> 15 5.4757e-01 0.0001376991 0 B66033

set.seed(1); online_fallback(sample.df)
#>          pval       alphai R     id
#> 1  2.9000e-08 0.0027407722 1 A15432
#> 2  6.7430e-02 0.0033368035 0 B90969
#> 3  1.5140e-02 0.0005076521 0 C18705
#> 4  8.1740e-02 0.0004221825 0 B49731
#> 5  1.7100e-03 0.0003579233 0 E99902
#> 6  2.7201e-01 0.0003096307 0 D46627
#> 7  3.6000e-05 0.0002724265 1 C38292
#> 8  7.9149e-01 0.0005154453 0 A30619
#> 9  7.5900e-08 0.0002192416 1 A41418
#> 10 2.8295e-01 0.0004188841 0 E29198
#> 11 6.9274e-01 0.0001832205 0 D51456
#> 12 7.2342e-01 0.0001692671 0 A63155
#> 13 3.0443e-01 0.0001572680 0 C88669
#> 14 5.4757e-01 0.0001468417 0 B66033
#> 15 1.3600e-03 0.0001376991 0 E03673

set.seed(1); online_fallback(sample.df, alpha=0.1)
#>          pval       alphai R     id
#> 1  2.9000e-08 0.0054815445 1 A15432
#> 2  6.7430e-02 0.0066736071 0 B90969
#> 3  1.5140e-02 0.0010153042 0 C18705
#> 4  8.1740e-02 0.0008443651 0 B49731
#> 5  1.7100e-03 0.0007158466 0 E99902
#> 6  2.7201e-01 0.0006192614 0 D46627
#> 7  3.6000e-05 0.0005448529 1 C38292
#> 8  7.9149e-01 0.0010308906 0 A30619
#> 9  7.5900e-08 0.0004384833 1 A41418
#> 10 2.8295e-01 0.0008377682 0 E29198
#> 11 6.9274e-01 0.0003664409 0 D51456
#> 12 7.2342e-01 0.0003385341 0 A63155
#> 13 3.0443e-01 0.0003145360 0 C88669
#> 14 5.4757e-01 0.0002936833 0 B66033
#> 15 1.3600e-03 0.0002753983 0 E03673
```
