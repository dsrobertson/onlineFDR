# Alpha-spending for online FWER control

Implements online FWER control using a Bonferroni-like test.

## Usage

``` r
Alpha_spending(
  d,
  alpha = 0.05,
  gammai,
  random = TRUE,
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

  Optional vector of \\\gamma_i\\, where hypothesis \\i\\ is rejected if
  the \\i\\-th p-value is less than or equal to \\\alpha \gamma_i\\. A
  default is provided as proposed by Javanmard and Montanari (2018),
  equation 31.

- random:

  Logical. If `TRUE` (the default), then the order of the p-values in
  each batch (i.e. those that have exactly the same date) is randomised.

- date.format:

  Optional string giving the format that is used for dates.

## Value

- out:

  A dataframe with the original data `d` (which will be reordered if
  there are batches and `random = TRUE`), the adjusted significance
  thresholds `alphai` and the indicator function of discoveries `R`,
  where `R[i] = 1` corresponds to hypothesis \\i\\ being rejected
  (otherwise `R[i] = 0`).

## Details

The function takes as its input either a vector of p-values, or a
dataframe with three columns: an identifier (\`id'), date (\`date') and
p-value (\`pval'). The case where p-values arrive in batches corresponds
to multiple instances of the same date. If no column of dates is
provided, then the p-values are treated as being ordered in sequence,
arriving one at a time.

Alpha-spending provides strong FWER control for a potentially infinite
stream of arbitrarily dependent p-values by using a Bonferroni-like
test. Given an overall significance level \\\alpha\\, we choose a
(potentially infinite) sequence of non-negative numbers \\\gamma_i\\
such that they sum to 1. Hypothesis \\i\\ is rejected if the \\i\\-th
p-value is less than or equal to \\\alpha \gamma_i\\.

Note that the procedure controls the generalised familywise error rate
(k-FWER) for \\k \> 1\\ if \\\alpha\\ is replaced by \\\min(1,
k\alpha)\\.

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
pval = c(2.90e-17, 0.06743, 0.01514, 0.08174, 0.00171,
        3.60e-05, 0.79149, 0.27201, 0.28295, 7.59e-08,
        0.69274, 0.30443, 0.00136, 0.72342, 0.54757))

set.seed(1); Alpha_spending(sample.df)
#>        id       date       pval       alphai R
#> 1  A15432 2014-12-01 2.9000e-17 0.0027407722 1
#> 2  B90969 2014-12-01 6.7430e-02 0.0005960313 0
#> 3  C18705 2014-12-01 1.5140e-02 0.0005076521 0
#> 4  B49731 2015-09-21 8.1740e-02 0.0004221825 0
#> 5  E99902 2015-09-21 1.7100e-03 0.0003579233 0
#> 6  D46627 2015-09-21 2.7201e-01 0.0003096307 0
#> 7  C38292 2015-09-21 3.6000e-05 0.0002724265 1
#> 8  A30619 2015-09-21 7.9149e-01 0.0002430188 0
#> 9  A41418 2016-05-19 7.5900e-08 0.0002192416 1
#> 10 E29198 2016-05-19 2.8295e-01 0.0001996425 0
#> 11 D51456 2016-11-12 6.9274e-01 0.0001832205 0
#> 12 A63155 2017-03-27 7.2342e-01 0.0001692671 0
#> 13 C88669 2017-03-27 3.0443e-01 0.0001572680 0
#> 14 B66033 2017-03-27 5.4757e-01 0.0001468417 0
#> 15 E03673 2017-03-27 1.3600e-03 0.0001376991 0

Alpha_spending(sample.df, random=FALSE)
#>        id       date       pval       alphai R
#> 1  A15432 2014-12-01 2.9000e-17 0.0027407722 1
#> 2  B90969 2014-12-01 6.7430e-02 0.0005960313 0
#> 3  C18705 2014-12-01 1.5140e-02 0.0005076521 0
#> 4  B49731 2015-09-21 8.1740e-02 0.0004221825 0
#> 5  E99902 2015-09-21 1.7100e-03 0.0003579233 0
#> 6  C38292 2015-09-21 3.6000e-05 0.0003096307 1
#> 7  A30619 2015-09-21 7.9149e-01 0.0002724265 0
#> 8  D46627 2015-09-21 2.7201e-01 0.0002430188 0
#> 9  E29198 2016-05-19 2.8295e-01 0.0002192416 0
#> 10 A41418 2016-05-19 7.5900e-08 0.0001996425 1
#> 11 D51456 2016-11-12 6.9274e-01 0.0001832205 0
#> 12 C88669 2017-03-27 3.0443e-01 0.0001692671 0
#> 13 E03673 2017-03-27 1.3600e-03 0.0001572680 0
#> 14 A63155 2017-03-27 7.2342e-01 0.0001468417 0
#> 15 B66033 2017-03-27 5.4757e-01 0.0001376991 0

set.seed(1); Alpha_spending(sample.df, alpha=0.1)
#>        id       date       pval       alphai R
#> 1  A15432 2014-12-01 2.9000e-17 0.0054815445 1
#> 2  B90969 2014-12-01 6.7430e-02 0.0011920626 0
#> 3  C18705 2014-12-01 1.5140e-02 0.0010153042 0
#> 4  B49731 2015-09-21 8.1740e-02 0.0008443651 0
#> 5  E99902 2015-09-21 1.7100e-03 0.0007158466 0
#> 6  D46627 2015-09-21 2.7201e-01 0.0006192614 0
#> 7  C38292 2015-09-21 3.6000e-05 0.0005448529 1
#> 8  A30619 2015-09-21 7.9149e-01 0.0004860376 0
#> 9  A41418 2016-05-19 7.5900e-08 0.0004384833 1
#> 10 E29198 2016-05-19 2.8295e-01 0.0003992849 0
#> 11 D51456 2016-11-12 6.9274e-01 0.0003664409 0
#> 12 A63155 2017-03-27 7.2342e-01 0.0003385341 0
#> 13 C88669 2017-03-27 3.0443e-01 0.0003145360 0
#> 14 B66033 2017-03-27 5.4757e-01 0.0002936833 0
#> 15 E03673 2017-03-27 1.3600e-03 0.0002753983 0

```
