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

  Overall significance level of the FDR procedure, the default is 0.05.

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
  there are batches and `random = TRUE`), the adjusted signifcance
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
stream of p-values by using a Bonferroni-like test. Given an overall
significance level \\\alpha\\, we choose a (potentially infinite)
sequence of non-negative numbers \\\gamma_i\\ such that they sum to 1.
Hypothesis \\i\\ is rejected if the \\i\\-th p-value is less than or
equal to \\\alpha \gamma_i\\.

Note that the procedure controls the generalised familywise error rate
(k-FWER) for \\k \> 1\\ if \\\alpha\\ is replaced by min(\\1,
k\alpha\\).

## References

Javanmard, A. and Montanari, A. (2018) Online Rules for Control of False
Discovery Rate and False Discovery Exceedance. *Annals of Statistics*,
46(2):526-554.

Tian, J. and Ramdas, A. (2021). Online control of the familywise error
rate. *Statistical Methods for Medical Research* (to appear),
<https://arxiv.org/abs/1910.04900>.

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
#> 1  A15432 2014-12-01 2.9000e-17 0.0026758385 1
#> 2  B90969 2014-12-01 6.7430e-02 0.0005819103 0
#> 3  C18705 2014-12-01 1.5140e-02 0.0004956249 0
#> 4  B49731 2015-09-21 8.1740e-02 0.0004121803 0
#> 5  E99902 2015-09-21 1.7100e-03 0.0003494435 0
#> 6  D46627 2015-09-21 2.7201e-01 0.0003022950 0
#> 7  C38292 2015-09-21 3.6000e-05 0.0002659722 1
#> 8  A30619 2015-09-21 7.9149e-01 0.0002372613 0
#> 9  A41418 2016-05-19 7.5900e-08 0.0002140474 1
#> 10 E29198 2016-05-19 2.8295e-01 0.0001949126 0
#> 11 D51456 2016-11-12 6.9274e-01 0.0001788796 0
#> 12 A63155 2017-03-27 7.2342e-01 0.0001652568 0
#> 13 C88669 2017-03-27 3.0443e-01 0.0001535420 0
#> 14 B66033 2017-03-27 5.4757e-01 0.0001433627 0
#> 15 E03673 2017-03-27 1.3600e-03 0.0001344368 0

Alpha_spending(sample.df, random=FALSE)
#>        id       date       pval       alphai R
#> 1  A15432 2014-12-01 2.9000e-17 0.0026758385 1
#> 2  B90969 2014-12-01 6.7430e-02 0.0005819103 0
#> 3  C18705 2014-12-01 1.5140e-02 0.0004956249 0
#> 4  B49731 2015-09-21 8.1740e-02 0.0004121803 0
#> 5  E99902 2015-09-21 1.7100e-03 0.0003494435 0
#> 6  C38292 2015-09-21 3.6000e-05 0.0003022950 1
#> 7  A30619 2015-09-21 7.9149e-01 0.0002659722 0
#> 8  D46627 2015-09-21 2.7201e-01 0.0002372613 0
#> 9  E29198 2016-05-19 2.8295e-01 0.0002140474 0
#> 10 A41418 2016-05-19 7.5900e-08 0.0001949126 1
#> 11 D51456 2016-11-12 6.9274e-01 0.0001788796 0
#> 12 C88669 2017-03-27 3.0443e-01 0.0001652568 0
#> 13 E03673 2017-03-27 1.3600e-03 0.0001535420 0
#> 14 A63155 2017-03-27 7.2342e-01 0.0001433627 0
#> 15 B66033 2017-03-27 5.4757e-01 0.0001344368 0

set.seed(1); Alpha_spending(sample.df, alpha=0.1)
#>        id       date       pval       alphai R
#> 1  A15432 2014-12-01 2.9000e-17 0.0053516771 1
#> 2  B90969 2014-12-01 6.7430e-02 0.0011638206 0
#> 3  C18705 2014-12-01 1.5140e-02 0.0009912499 0
#> 4  B49731 2015-09-21 8.1740e-02 0.0008243606 0
#> 5  E99902 2015-09-21 1.7100e-03 0.0006988870 0
#> 6  D46627 2015-09-21 2.7201e-01 0.0006045900 0
#> 7  C38292 2015-09-21 3.6000e-05 0.0005319444 1
#> 8  A30619 2015-09-21 7.9149e-01 0.0004745225 0
#> 9  A41418 2016-05-19 7.5900e-08 0.0004280949 1
#> 10 E29198 2016-05-19 2.8295e-01 0.0003898252 0
#> 11 D51456 2016-11-12 6.9274e-01 0.0003577593 0
#> 12 A63155 2017-03-27 7.2342e-01 0.0003305137 0
#> 13 C88669 2017-03-27 3.0443e-01 0.0003070841 0
#> 14 B66033 2017-03-27 5.4757e-01 0.0002867254 0
#> 15 E03673 2017-03-27 1.3600e-03 0.0002688736 0

```
