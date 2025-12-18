# ADDIS: Adaptive discarding algorithm for online FDR control

Implements the ADDIS algorithm for online FDR control, where ADDIS
stands for an ADaptive algorithm that DIScards conservative nulls, as
presented by Tian and Ramdas (2019). The algorithm compensates for the
power loss of SAFFRON with conservative nulls, by including both
adaptivity in the fraction of null hypotheses (like SAFFRON) and the
conservativeness of nulls (unlike SAFFRON).

## Usage

``` r
ADDIS(
  d,
  alpha = 0.05,
  gammai,
  w0,
  lambda = 0.25,
  tau = 0.5,
  async = FALSE,
  random = TRUE,
  display_progress = FALSE,
  date.format = "%Y-%m-%d"
)
```

## Arguments

- d:

  Either a vector of p-values, or a dataframe with three columns: an
  identifier (\`id'), p-value (\`pval'), and decision times
  (\`decision.times').

- alpha:

  Overall significance level of the procedure, the default is 0.05.

- gammai:

  Optional vector of \\\gamma_i\\. A default is provided with
  \\\gamma_j\\ proportional to \\1/j^(1.6)\\.

- w0:

  Initial \`wealth' of the procedure, defaults to \\\alpha/2\\.

- lambda:

  Optional parameter that sets the threshold for \`candidate'
  hypotheses. Must be between 0 and tau, defaults to 0.25.

- tau:

  Optional threshold for hypotheses to be selected for testing. Must be
  between 0 and 1, defaults to 0.5.

- async:

  Logical. If `TRUE` runs the version for an asynchronous testing
  process. Defaults to FALSE.

- random:

  Logical. If `TRUE` (the default), then the order of the p-values in
  each batch (i.e. those that have exactly the same date) is randomised.
  Only needed if async=FALSE.

- display_progress:

  Logical. If `TRUE` prints out a progress bar for the algorithm
  runtime.

- date.format:

  Optional string giving the format that is used for dates.

## Value

- out:

  A dataframe with the original p-values `pval`, the adjusted testing
  levels \\\alpha_i\\ and the indicator function of discoveries `R`.
  Hypothesis \\i\\ is rejected if the \\i\\-th p-value is less than or
  equal to \\\alpha_i\\, in which case `R[i] = 1` (otherwise
  `R[i] = 0`).

## Details

The function takes as its input either a vector of p-values, or a
dataframe with three columns. The dataframe requires an identifier
(\`id'), date (\`date') and p-value (\`pval'). If the asynchronous
version is specified (see below), then the column date should be
replaced by the decision times.

Given an overall significance level \\\alpha\\, ADDIS depends on
constants \\w_0\\, \\\lambda\\ and \\\tau\\. Here \\w_0\\ represents the
initial \`wealth' of the procedure and satisfies \\0 \le w_0 \le
\alpha\\. \\\tau \in (0,1)\\ represents the threshold for a hypothesis
to be selected for testing: p-values greater than \\\tau\\ are
implicitly \`discarded' by the procedure. Finally, \\\lambda \in
\[0,\tau)\\ sets the threshold for a p-value to be a candidate for
rejection: ADDIS will never reject a p-value larger than \\\lambda\\.
The algorithm also require a sequence of non-negative non-increasing
numbers \\\gamma_i\\ that sum to 1.

The ADDIS procedure provably controls the FDR for independent p-values.
Tian and Ramdas (2019) also presented a version for an asynchronous
testing process, consisting of tests that start and finish at
(potentially) random times. The discretised finish times of the test
correspond to the decision times. These decision times are given as the
input `decision.times`. Note that this asynchronous version controls a
modified version of the FDR.

Further details of the ADDIS algorithms can be found in Tian and Ramdas
(2019).

## References

Tian, J. and Ramdas, A. (2019). ADDIS: an adaptive discarding algorithm
for online FDR control with conservative nulls. *Advances in Neural
Information Processing Systems*, 9388-9396.

## See also

ADDIS is identical to
[`SAFFRON`](https://dsrobertson.github.io/onlineFDR/reference/SAFFRON.md)
with option `discard=TRUE`.

ADDIS with option `async=TRUE` is identical to
[`SAFFRONstar`](https://dsrobertson.github.io/onlineFDR/reference/SAFFRONstar.md)
with option `discard=TRUE`.

## Examples

``` r
sample.df1 <- data.frame(
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

ADDIS(sample.df1, random=FALSE)
#>          pval      alphai R     id
#> 1  2.9000e-08 0.002734314 1 A15432
#> 2  6.7430e-02 0.005468627 0 B90969
#> 3  1.5140e-02 0.005468627 0 C18705
#> 4  8.1740e-02 0.005468627 0 B49731
#> 5  1.7100e-03 0.005468627 1 E99902
#> 6  3.6000e-05 0.010937254 1 C38292
#> 7  7.9149e-01 0.016405881 0 A30619
#> 8  2.7201e-01 0.016405881 0 D46627
#> 9  2.8295e-01 0.005411923 0 E29198
#> 10 7.5900e-08 0.002828822 1 A41418
#> 11 6.9274e-01 0.008297449 0 D51456
#> 12 3.0443e-01 0.008297449 0 C88669
#> 13 1.3600e-03 0.003589243 1 E03673
#> 14 7.2342e-01 0.009057870 0 A63155
#> 15 5.4757e-01 0.009057870 0 B66033


sample.df2 <- data.frame(
id = c('A15432', 'B90969', 'C18705', 'B49731', 'E99902',
    'C38292', 'A30619', 'D46627', 'E29198', 'A41418',
    'D51456', 'C88669', 'E03673', 'A63155', 'B66033'),
pval = c(2.90e-08, 0.06743, 0.01514, 0.08174, 0.00171,
        3.60e-05, 0.79149, 0.27201, 0.28295, 7.59e-08,
        0.69274, 0.30443, 0.00136, 0.72342, 0.54757),
decision.times = seq_len(15) + 1)

ADDIS(sample.df2, async = TRUE) # Asynchronous
#>          pval       alphai R     id
#> 1  2.9000e-08 0.0027343135 1 A15432
#> 2  6.7430e-02 0.0009019871 0 B90969
#> 3  1.5140e-02 0.0018039742 0 C18705
#> 4  8.1740e-02 0.0018039742 0 B49731
#> 5  1.7100e-03 0.0018039742 1 E99902
#> 6  3.6000e-05 0.0018039742 1 C38292
#> 7  7.9149e-01 0.0036079483 0 A30619
#> 8  2.7201e-01 0.0054119225 0 D46627
#> 9  2.8295e-01 0.0054119225 0 E29198
#> 10 7.5900e-08 0.0028288216 1 A41418
#> 11 6.9274e-01 0.0017852686 0 D51456
#> 12 3.0443e-01 0.0035892428 0 C88669
#> 13 1.3600e-03 0.0035892428 1 E03673
#> 14 7.2342e-01 0.0021921853 0 A63155
#> 15 5.4757e-01 0.0039961595 0 B66033

```
