# LONDstar: Asynchronous online mFDR control based on number of discoveries

Implements the LOND algorithm for asynchronous online testing, as
presented by Zrnic et al. (2021).

## Usage

``` r
LONDstar(
  d,
  alpha = 0.05,
  version,
  betai,
  batch.sizes,
  display_progress = FALSE
)
```

## Arguments

- d:

  Either a vector of p-values, or a dataframe with three columns: an
  identifier (\`id'), p-value (\`pval'), and either \`decision.times',
  or \`lags', depending on which version you're using. See version for
  more details.

- alpha:

  Overall significance level of the procedure, the default is 0.05.

- version:

  Takes values 'async', 'dep' or 'batch'. This specifies the version of
  LONDstar to use. `version='async'` requires a column of decision times
  (\`decision.times'). `version='dep'` requires a column of lags
  (\`lags'). `version='batch'` requires a vector of batch sizes
  (\`batch.sizes').

- betai:

  Optional vector of \\\beta_i\\. A default is provided as proposed by
  Javanmard and Montanari (2018), equation 31.

- batch.sizes:

  A vector of batch sizes, this is required for `version='batch'`.

- display_progress:

  Logical. If `TRUE` prints out a progress bar for the algorithm
  runtime.

## Value

- out:

  A dataframe with the original p-values `pval`, the adjusted testing
  levels \\\alpha_i\\ and the indicator function of discoveries `R`.
  Hypothesis \\i\\ is rejected if the \\i\\-th p-value is less than or
  equal to \\\alpha_i\\, in which case `R[i] = 1` (otherwise
  `R[i] = 0`).

## Details

The function takes as its input either a vector of p-values, or a
dataframe with three columns: an identifier (\`id'), p-value (\`pval'),
or a column describing the conflict sets for the hypotheses. This takes
the form of a vector of decision times or lags. Batch sizes can be
specified as a separate argument (see below).

Zrnic et al. (2021) present explicit three versions of LONDstar:

1\) `version='async'` is for an asynchronous testing process, consisting
of tests that start and finish at (potentially) random times. The
discretised finish times of the test correspond to the decision times.
These decision times are given as the input `decision.times` for this
version of the LONDstar algorithm.

2\) `version='dep'` is for online testing under local dependence of the
p-values. More precisely, for any \\t\>0\\ we allow the p-value \\p_t\\
to have arbitrary dependence on the previous \\L_t\\ p-values. The fixed
sequence \\L_t\\ is referred to as \`lags', and is given as the input
`lags` for this version of the LONDstar algorithm.

3\) `version='batch'` is for controlling the mFDR in mini-batch testing,
where a mini-batch represents a grouping of tests run asynchronously
which result in dependent p-values. Once a mini-batch of tests is fully
completed, a new one can start, testing hypotheses independent of the
previous batch. The batch sizes are given as the input `batch.sizes` for
this version of the LONDstar algorithm.

Given an overall significance level \\\alpha\\, LONDstar requires a
sequence of non-negative non-increasing numbers \\\beta_i\\ that sum to
\\\alpha\\.

Note that these LONDstar algorithms control the *modified* FDR (mFDR).

Further details of the LONDstar algorithms can be found in Zrnic et al.
(2021).

## References

Javanmard, A. and Montanari, A. (2018) Online Rules for Control of False
Discovery Rate and False Discovery Exceedance. *Annals of Statistics*,
46(2):526-554.

Zrnic, T., Ramdas, A. and Jordan, M.I. (2021). Asynchronous Online
Testing of Multiple Hypotheses. *Journal of Machine Learning Research*
(to appear), <https://arxiv.org/abs/1812.05068>.

## See also

[`LOND`](https://dsrobertson.github.io/onlineFDR/reference/LOND.md)
presents versions of LOND for *synchronous* p-values, i.e. where each
test can only start when the previous test has finished.

## Examples

``` r
sample.df <- data.frame(
id = c('A15432', 'B90969', 'C18705', 'B49731', 'E99902',
    'C38292', 'A30619', 'D46627', 'E29198', 'A41418',
    'D51456', 'C88669', 'E03673', 'A63155', 'B66033'),
pval = c(2.90e-08, 0.06743, 0.01514, 0.08174, 0.00171,
        3.60e-05, 0.79149, 0.27201, 0.28295, 7.59e-08,
        0.69274, 0.30443, 0.00136, 0.72342, 0.54757),
decision.times = seq_len(15) + 1)

LONDstar(sample.df, version='async')
#>          pval       alphai R
#> 1  2.9000e-08 0.0026758385 1
#> 2  6.7430e-02 0.0005819103 0
#> 3  1.5140e-02 0.0004956249 0
#> 4  8.1740e-02 0.0004121803 0
#> 5  1.7100e-03 0.0003494435 0
#> 6  3.6000e-05 0.0003022950 1
#> 7  7.9149e-01 0.0002659722 0
#> 8  2.7201e-01 0.0004745225 0
#> 9  2.8295e-01 0.0004280949 0
#> 10 7.5900e-08 0.0003898252 1
#> 11 6.9274e-01 0.0003577593 0
#> 12 3.0443e-01 0.0004957705 0
#> 13 1.3600e-03 0.0004606261 0
#> 14 7.2342e-01 0.0004300881 0
#> 15 5.4757e-01 0.0004033104 0

sample.df2 <- data.frame(
id = c('A15432', 'B90969', 'C18705', 'B49731', 'E99902',
    'C38292', 'A30619', 'D46627', 'E29198', 'A41418',
    'D51456', 'C88669', 'E03673', 'A63155', 'B66033'),
pval = c(2.90e-08, 0.06743, 0.01514, 0.08174, 0.00171,
        3.60e-05, 0.79149, 0.27201, 0.28295, 7.59e-08,
        0.69274, 0.30443, 0.00136, 0.72342, 0.54757),
lags = rep(1,15))

LONDstar(sample.df2, version='dep')
#>          pval lag       alphai R
#> 1  2.9000e-08   1 0.0026758385 1
#> 2  6.7430e-02   1 0.0005819103 0
#> 3  1.5140e-02   1 0.0004956249 0
#> 4  8.1740e-02   1 0.0004121803 0
#> 5  1.7100e-03   1 0.0003494435 0
#> 6  3.6000e-05   1 0.0003022950 1
#> 7  7.9149e-01   1 0.0002659722 0
#> 8  2.7201e-01   1 0.0004745225 0
#> 9  2.8295e-01   1 0.0004280949 0
#> 10 7.5900e-08   1 0.0003898252 1
#> 11 6.9274e-01   1 0.0003577593 0
#> 12 3.0443e-01   1 0.0004957705 0
#> 13 1.3600e-03   1 0.0004606261 0
#> 14 7.2342e-01   1 0.0004300881 0
#> 15 5.4757e-01   1 0.0004033104 0

sample.df3 <- data.frame(
id = c('A15432', 'B90969', 'C18705', 'B49731', 'E99902',
    'C38292', 'A30619', 'D46627', 'E29198', 'A41418',
    'D51456', 'C88669', 'E03673', 'A63155', 'B66033'),
pval = c(2.90e-08, 0.06743, 0.01514, 0.08174, 0.00171,
        3.60e-05, 0.79149, 0.27201, 0.28295, 7.59e-08,
        0.69274, 0.30443, 0.00136, 0.72342, 0.54757))

LONDstar(sample.df3, version='batch', batch.sizes = c(4,6,5))
```
