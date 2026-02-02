# LORDstar: Asynchronous online mFDR control based on recent discovery

Implements LORD algorithms for asynchronous online testing, as presented
by Zrnic et al. (2021).

## Usage

``` r
LORDstar(
  d,
  alpha = 0.05,
  version,
  gammai,
  w0,
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
  LORDstar to use. `version='async'` requires a column of decision times
  (\`decision.times'). `version='dep'` requires a column of lags
  (\`lags'). `version='batch'` requires a vector of batch sizes
  (\`batch.sizes').

- gammai:

  Optional vector of \\\gamma_i\\. A default is provided as proposed by
  Javanmard and Montanari (2018), equation 31.

- w0:

  Initial \`wealth' of the procedure, defaults to \\\alpha/10\\.

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
and a column describing the conflict sets for the hypotheses. This takes
the form of a vector of decision times or lags. Batch sizes can be
specified as a separate argument (see below).

Zrnic et al. (2021) present explicit three versions of LORDstar:

1\) `version='async'` is for an asynchronous testing process, consisting
of tests that start and finish at (potentially) random times. The
discretised finish times of the test correspond to the decision times.
These decision times are given as the input `decision.times` for this
version of the LORDstar algorithm.

2\) `version='dep'` is for online testing under local dependence of the
p-values. More precisely, for any \\t\>0\\ we allow the p-value \\p_t\\
to have arbitrary dependence on the previous \\L_t\\ p-values. The fixed
sequence \\L_t\\ is referred to as \`lags', and is given as the input
`lags` for this version of the LORDstar algorithm.

3\) `version='batch'` is for controlling the mFDR in mini-batch testing,
where a mini-batch represents a grouping of tests run asynchronously
which result in dependent p-values. Once a mini-batch of tests is fully
completed, a new one can start, testing hypotheses independent of the
previous batch. The batch sizes are given as the input `batch.sizes` for
this version of the LORDstar algorithm.

Given an overall significance level \\\alpha\\, LORDstar depends on
\\w_0\\ (where \\0 \le w_0 \le \alpha\\), which represents the intial
\`wealth' of the procedure. The algorithms also require a sequence of
non-negative non-increasing numbers \\\gamma_i\\ that sum to 1.

Note that these LORDstar algorithms control the *modified* FDR (mFDR).
The \`async' version also controls the usual FDR if the p-values are
assumed to be independent.

Further details of the LORDstar algorithms can be found in Zrnic et al.
(2021).

## References

Javanmard, A. and Montanari, A. (2018) Online Rules for Control of False
Discovery Rate and False Discovery Exceedance. *Annals of Statistics*,
46(2):526-554.

Zrnic, T., Ramdas, A. and Jordan, M.I. (2021). Asynchronous Online
Testing of Multiple Hypotheses. *Journal of Machine Learning Research*
22:1-33.

## See also

[`LORD`](https://dsrobertson.github.io/onlineFDR/reference/LORD.md)
presents versions of LORD for *synchronous* p-values, i.e. where each
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

LORDstar(sample.df, version='async')
#>          pval       alphai R
#> 1  2.9000e-08 2.675839e-04 1
#> 2  6.7430e-02 5.819103e-05 0
#> 3  1.5140e-02 2.457817e-03 0
#> 4  8.1740e-02 5.649373e-04 0
#> 5  1.7100e-03 4.810068e-04 0
#> 6  3.6000e-05 4.011918e-04 1
#> 7  7.9149e-01 3.410964e-04 0
#> 8  2.7201e-01 2.971630e-03 0
#> 9  2.8295e-01 8.426900e-04 0
#> 10 7.5900e-08 7.286513e-04 1
#> 11 6.9274e-01 6.227110e-04 0
#> 12 3.0443e-01 3.217229e-03 0
#> 13 1.3600e-03 1.060551e-03 0
#> 14 7.2342e-01 9.246646e-04 0
#> 15 5.4757e-01 8.010731e-04 0

sample.df2 <- data.frame(
id = c('A15432', 'B90969', 'C18705', 'B49731', 'E99902',
    'C38292', 'A30619', 'D46627', 'E29198', 'A41418',
    'D51456', 'C88669', 'E03673', 'A63155', 'B66033'),
pval = c(2.90e-08, 0.06743, 0.01514, 0.08174, 0.00171,
        3.60e-05, 0.79149, 0.27201, 0.28295, 7.59e-08,
        0.69274, 0.30443, 0.00136, 0.72342, 0.54757),
lags = rep(1,15))

LORDstar(sample.df2, version='dep')
#>          pval lag       alphai R
#> 1  2.9000e-08   1 2.675839e-04 1
#> 2  6.7430e-02   1 5.819103e-05 0
#> 3  1.5140e-02   1 2.457817e-03 0
#> 4  8.1740e-02   1 5.649373e-04 0
#> 5  1.7100e-03   1 4.810068e-04 0
#> 6  3.6000e-05   1 4.011918e-04 1
#> 7  7.9149e-01   1 3.410964e-04 0
#> 8  2.7201e-01   1 2.971630e-03 0
#> 9  2.8295e-01   1 8.426900e-04 0
#> 10 7.5900e-08   1 7.286513e-04 1
#> 11 6.9274e-01   1 6.227110e-04 0
#> 12 3.0443e-01   1 3.217229e-03 0
#> 13 1.3600e-03   1 1.060551e-03 0
#> 14 7.2342e-01   1 9.246646e-04 0
#> 15 5.4757e-01   1 8.010731e-04 0
```
