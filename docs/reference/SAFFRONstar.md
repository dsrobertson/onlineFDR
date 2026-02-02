# SAFFRONstar: Adaptive online mFDR control for asynchronous testing

Implements the SAFFRON algorithm for asynchronous online testing, as
presented by Zrnic et al. (2021).

## Usage

``` r
SAFFRONstar(
  d,
  alpha = 0.05,
  version,
  gammai,
  w0,
  lambda = 0.5,
  batch.sizes,
  display_progress = FALSE
)
```

## Arguments

- d:

  Either a vector of p-values, or a dataframe with three columns: an
  identifier (\`id'), p-value (\`pval'), and either decision.times', or
  \`lags', depending on which version you're using. See version for more
  details.

- alpha:

  Overall significance level of the procedure, the default is 0.05.

- version:

  Takes values 'async', 'dep' or 'batch'. This specifies the version of
  SAFFRONstar to use. `version='async'` requires a column of decision
  times (\`decision.times'). `version='dep'` requires a column of lags
  (\`lags'). `version='batch'` requires a vector of batch sizes
  (\`batch.sizes').

- gammai:

  Optional vector of \\\gamma_i\\. A default is provided with
  \\\gamma_j\\ proportional to \\1/j^(1.6)\\.

- w0:

  Initial \`wealth' of the procedure, defaults to \\\alpha/10\\.

- lambda:

  Optional threshold for a \`candidate' hypothesis, must be between 0
  and 1. Defaults to 0.5.

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

Zrnic et al. (2021) present explicit three versions of SAFFRONstar:

1\) `version='async'` is for an asynchronous testing process, consisting
of tests that start and finish at (potentially) random times. The
discretised finish times of the test correspond to the decision times.
These decision times are given as the input `decision.times` for this
version of the SAFFRONstar algorithm. For this version of SAFFRONstar,
Tian and Ramdas (2019) presented an algorithm that can improve the power
of the procedure in the presence of conservative nulls by adaptively
\`discarding' these p-values. This can be called by setting the option
`discard=TRUE`.

2\) `version='dep'` is for online testing under local dependence of the
p-values. More precisely, for any \\t\>0\\ we allow the p-value \\p_t\\
to have arbitrary dependence on the previous \\L_t\\ p-values. The fixed
sequence \\L_t\\ is referred to as \`lags', and is given as the input
`lags` for this version of the SAFFRONstar algorithm.

3\) `version='batch'` is for controlling the mFDR in mini-batch testing,
where a mini-batch represents a grouping of tests run asynchronously
which result in dependent p-values. Once a mini-batch of tests is fully
completed, a new one can start, testing hypotheses independent of the
previous batch. The batch sizes are given as the input `batch.sizes` for
this version of the SAFFRONstar algorithm.

Given an overall significance level \\\alpha\\, SAFFRONstar depends on
constants \\w_0\\ and \\\lambda\\, where \\w_0\\ satisfies \\0 \le w_0
\le \alpha\\ and represents the intial \`wealth' of the procedure, and
\\0 \< \lambda \< 1\\ represents the threshold for a \`candidate'
hypothesis. A \`candidate' refers to p-values smaller than \\\lambda\\,
since SAFFRONstar will never reject a p-value larger than \\\lambda\\.
The algorithms also require a sequence of non-negative non-increasing
numbers \\\gamma_i\\ that sum to 1.

Note that these SAFFRONstar algorithms control the *modified* FDR
(mFDR). The \`async' version also controls the usual FDR if the p-values
are assumed to be independent.

Further details of the SAFFRONstar algorithms can be found in Zrnic et
al. (2021).

## References

Zrnic, T., Ramdas, A. and Jordan, M.I. (2021). Asynchronous Online
Testing of Multiple Hypotheses. *Journal of Machine Learning Research*,
22:1-33.

## See also

[`SAFFRON`](https://dsrobertson.github.io/onlineFDR/reference/SAFFRON.md)
presents versions of SAFFRON for *synchronous* p-values, i.e. where each
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

SAFFRONstar(sample.df, version='async')
#>          pval      alphai R
#> 1  2.9000e-08 0.005468627 1
#> 2  6.7430e-02 0.001803974 0
#> 3  1.5140e-02 0.007272601 0
#> 4  8.1740e-02 0.003607948 0
#> 5  1.7100e-03 0.003607948 1
#> 6  3.6000e-05 0.003607948 1
#> 7  7.9149e-01 0.014545202 0
#> 8  2.7201e-01 0.018153151 0
#> 9  2.8295e-01 0.007379710 0
#> 10 7.5900e-08 0.007379710 1
#> 11 6.9274e-01 0.007379710 0
#> 12 3.0443e-01 0.018316965 0
#> 13 1.3600e-03 0.007874188 1
#> 14 7.2342e-01 0.007874188 0
#> 15 5.4757e-01 0.018811442 0

sample.df2 <- data.frame(
id = c('A15432', 'B90969', 'C18705', 'B49731', 'E99902',
    'C38292', 'A30619', 'D46627', 'E29198', 'A41418',
    'D51456', 'C88669', 'E03673', 'A63155', 'B66033'),
pval = c(2.90e-08, 0.06743, 0.01514, 0.08174, 0.00171,
        3.60e-05, 0.79149, 0.27201, 0.28295, 7.59e-08,
        0.69274, 0.30443, 0.00136, 0.72342, 0.54757),
lags = rep(1,15))

SAFFRONstar(sample.df2, version='dep')
#>          pval lag      alphai R
#> 1  2.9000e-08   1 0.005468627 1
#> 2  6.7430e-02   1 0.001803974 0
#> 3  1.5140e-02   1 0.007272601 0
#> 4  8.1740e-02   1 0.003607948 0
#> 5  1.7100e-03   1 0.003607948 1
#> 6  3.6000e-05   1 0.003607948 1
#> 7  7.9149e-01   1 0.014545202 0
#> 8  2.7201e-01   1 0.018153151 0
#> 9  2.8295e-01   1 0.007379710 0
#> 10 7.5900e-08   1 0.007379710 1
#> 11 6.9274e-01   1 0.007379710 0
#> 12 3.0443e-01   1 0.018316965 0
#> 13 1.3600e-03   1 0.007874188 1
#> 14 7.2342e-01   1 0.007874188 0
#> 15 5.4757e-01   1 0.018811442 0
```
