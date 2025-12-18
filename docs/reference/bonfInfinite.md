# Online FDR control based on a Bonferroni-like test

This funcion is deprecated, please use
[`Alpha_spending`](https://dsrobertson.github.io/onlineFDR/reference/Alpha_spending.md)
instead.

## Usage

``` r
bonfInfinite(
  d,
  alpha = 0.05,
  alphai,
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

- alphai:

  Optional vector of \\\alpha_i\\, where hypothesis \\i\\ is rejected if
  the \\i\\-th p-value is less than or equal to \\\alpha_i\\. A default
  is provided as proposed by Javanmard and Montanari (2018), equation
  31.

- random:

  Logical. If `TRUE` (the default), then the order of the p-values in
  each batch (i.e. those that have exactly the same date) is randomised.

- date.format:

  Optional string giving the format that is used for dates.

## Value

- d.out:

  A dataframe with the original data `d` (which will be reordered if
  there are batches and `random = TRUE`), the adjusted signifcance
  thresholds `alphai` and the indicator function of discoveries `R`,
  where `R[i] = 1` corresponds to hypothesis \\i\\ being rejected
  (otherwise `R[i] = 0`).

## Details

Implements online FDR control using a Bonferroni-like test.

The function takes as its input either a vector of p-values, or a
dataframe with three columns: an identifier (\`id'), date (\`date') and
p-value (\`pval'). The case where p-values arrive in batches corresponds
to multiple instances of the same date. If no column of dates is
provided, then the p-values are treated as being ordered in sequence,
arriving one at a time.

The procedure controls FDR for a potentially infinite stream of p-values
by using a Bonferroni-like test. Given an overall significance level
\\\alpha\\, we choose a (potentially infinite) sequence of non-negative
numbers \\\alpha_i\\ such that they sum to \\\alpha\\. Hypothesis \\i\\
is rejected if the \\i\\-th p-value is less than or equal to
\\\alpha_i\\.

## References

Javanmard, A. and Montanari, A. (2018) Online Rules for Control of False
Discovery Rate and False Discovery Exceedance. *Annals of Statistics*,
46(2):526-554.
