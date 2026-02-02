# LORD (dep): Online FDR control based on recent discovery for dependent p-values

This funcion is deprecated, please use
[`LORD`](https://dsrobertson.github.io/onlineFDR/reference/LORD.md)
instead with `version = 'dep'`.

## Usage

``` r
LORDdep(
  d,
  alpha = 0.05,
  xi,
  w0 = alpha/10,
  b0 = alpha - w0,
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

- xi:

  Optional vector of \\\xi_i\\. A default is provided to satisfy the
  condition given in Javanmard and Montanari (2018), example 3.7.

- w0:

  Initial \`wealth' of the procedure. Defaults to \\\alpha/10\\.

- b0:

  The \`payout' for rejecting a hypothesis. Defaults to \\\alpha -
  w_0\\.

- random:

  Logical. If `TRUE` (the default), then the order of the p-values in
  each batch (i.e. those that have exactly the same date) is randomised.

- date.format:

  Optional string giving the format that is used for dates.

## Value

- d.out:

  A dataframe with the original data `d` (which will be reordered if
  there are batches and `random = TRUE`), the LORD-adjusted significance
  thresholds \\\alpha_i\\ and the indicator function of discoveries `R`.
  Hypothesis \\i\\ is rejected if the \\i\\-th p-value is less than or
  equal to \\\alpha_i\\, in which case `R[i] = 1` (otherwise
  `R[i] = 0`).

## Details

LORDdep implements the LORD procedure for online FDR control for
dependent p-values, where LORD stands for (significance) Levels based On
Recent Discovery, as presented by Javanmard and Montanari (2018).

The function takes as its input either a vector of p-values or a
dataframe with three columns: an identifier (\`id'), date (\`date') and
p-value (\`pval'). The case where p-values arrive in batches corresponds
to multiple instances of the same date. If no column of dates is
provided, then the p-values are treated as being ordered in sequence,
arriving one at a time.

This modified LORD procedure controls FDR for dependent p-values. Given
an overall significance level \\\alpha\\, we choose a sequence of
non-negative numbers \\\xi_i\\ such that they satisfy a condition given
in Javanmard and Montanari (2018), example 3.8.

The procedure depends on constants \\w_0\\ and \\b_0\\, where \\w_0 \ge
0\\ represents the intial \`wealth' and \\b_0 \> 0\\ represents the
\`payout' for rejecting a hypothesis. We require \\w_0+b_0 \le \alpha\\
for FDR control to hold.

Further details of the modified LORD procedure can be found in Javanmard
and Montanari (2018).

## References

Javanmard, A. and Montanari, A. (2018) Online Rules for Control of False
Discovery Rate and False Discovery Exceedance. *Annals of Statistics*,
46(2):526-554.
