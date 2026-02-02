# LORD: Online FDR control based on recent discovery

Implements the LORD procedure for online FDR control, where LORD stands
for (significance) Levels based On Recent Discovery, as presented by
Javanmard and Montanari (2018) and Ramdas et al. (2017).

## Usage

``` r
LORD(
  d,
  alpha = 0.05,
  gammai,
  version = "++",
  w0,
  b0,
  tau.discard = 0.5,
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

  Overall significance level of the FDR procedure, the default is 0.05.

- gammai:

  Optional vector of \\\gamma_i\\. A default is provided as proposed by
  Javanmard and Montanari (2018), equation 31 for all versions of LORD
  except 'dep'. The latter is provided a default to satisfy a condition
  given in Javanmard and Montanari (2018), example 3.8.

- version:

  Takes values '++', 3, 'discard', or 'dep'. This specifies the version
  of LORD to use, and defaults to '++'.

- w0:

  Initial \`wealth' of the procedure, defaults to \\\alpha/10\\.

- b0:

  The 'payout' for rejecting a hypothesis in all versions of LORD except
  for '++'. Defaults to \\\alpha - w_0\\.

- tau.discard:

  Optional threshold for hypotheses to be selected for testing. Must be
  between 0 and 1, defaults to 0.5. This is required if
  `version='discard'`.

- random:

  Logical. If `TRUE` (the default), then the order of the p-values in
  each batch (i.e. those that have exactly the same date) is randomised.

- display_progress:

  Logical. If `TRUE` prints out a progress bar for the algorithm
  runtime.

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

The function takes as its input either a vector of p-values or a
dataframe with three columns: an identifier (\`id'), date (\`date') and
p-value (\`pval'). The case where p-values arrive in batches corresponds
to multiple instances of the same date. If no column of dates is
provided, then the p-values are treated as being ordered in sequence,
arriving one at a time..

The LORD procedure provably controls FDR for independent p-values (see
below for dependent p-values). Given an overall significance level
\\\alpha\\, we choose a sequence of non-negative non-increasing numbers
\\\gamma_i\\ that sum to 1.

Javanmard and Montanari (2018) presented versions of LORD which differ
in the way the adjusted significance thresholds \\\alpha_i\\ are
calculated. The significance thresholds for LORD 2 are based on all
previous discovery times. LORD 2 has been superseded by the algorithm
given in Ramdas et al. (2017), LORD++ (`version='++'`), which is the
default version. The significance thresholds for LORD 3 (`version=3`)
are based on the time of the last discovery as well as the 'wealth'
accumulated at that time. Finally, Tian and Ramdas (2019) presented a
version of LORD (`version='discard'`) that can improve the power of the
procedure in the presence of conservative nulls by adaptively
\`discarding' these p-values.

LORD depends on constants \\w_0\\ and (for versions 3 and 'dep')
\\b_0\\, where \\0 \le w_0 \le \alpha\\ represents the initial \`wealth'
of the procedure and \\b_0 \> 0\\ represents the \`payout' for rejecting
a hypothesis. We require \\w_0+b_0 \le \alpha\\ for FDR control to hold.
Version 'discard' also depends on a constant \\\tau\\, where \\\tau \in
(0,1)\\ represents the threshold for a hypothesis to be selected for
testing: p-values greater than \\\tau\\ are implicitly \`discarded' by
the procedure.

Note that FDR control also holds for the LORD procedure if only the
p-values corresponding to true nulls are mutually independent, and
independent from the non-null p-values.

For dependent p-values, a modified LORD procedure was proposed in
Javanmard and Montanari (2018), which is called be setting
`version='dep'`. Given an overall significance level \\\alpha\\, we
choose a sequence of non-negative numbers \\\xi_i\\ such that they
satisfy a condition given in Javanmard and Montanari (2018), example
3.8.

Further details of the LORD procedures can be found in Javanmard and
Montanari (2018), Ramdas et al. (2017) and Tian and Ramdas (2019).

## References

Javanmard, A. and Montanari, A. (2018) Online Rules for Control of False
Discovery Rate and False Discovery Exceedance. *Annals of Statistics*,
46(2):526-554.

Ramdas, A., Yang, F., Wainwright M.J. and Jordan, M.I. (2017). Online
control of the false discovery rate with decaying memory. *Advances in
Neural Information Processing Systems 30*, 5650-5659.

Tian, J. and Ramdas, A. (2019). ADDIS: an adaptive discarding algorithm
for online FDR control with conservative nulls. *Advances in Neural
Information Processing Systems*, 9388-9396.

## See also

[`LORDstar`](https://dsrobertson.github.io/onlineFDR/reference/LORDstar.md)
presents versions of LORD for *asynchronous* testing, i.e. where each
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

LORD(sample.df, random=FALSE)
#>          pval       alphai R     id
#> 1  2.9000e-08 0.0002675839 1 A15432
#> 2  6.7430e-02 0.0024664457 0 B90969
#> 3  1.5140e-02 0.0005732818 0 C18705
#> 4  8.1740e-02 0.0004872805 0 B49731
#> 5  1.7100e-03 0.0004059066 0 E99902
#> 6  3.6000e-05 0.0003447286 1 C38292
#> 7  7.9149e-01 0.0029745013 0 A30619
#> 8  2.7201e-01 0.0008450114 0 D46627
#> 9  2.8295e-01 0.0007305648 0 E29198
#> 10 7.5900e-08 0.0006243143 1 A41418
#> 11 6.9274e-01 0.0032185913 0 D51456
#> 12 3.0443e-01 0.0010617227 0 C88669
#> 13 1.3600e-03 0.0009256825 0 E03673
#> 14 7.2342e-01 0.0008019657 0 A63155
#> 15 5.4757e-01 0.0007059610 0 B66033

set.seed(1); LORD(sample.df, version='dep')
#>          pval       alphai R     id
#> 1  2.9000e-08 2.091542e-03 1 A15432
#> 2  6.7430e-02 1.002025e-02 0 B90969
#> 3  1.5140e-02 1.677763e-03 0 C18705
#> 4  8.1740e-02 6.262659e-04 0 B49731
#> 5  1.7100e-03 3.201787e-04 0 E99902
#> 6  2.7201e-01 1.933725e-04 0 D46627
#> 7  3.6000e-05 1.293954e-04 1 C38292
#> 8  7.9149e-01 1.548152e-04 0 A30619
#> 9  7.5900e-08 1.166482e-04 1 A41418
#> 10 2.8295e-01 1.422614e-04 0 E29198
#> 11 6.9274e-01 1.145119e-04 0 D51456
#> 12 7.2342e-01 9.432408e-05 0 A63155
#> 13 3.0443e-01 7.916885e-05 0 C88669
#> 14 5.4757e-01 6.749313e-05 0 B66033
#> 15 1.3600e-03 5.830055e-05 0 E03673

set.seed(1); LORD(sample.df, version='discard')
#>          pval       alphai R     id
#> 1  2.9000e-08 0.0002675839 1 A15432
#> 2  6.7430e-02 0.0011285264 0 B90969
#> 3  1.5140e-02 0.0002823266 0 C18705
#> 4  8.1740e-02 0.0002394680 0 B49731
#> 5  1.7100e-03 0.0001998165 0 E99902
#> 6  2.7201e-01 0.0001700069 0 D46627
#> 7  3.6000e-05 0.0001475152 1 C38292
#> 8  7.9149e-01 0.0014680343 0 A30619
#> 9  7.5900e-08 0.0014680343 1 A41418
#> 10 2.8295e-01 0.0017451837 0 E29198
#> 11 6.9274e-01 0.0006438778 0 D51456
#> 12 7.2342e-01 0.0006438778 0 A63155
#> 13 3.0443e-01 0.0006438778 0 C88669
#> 14 5.4757e-01 0.0005497556 0 B66033
#> 15 1.3600e-03 0.0005497556 0 E03673

set.seed(1); LORD(sample.df, alpha=0.1, w0=0.05)
#>          pval       alphai R     id
#> 1  2.9000e-08 0.0026758385 1 A15432
#> 2  6.7430e-02 0.0032577488 0 B90969
#> 3  1.5140e-02 0.0010775352 0 C18705
#> 4  8.1740e-02 0.0009078052 0 B49731
#> 5  1.7100e-03 0.0007616238 0 E99902
#> 6  2.7201e-01 0.0006517385 0 D46627
#> 7  3.6000e-05 0.0005682672 1 C38292
#> 8  7.9149e-01 0.0058549106 0 A30619
#> 9  7.5900e-08 0.0016151293 1 A41418
#> 10 2.8295e-01 0.0067518870 0 E29198
#> 11 6.9274e-01 0.0023619734 0 D51456
#> 12 7.2342e-01 0.0020342733 0 A63155
#> 13 3.0443e-01 0.0017477495 0 C88669
#> 14 5.4757e-01 0.0015277362 0 B66033
#> 15 1.3600e-03 0.0013569121 0 E03673

```
