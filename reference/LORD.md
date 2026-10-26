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
  given in Javanmard and Montanari (2018), example 3.7.

- version:

  Takes values '++', 3, 'discard', or 'dep'. This specifies the version
  of LORD to use, and defaults to '++'.

- w0:

  Initial \`wealth' of the procedure, defaults to \\\alpha/10\\
  (\\\min(\alpha/10, \tau\alpha/2)\\ for `version='discard'`). Version
  'discard' requires \\w_0 \le \tau\alpha\\.

- b0:

  The \`payout' for rejecting a hypothesis in versions 3 and 'dep' of
  LORD. Defaults to \\\alpha - w_0\\.

- tau.discard:

  Optional threshold for hypotheses to be selected for testing. Must be
  in (0, 1\], defaults to 0.5. This is only used if `version='discard'`.

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

  A dataframe with the p-values `pval` (reordered if there are batches
  and `random = TRUE`), the LORD-adjusted significance thresholds
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
arriving one at a time.

LORD++ and the discarding version of LORD provably control the FDR for
independent p-values; for LORD 3 this has only been shown empirically
(see below for dependent p-values). Given an overall significance level
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
(0,1\]\\ represents the threshold for a hypothesis to be selected for
testing: p-values greater than \\\tau\\ are implicitly \`discarded' by
the procedure.

Note that FDR control also holds for the LORD procedure if only the
p-values corresponding to true nulls are mutually independent, and
independent from the non-null p-values.

FDR control is not guaranteed for positively dependent p-values in
general, even with the default settings. For example, the FDR of LORD++
can exceed \\\alpha\\ when long runs of positively correlated null
p-values follow strong signals; see
[`vignette("theory")`](https://dsrobertson.github.io/onlineFDR/articles/theory.md).

For dependent p-values, a modified LORD procedure was proposed in
Javanmard and Montanari (2018), which is called by setting
`version='dep'`. Given an overall significance level \\\alpha\\, we
choose a sequence of non-negative numbers \\\xi_i\\ (supplied through
the argument `gammai`) such that \\\sum_i \xi_i (1 + \log i) \le
\alpha/b_0\\ if \\w_0 \le b_0\\, which is the condition in Javanmard and
Montanari (2018), example 3.7. If \\w_0 \> b_0\\, the condition is
instead \\\sum_i (w_0 + b_0 \log i) \xi_i \le \alpha\\, which follows
from their Theorem 3.6.

Further details of the LORD procedures can be found in Javanmard and
Montanari (2018), Ramdas et al. (2017) and Tian and Ramdas (2019).

## References

Javanmard, A. and Montanari, A. (2018). Online Rules for Control of
False Discovery Rate and False Discovery Exceedance. *Annals of
Statistics*, 46(2):526-554.

Ramdas, A., Yang, F., Wainwright, M.J. and Jordan, M.I. (2017). Online
control of the false discovery rate with decaying memory. *Advances in
Neural Information Processing Systems 30*, 5650-5659.

Tian, J. and Ramdas, A. (2019). ADDIS: an adaptive discarding algorithm
for online FDR control with conservative nulls. *Advances in Neural
Information Processing Systems 32*, 9388-9396.

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
#> 1  2.9000e-08 0.0002740772 1 A15432
#> 2  6.7430e-02 0.0025262981 0 B90969
#> 3  1.5140e-02 0.0005871934 0 C18705
#> 4  8.1740e-02 0.0004991051 0 B49731
#> 5  1.7100e-03 0.0004157566 0 E99902
#> 6  3.6000e-05 0.0003530941 1 C38292
#> 7  7.9149e-01 0.0030466825 0 A30619
#> 8  2.7201e-01 0.0008655170 0 D46627
#> 9  2.8295e-01 0.0007482932 0 E29198
#> 10 7.5900e-08 0.0006394643 1 A41418
#> 11 6.9274e-01 0.0032966958 0 D51456
#> 12 3.0443e-01 0.0010874871 0 C88669
#> 13 1.3600e-03 0.0009481457 0 E03673
#> 14 7.2342e-01 0.0008214267 0 A63155
#> 15 5.4757e-01 0.0007230924 0 B66033

set.seed(1); LORD(sample.df, version='dep')
#>          pval       alphai R     id
#> 1  2.9000e-08 2.091536e-03 1 A15432
#> 2  6.7430e-02 1.002023e-02 0 B90969
#> 3  1.5140e-02 1.677758e-03 0 C18705
#> 4  8.1740e-02 6.262642e-04 0 B49731
#> 5  1.7100e-03 3.201778e-04 0 E99902
#> 6  2.7201e-01 1.933720e-04 0 D46627
#> 7  3.6000e-05 1.293950e-04 1 C38292
#> 8  7.9149e-01 1.548148e-04 0 A30619
#> 9  7.5900e-08 1.166479e-04 1 A41418
#> 10 2.8295e-01 1.422610e-04 0 E29198
#> 11 6.9274e-01 1.145116e-04 0 D51456
#> 12 7.2342e-01 9.432385e-05 0 A63155
#> 13 3.0443e-01 7.916865e-05 0 C88669
#> 14 5.4757e-01 6.749296e-05 0 B66033
#> 15 1.3600e-03 5.830041e-05 0 E03673

set.seed(1); LORD(sample.df, version='discard')
#>          pval       alphai R     id
#> 1  2.9000e-08 0.0002740772 1 A15432
#> 2  6.7430e-02 0.0011559120 0 B90969
#> 3  1.5140e-02 0.0002891777 0 C18705
#> 4  8.1740e-02 0.0002452791 0 B49731
#> 5  1.7100e-03 0.0002046653 0 E99902
#> 6  2.7201e-01 0.0001741324 0 D46627
#> 7  3.6000e-05 0.0001510949 1 C38292
#> 8  7.9149e-01 0.0015036586 0 A30619
#> 9  7.5900e-08 0.0015036586 1 A41418
#> 10 2.8295e-01 0.0017875335 0 E29198
#> 11 6.9274e-01 0.0006595026 0 D51456
#> 12 7.2342e-01 0.0006595026 0 A63155
#> 13 3.0443e-01 0.0006595026 0 C88669
#> 14 5.4757e-01 0.0005630964 0 B66033
#> 15 1.3600e-03 0.0005630964 0 E03673

set.seed(1); LORD(sample.df, alpha=0.1, w0=0.05)
#>          pval       alphai R     id
#> 1  2.9000e-08 0.0027407722 1 A15432
#> 2  6.7430e-02 0.0033368035 0 B90969
#> 3  1.5140e-02 0.0011036834 0 C18705
#> 4  8.1740e-02 0.0009298346 0 B49731
#> 5  1.7100e-03 0.0007801059 0 E99902
#> 6  2.7201e-01 0.0006675540 0 D46627
#> 7  3.6000e-05 0.0005820572 1 C38292
#> 8  7.9149e-01 0.0059969897 0 A30619
#> 9  7.5900e-08 0.0016543231 1 A41418
#> 10 2.8295e-01 0.0069157328 0 E29198
#> 11 6.9274e-01 0.0024192906 0 D51456
#> 12 7.2342e-01 0.0020836384 0 A63155
#> 13 3.0443e-01 0.0017901616 0 C88669
#> 14 5.4757e-01 0.0015648092 0 B66033
#> 15 1.3600e-03 0.0013898398 1 E03673

```
