# ADDIS-spending: Adaptive discarding algorithm for online FWER control

Implements the ADDIS algorithm for online FWER control, where ADDIS
stands for an ADaptive algorithm that DIScards conservative nulls, as
presented by Tian and Ramdas (2021). The procedure compensates for the
power loss of Alpha-spending, by including both adaptivity in the
fraction of null hypotheses and the conservativeness of nulls.

## Usage

``` r
ADDIS_spending(
  d,
  alpha = 0.05,
  gammai,
  lambda = 0.25,
  tau = 0.5,
  dep = FALSE,
  display_progress = FALSE
)
```

## Arguments

- d:

  Either a vector of p-values, or a dataframe with three columns: an
  identifier (\`id'), p-value (\`pval'), and lags (\`lags').

- alpha:

  Overall significance level of the procedure, the default is 0.05.

- gammai:

  Optional vector of \\\gamma_i\\. A default is provided with
  \\\gamma_j\\ proportional to \\1/j^(1.6)\\.

- lambda:

  Optional parameter that sets the threshold for \`candidate'
  hypotheses. Must be between 0 and 1, defaults to 0.25.

- tau:

  Optional threshold for hypotheses to be selected for testing. Must be
  between 0 and 1, defaults to 0.5.

- dep:

  Logical. If `TRUE` runs the version for locally dependent p-values

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
and lags, if the dependent version is specified (see below). Given an
overall significance level \\\alpha\\, ADDIS depends on constants
\\\lambda\\ and \\\tau\\, where \\\lambda \< \tau\\. Here \\\tau \in
(0,1)\\ represents the threshold for a hypothesis to be selected for
testing: p-values greater than \\\tau\\ are implicitly \`discarded' by
the procedure, while \\\lambda \in (0,1)\\ sets the threshold for a
p-value to be a candidate for rejection: ADDIS-spending will never
reject a p-value larger than \\\lambda\\. The algorithms also require a
sequence of non-negative non-increasing numbers \\\gamma_i\\ that sum to
1.

The ADDIS-spending procedure provably controls the FWER in the strong
sense for independent p-values. Note that the procedure also controls
the generalised familywise error rate (k-FWER) for \\k \> 1\\ if
\\\alpha\\ is replaced by min(\\1, k\alpha\\).

Tian and Ramdas (2021) also presented a version for handling local
dependence. More precisely, for any \\t\>0\\ we allow the p-value
\\p_t\\ to have arbitrary dependence on the previous \\L_t\\ p-values.
The fixed sequence \\L_t\\ is referred to as \`lags', and is given as
the input `lags` for this version of the ADDIS-spending algorithm.

Further details of the ADDIS-spending algorithms can be found in Tian
and Ramdas (2021).

## References

Tian, J. and Ramdas, A. (2021). Online control of the familywise error
rate. *Statistical Methods for Medical Research* 30(4):976–993.

## See also

[`ADDIS`](https://dsrobertson.github.io/onlineFDR/reference/ADDIS.md)
provides online control of the FDR.

## Examples

``` r
sample.df <- data.frame(
id = c('A15432', 'B90969', 'C18705', 'B49731', 'E99902',
    'C38292', 'A30619', 'D46627', 'E29198', 'A41418',
    'D51456', 'C88669', 'E03673', 'A63155', 'B66033'),
pval = c(2.90e-08, 0.06743, 0.01514, 0.08174, 0.00171,
        3.60e-05, 0.79149, 0.27201, 0.28295, 7.59e-08,
        0.69274, 0.30443, 0.00136, 0.72342, 0.54757),
lags = rep(1,15))

ADDIS_spending(sample.df) #independent
#>          pval       alphai R
#> 1  2.9000e-08 0.0054686271 1
#> 2  6.7430e-02 0.0054686271 0
#> 3  1.5140e-02 0.0054686271 0
#> 4  8.1740e-02 0.0054686271 0
#> 5  1.7100e-03 0.0054686271 1
#> 6  3.6000e-05 0.0054686271 1
#> 7  7.9149e-01 0.0054686271 0
#> 8  2.7201e-01 0.0054686271 0
#> 9  2.8295e-01 0.0018039742 0
#> 10 7.5900e-08 0.0009429405 1
#> 11 6.9274e-01 0.0009429405 0
#> 12 3.0443e-01 0.0009429405 0
#> 13 1.3600e-03 0.0005950895 0
#> 14 7.2342e-01 0.0005950895 0
#> 15 5.4757e-01 0.0005950895 0

ADDIS_spending(sample.df, dep = TRUE) #Locally dependent
#>          pval       alphai R
#> 1  2.9000e-08 0.0054686271 1
#> 2  6.7430e-02 0.0018039742 0
#> 3  1.5140e-02 0.0018039742 0
#> 4  8.1740e-02 0.0018039742 0
#> 5  1.7100e-03 0.0018039742 1
#> 6  3.6000e-05 0.0018039742 1
#> 7  7.9149e-01 0.0018039742 0
#> 8  2.7201e-01 0.0018039742 0
#> 9  2.8295e-01 0.0009429405 0
#> 10 7.5900e-08 0.0005950895 1
#> 11 6.9274e-01 0.0005950895 0
#> 12 3.0443e-01 0.0005950895 0
#> 13 1.3600e-03 0.0004164149 0
#> 14 7.2342e-01 0.0004164149 0
#> 15 5.4757e-01 0.0004164149 0
```
