# ADDIS_exhaustive: Exhaustive ADDIS procedure for online FDR control

Implements an exhaustive variant of the ADDIS algorithm for online FDR
control by adapting code from the Fischer, L.: Exhaustive ADDIS
procedures for online FWER control.

## Usage

``` r
ADDIS_exhaustive(d, alpha = 0.05, tau = 0.5, lambda = 0.25, gamma = NULL)
```

## Arguments

- d:

  Either a vector of p-values, or a dataframe with at least a \`pval\`
  column (and optionally \`id\`).

- alpha:

  Overall significance level of the procedure, default 0.05.

- tau:

  Optional threshold for hypotheses to be selected for testing. Must be
  between 0 and 1, defaults to 0.5.

- lambda:

  Optional parameter that sets the threshold for \`candidate'
  hypotheses. Must be between 0 and tau, defaults to 0.25.

- gamma:

  Optional vector of initial weights. If \`NULL\` (the default), a
  decreasing sequence proportional to j^(-1.6) is used, as in ADDIS().

## Value

A dataframe with the original p-values \`pval\`, the per-hypothesis
testing levels \`alphai\`, and the indicator of discoveries \`R\`.

## References

Fischer, L.: Exhaustive ADDIS procedures for online FWER control.
arXiv:2308.13827 \<https://arxiv.org/abs/2308.13827\>

## Author

Lasse Fischer
