# ADDIS-exhaustive: Exhaustive ADDIS-spending procedure for online FWER control

Implements an exhaustive variant of the ADDIS-spending algorithm for
online FWER control, as presented by Fischer et al. (2023). The
procedure is a uniform improvement of ADDIS-spending, and no other FWER
controlling procedure can enlarge the event of rejecting any hypothesis.

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

## Details

The function takes as its input either a vector of p-values, or a
dataframe with two columns: an identifier (\`id') and p-value (\`pval').
Given an overall significance level \\\alpha\\, ADDIS-exhaustive depends
on constants \\\lambda\\ and \\\tau\\, where \\\lambda \< \tau\\. Here
\\\tau \in (0,1)\\ represents the threshold for a hypothesis to be
selected for testing: p-values greater than \\\tau\\ are implicitly
\`discarded' by the procedure, while \\\lambda \in (0,1)\\ sets the
threshold for a p-value to be a candidate for rejection:
ADDIS-exhaustive will never reject a p-value larger than \\\lambda\\.
The algorithms also require a sequence of non-negative non-increasing
numbers \\\gamma_i\\ that sum to 1.

The ADDIS-exhaustive procedure provably controls the FWER in the strong
sense for independent p-values.

## References

Fischer, L., Bofill Roig, M. and Brannath W. (2024). An exhaustive ADDIS
principle for online FWER control. *Biometrical Journal* 66(3) 2300237.

## Author

Lasse Fischer
