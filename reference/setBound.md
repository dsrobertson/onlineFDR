# setBound

Calculates a default sequence of non-negative numbers \\\gamma_i\\ that
sum to 1 (for LOND and LONDstar, a sequence \\\beta_i\\ that sums to
\\\alpha\\), given an upper bound \\N\\ on the number of hypotheses to
be tested. For LORDdep the sequence \\\xi_i\\ instead satisfies
\\\sum\_{i=1}^N \xi_i (1 + \log i) = 1\\, the condition in Javanmard and
Montanari (2018), example 3.7, for `LORD(version='dep')`.

## Usage

``` r
setBound(alg, alpha = 0.05, N)
```

## Arguments

- alg:

  A string that takes the value of one of the following: LOND, LORD,
  LORDdep, SAFFRON, ADDIS, LONDstar, LORDstar, SAFFRONstar,
  Alpha_investing, Alpha_spending, ADDIS_spending or online_fallback

- alpha:

  Overall significance level of the FDR procedure, the default is 0.05.
  The bounds for LOND and LONDstar depend on alpha.

- N:

  An upper bound on the number of hypotheses to be tested

## Value

- bound:

  A vector giving the values of a default sequence \\\gamma_i\\ of
  non-negative numbers.
