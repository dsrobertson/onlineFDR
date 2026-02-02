# setBound

Calculates a default sequence of non-negative numbers \\\gamma_i\\ that
sum to 1, given an upper bound \\N\\ on the number of hypotheses to be
tested.

## Usage

``` r
setBound(alg, alpha = 0.05, N)
```

## Arguments

- alg:

  A string that takes the value of one of the following: LOND, LORD,
  LORDdep, SAFFRON, ADDIS, LONDstar, LORDstar, SAFFRONstar, or
  Alpha_investing

- alpha:

  Overall significance level of the FDR procedure, the default is 0.05.
  The bounds for LOND and LORDdep depend on alpha.

- N:

  An upper bound on the number of hypotheses to be tested

## Value

- bound:

  A vector giving the values of a default sequence \\\gamma_i\\ of
  nonnegative numbers.
