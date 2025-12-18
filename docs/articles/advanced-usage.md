# Advanced usage of onlineFDR

## Brief Background of the `onlineFDR` algorithms

Javanmard and Montanari proposed two procedures, LOND and LORD, to
control the FDR in an online manner (Javanmard and Montanari (2015,
2018)), with the latter extended by Ramdas *et al.* (2017). The LOND
procedure sets the adjusted significance thresholds based on the number
of discoveries made so far, while LORD sets them according to the time
of the most recent discovery. Ramdas *et al.* (2018) then proposed the
SAFFRON procedure, which provides an adaptive method of online FDR
control. They also proposed a variant of the Alpha-investing algorithm
of Foster and Stine (2008) that guarantees FDR control, using SAFFRON’s
update rule.

Subsequently, Zrnic *et al.* (2021) proposed procedures to control the
modified FDR (mFDR) in the context of *asynchronous* testing, i.e. where
each hypothesis test can itself be a sequential process and the tests
can overlap in time. They presented asynchronous versions of the LOND,
LORD and SAFFRON procedures for a variety of trial settings. For both
synchronous and asynchronous testing, Tian & Ramdas (2019) proposed the
ADDIS algorithms which compensate for the loss in power in the presence
of conservative nulls by adaptively ‘discarding’ these p-values.

Finally, Tian & Ramdas (2021) proposed procedures that provide online
control of the FWER. One procedure, online fallback, gives a uniform
improvement to the naive Alpha-spending procedure (see below). The
ADDIS-spending procedure compensates for the power loss of these
procedures by including both adapativity in the fraction of null
hypotheses and the conservativeness of nulls.

## Variations to the default options

In the following section, we consider the arguments that a typical user
might consider amending for their analysis.

### Common arguments

As a default, the `alpha` argument is set to 0.05, where `alpha` sets
the overall significance level of the FDR of FWER controlling procedure.
By convention, the standard significance level utilised is the 5%.
However, there are applications where an alternate threshold could be
considered. For example, a more stringent threshold might be appropriate
when there are limited resources to follow up significant findings. A
less stringent threshold might be appropriate when the downstream
analysis is a global analysis which can tolerate a higher proportion of
false positives.

To ensure correct interpretation of the dates provided there is a
date.format argument. As a default, the date format is set to receive
dates as year-month(00-12)-day(number). The following website provides
clear guidance on symbols used to interpret the date information:
<https://stat.ethz.ch/R-manual/R-devel/library/base/html/strptime.html>

As a default, the `random` argument is set to `TRUE`. In this situation,
the order of p-values in each batch (i.e. with the same date) are
randomised. This is to avoid the risk of p-values being ordered
post-hoc, which can lead to an inflation of the FDR. As the dataset
grows the data is reprocessed. To ensure the consistency of the output
(with the randomisation within the previous batches remaining the same),
it is necessary to set the same `seed` for all analyses.

The user also has the option to turn off the randomisation step, by
setting the `random` argument to `FALSE`. This approach would be
appropriate if the user has both a date *and* a time stamp for the
p-values, in which case the data should be ordered by date and time
beforehand and then passed to a wrapper function. Another scenario would
be when p-values within the batches are ordered using *independent* side
information, so that hypotheses most likely to be rejected come first,
which would potentially increase the power of the procedure (see
Javanmard and Montanari (2018) and Li and Barber (2017)).

### LOND

As a default, the `dep` argument is set to `FALSE`. Alternatively, this
can be set to `TRUE` and will implement the LOND procedure to guarantee
FDR control for arbitrarily dependent p-values. This method will in
general be more conservative.

``` r
set.seed(1); results.indep <- LOND(sample.df)    # for independent p-values
set.seed(1); results.dep <- LOND(sample.df, dep=TRUE)   # for dependent p-values

# compare adjusted significance thresholds
cbind(independent = results.indep$alphai, dependent = results.dep$alphai)
#>        independent    dependent
#>  [1,] 0.0026758385 0.0026758385
#>  [2,] 0.0011638206 0.0007758804
#>  [3,] 0.0009912499 0.0005406818
#>  [4,] 0.0008243606 0.0003956931
#>  [5,] 0.0006988870 0.0003060819
#>  [6,] 0.0006045900 0.0002467714
#>  [7,] 0.0005319444 0.0002051576
#>  [8,] 0.0007117838 0.0002618915
#>  [9,] 0.0006421423 0.0002269882
#> [10,] 0.0007796504 0.0002661860
#> [11,] 0.0007155186 0.0002369363
#> [12,] 0.0006610273 0.0002130140
#> [13,] 0.0006141682 0.0001931265
#> [14,] 0.0005734509 0.0001763616
#> [15,] 0.0005377472 0.0001620585
```

The vector `betai` is supplied by default, but can optionally be
specified by the user (as described above, see the formula for
$\beta_{j}$[here](#LOND_beta)).

### LORD

The default version of LORD used is version ‘++’, but the user can
optionally specify versions 3, ‘discard’ and ‘dep’ using the `version`
argument (see [here](#LORD) for further details about the different
versions).

``` r
set.seed(1); results.LORD.plus <- LORD(sample.df)
set.seed(1); results.LORD3 <- LORD(sample.df, version=3)
set.seed(1); results.LORD.discard <- LORD(sample.df, version='discard')
set.seed(1); results.LORD.dep <- LORD(sample.df, version='dep') 

# compare adjusted significance thresholds
cbind(LORD.plus = results.LORD.plus$alphai,
    LORD3 = results.LORD3$alphai,
    LORD.discard  = results.LORD.discard$alphai,
    LORD.dep = results.LORD.dep$alphai)
#>          LORD.plus        LORD3 LORD.discard     LORD.dep
#>  [1,] 0.0002675839 0.0002675839 0.0002675839 2.091542e-03
#>  [2,] 0.0024664457 0.0026615183 0.0011285264 1.002025e-02
#>  [3,] 0.0005732818 0.0005787961 0.0002823266 1.677763e-03
#>  [4,] 0.0004872805 0.0004929725 0.0002394680 6.262659e-04
#>  [5,] 0.0004059066 0.0004099744 0.0001998165 3.201787e-04
#>  [6,] 0.0003447286 0.0003475734 0.0001700069 1.933725e-04
#>  [7,] 0.0002986627 0.0003006772 0.0001475152 1.293954e-04
#>  [8,] 0.0029389397 0.0072216015 0.0014680343 1.548152e-04
#>  [9,] 0.0008168502 0.0015704700 0.0014680343 1.166482e-04
#> [10,] 0.0033835974 0.0091593329 0.0017451837 1.422614e-04
#> [11,] 0.0011873999 0.0019918653 0.0006438778 1.145119e-04
#> [12,] 0.0010225858 0.0016965126 0.0006438778 9.432408e-05
#> [13,] 0.0008785607 0.0014108836 0.0006438778 7.916885e-05
#> [14,] 0.0007679398 0.0011961369 0.0005497556 6.749313e-05
#> [15,] 0.0006820264 0.0010347488 0.0005497556 5.830055e-05
```

By default $w_{0} = \alpha/10$ and (for LORD 3 and LORD dep)
$b0 = alpha - w0$, but these parameters can optionally be specified by
the user subject to the requirements that $0 \leq w_{0} \leq \alpha$,
$b_{0} > 0$ and $w_{0} + b_{0} \leq \alpha$.

The value of `gammai` is also supplied by default, but can optionally be
specified by the user (as described above, see the formula for
$\gamma_{j}$[here](#LORDdep_xi) for version=‘dep’ and
[here](#LORD_gamma) for all other versions of LORD).

### SAFFRON

By default $w_{0} = \alpha/2$ and $\lambda = 0.5$, but these parameters
can optionally be specified by the user subject to the requirements that
$0 \leq w_{0} \leq \alpha$ and $0 < \lambda < 1$. The values of `gammai`
are also supplied by default, but can optionally be specified by the
user (as described above, see the formula for
$\gamma_{j}$[here](#SAFFRON_gamma)).

### ADDIS

By default $w_{0} = \alpha/2$, $\tau = 0.5$ and $\lambda = 0.25$, but
these parameters can optionally be specified by the user subject to the
requirements that $0 \leq w_{0} < \alpha$, $0 < \tau < 1$ and
$0 < \lambda < \tau$. The values of `gammai` are also supplied by
default, but can optionally be specified by the user.

### Alpha-spending and online fallback

The values of `gammai` are supplied by default, but can optionally be
specified by the user.

### ADDIS-spending

By default $\lambda = 0.25$ and $\tau = 0.5$, but these parameters can
optionally be specified by the user subject to the requirements that
$\lambda < \tau$, $0 < \lambda < 1$ and $0 < \tau < 1$. The values of
`gammai` are also supplied by default, but can optionally be specified
by the user.

### Asynchronous testing

Zrnic *et al.* (2021) proposed procedures to control the modified FDR
(mFDR) in the context of *asynchronous* testing, i.e. where each
hypothesis test can itself be a sequential process and the tests can
overlap in time. They presented asynchronous versions of the LOND, LORD
and SAFFRON procedures for a variety of trial settings, including the
following:

1: **Asynchronous online mFDR control**: This is for an asynchronous
testing process, consisting of tests that start and finish at
(potentially) random times. The discretised finish times of the test
correspond to the decision times.

2: **Online mFDR control under local dependence**: For any $t > 0$ we
allow the p-value $p_{t}$ to have arbitrary dependence on the previous
$L_{t}$ p-values. The fixed sequence $L_{t}$ is referred to as \`lags’.

3: **mFDR control in asynchronous mini-batch testing**: A mini-batch
represents a grouping of tests run asynchronously which result in
dependent p-values. Once a mini-batch of tests is fully completed, a new
one can start, testing hypotheses independent of the previous batch.
