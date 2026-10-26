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

Finally, Tian and Ramdas (2021) proposed procedures that provide online
control of the FWER. One procedure, online fallback, gives a uniform
improvement to the naive Alpha-spending procedure (see below). The
ADDIS-spending procedure compensates for the power loss of these
procedures by including both adaptivity in the fraction of null
hypotheses and the conservativeness of nulls.

## Variations to the default options

In the following section, we consider the arguments that a typical user
might consider amending for their analysis.

### Common arguments

As a default, the `alpha` argument is set to 0.05, where `alpha` sets
the overall significance level of the FDR or FWER controlling procedure.
By convention, the standard significance level utilised is 5%. However,
there are applications where an alternate threshold could be considered.
For example, a more stringent threshold might be appropriate when there
are limited resources to follow up significant findings. A less
stringent threshold might be appropriate when the downstream analysis is
a global analysis which can tolerate a higher proportion of false
positives.

To ensure correct interpretation of the dates provided there is a
date.format argument. As a default, the date format is set to receive
dates as year-month(01-12)-day(number). The following website provides
clear guidance on symbols used to interpret the date information:
<https://stat.ethz.ch/R-manual/R-devel/library/base/html/strptime.html>

As a default, the `random` argument is set to `TRUE`. In this situation,
the order of p-values in each batch (i.e. with the same date) is
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

set.seed(1); results.indep <- LOND(sample.df)    # independent or PRDS
set.seed(1); results.dep <- LOND(sample.df, dep=TRUE)   # arbitrary dependence

# compare adjusted significance thresholds
cbind(default = results.indep$alphai, dep = results.dep$alphai)
#>            default          dep
#>  [1,] 0.0027407722 0.0027407722
#>  [2,] 0.0011920626 0.0007947084
#>  [3,] 0.0010153042 0.0005538023
#>  [4,] 0.0008443651 0.0004052952
#>  [5,] 0.0007158466 0.0003135095
#>  [6,] 0.0006192614 0.0002527598
#>  [7,] 0.0005448529 0.0002101361
#>  [8,] 0.0007290564 0.0002682468
#>  [9,] 0.0006577249 0.0002324964
#> [10,] 0.0007985699 0.0002726455
#> [11,] 0.0007328818 0.0002426860
#> [12,] 0.0006770682 0.0002181831
#> [13,] 0.0006290720 0.0001978131
#> [14,] 0.0005873666 0.0001806413
#> [15,] 0.0005507966 0.0001659911
```

The vector `betai` is supplied by default, but can optionally be
specified by the user (see the formula for $`\beta_j`$ in the [theory
vignette](https://dsrobertson.github.io/onlineFDR/articles/theory.html#LOND_beta)).

### LORD

The default version of LORD used is version ‘++’, but the user can
optionally specify versions 3, ‘discard’ and ‘dep’ using the `version`
argument (see
[here](https://dsrobertson.github.io/onlineFDR/articles/theory.html#LORD)
for further details about the different versions).

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
#>  [1,] 0.0002740772 0.0002740772 0.0002740772 2.091536e-03
#>  [2,] 0.0025262981 0.0027257486 0.0011559120 1.002023e-02
#>  [3,] 0.0005871934 0.0005927641 0.0002891777 1.677758e-03
#>  [4,] 0.0004991051 0.0005048694 0.0002452791 6.262642e-04
#>  [5,] 0.0004157566 0.0004198683 0.0002046653 3.201778e-04
#>  [6,] 0.0003530941 0.0003559613 0.0001741324 1.933720e-04
#>  [7,] 0.0003059103 0.0003079335 0.0001510949 1.293950e-04
#>  [8,] 0.0030102579 0.0049234562 0.0015036586 1.548148e-04
#>  [9,] 0.0008366724 0.0010706961 0.0015036586 1.166479e-04
#> [10,] 0.0034657061 0.0070615791 0.0017875335 1.422610e-04
#> [11,] 0.0012162141 0.0015356702 0.0006595026 1.145116e-04
#> [12,] 0.0010474005 0.0013079618 0.0006595026 9.432385e-05
#> [13,] 0.0008998804 0.0010877502 0.0006595026 7.916865e-05
#> [14,] 0.0007865752 0.0009221867 0.0005630964 6.749296e-05
#> [15,] 0.0006985769 0.0007977612 0.0005630964 5.830041e-05
```

By default $`w_0 = \alpha/10`$ and (for LORD 3 and LORD dep)
$`b_0 = \alpha - w_0`$, but these parameters can optionally be specified
by the user subject to the requirements that $`0 \leq w_0 \leq \alpha`$,
$`b_0 > 0`$ and $`w_0+b_0 \leq \alpha`$. For LORD with discarding
(`version = 'discard'`), the default is
$`w_0 = \min(\alpha/10, \tau\alpha/2)`$ and the requirement is
$`0 \leq w_0 \leq \tau\alpha`$, where $`\tau`$ is the discarding
threshold `tau.discard`.

The value of `gammai` is also supplied by default, but can optionally be
specified by the user (see the formula for
$`\xi_j`$[here](https://dsrobertson.github.io/onlineFDR/articles/theory.html#LORDdep_xi)
for version=‘dep’ and for
$`\gamma_j`$[here](https://dsrobertson.github.io/onlineFDR/articles/theory.html#LORD_gamma)
for all other versions of LORD).

### SAFFRON

By default $`w_0 = \alpha/2`$ and $`\lambda = 0.5`$, but these
parameters can optionally be specified by the user subject to the
requirements that $`0 \leq w_0 \leq \alpha`$ and $`0 < \lambda < 1`$.
The values of `gammai` are also supplied by default, but can optionally
be specified by the user (see the formula for
$`\gamma_j`$[here](https://dsrobertson.github.io/onlineFDR/articles/theory.html#SAFFRON_gamma)).

### ADDIS

By default $`w_0 = \alpha/2`$, $`\tau = 0.5`$ and $`\lambda = 0.25`$,
but these parameters can optionally be specified by the user subject to
the requirements that $`0 \leq w_0 \leq \alpha`$, $`0 < \tau \leq 1`$
and $`0 < \lambda < \tau`$. The values of `gammai` are also supplied by
default, but can optionally be specified by the user.

### Alpha-spending and online fallback

The values of `gammai` are supplied by default, but can optionally be
specified by the user.

### ADDIS-spending

By default $`\lambda = 0.25`$ and $`\tau = 0.5`$, but these parameters
can optionally be specified by the user subject to the requirements that
$`0 < \lambda < \tau \leq 1`$. The values of `gammai` are also supplied
by default, but can optionally be specified by the user.

### Asynchronous testing

Zrnic *et al.* (2021) proposed procedures to control the modified FDR
(mFDR) in the context of *asynchronous* testing, i.e. where each
hypothesis test can itself be a sequential process and the tests can
overlap in time. They presented asynchronous versions of the LOND, LORD
and SAFFRON procedures for a variety of trial settings, including the
following:

1: **Asynchronous online mFDR control**: This is for an asynchronous
testing process, consisting of tests that start and finish at
(potentially) random times. The discretised finish times of the tests
correspond to the decision times.

2: **Online mFDR control under local dependence**: For any $`t>0`$ we
allow the p-value $`p_t`$ to have arbitrary dependence on the previous
$`L_t`$ p-values. The fixed sequence $`L_t`$ is referred to as ‘lags’.

3: **mFDR control in asynchronous mini-batch testing**: A mini-batch
represents a grouping of tests run asynchronously which result in
dependent p-values. Once a mini-batch of tests is fully completed, a new
one can start, testing hypotheses independent of the previous batch.

## References

Foster, D. and Stine, R. (2008). $`\alpha`$-investing: a procedure for
sequential control of expected false discoveries. *Journal of the Royal
Statistical Society (Series B)*, 70(2):429-444.

Javanmard, A., and Montanari, A. (2015). On Online Control of False
Discovery Rate. *arXiv preprint*, <https://arxiv.org/abs/1502.06197>.

Javanmard, A., and Montanari, A. (2018). Online Rules for Control of
False Discovery Rate and False Discovery Exceedance. *Annals of
Statistics*, 46(2):526-554.

Li, A., and Barber, R.F. (2017). Accumulation Tests for FDR Control in
Ordered Hypothesis Testing. *Journal of the American Statistical
Association*, 112(518):837-849.

Ramdas, A., Yang, F., Wainwright, M.J. and Jordan, M.I. (2017). Online
control of the false discovery rate with decaying memory. *Advances in
Neural Information Processing Systems 30*, 5650-5659.

Ramdas, A., Zrnic, T., Wainwright, M.J. and Jordan, M.I. (2018).
SAFFRON: an adaptive algorithm for online control of the false discovery
rate. *Proceedings of the 35th International Conference on Machine
Learning*, 80:4286-4294.

Tian, J. and Ramdas, A. (2019). ADDIS: an adaptive discarding algorithm
for online FDR control with conservative nulls. *Advances in Neural
Information Processing Systems 32*, 9388-9396.

Tian, J. and Ramdas, A. (2021). Online control of the familywise error
rate. *Statistical Methods in Medical Research*, 30(4):976-993.

Zrnic, T., Ramdas, A. and Jordan, M.I. (2021). Asynchronous Online
Testing of Multiple Hypotheses. *Journal of Machine Learning Research*,
22(33):1-39.
