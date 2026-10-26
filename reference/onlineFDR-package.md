# onlineFDR: A package for online error rate control

The onlineFDR package provides methods to control the false discovery
rate (FDR) or familywise error rate (FWER) for online hypothesis
testing, where hypotheses arrive in a stream. A null hypothesis is
rejected based on the evidence against it and on the previous rejection
decisions.

## Details

|          |           |
|----------|-----------|
| Package: | onlineFDR |
| Type:    | Package   |
| License: | GPL-3     |

Javanmard and Montanari (2015, 2018) proposed two methods for online FDR
control. The first is LORD, which stands for (significance) Levels based
On Recent Discovery and is implemented by the function
[`LORD`](https://dsrobertson.github.io/onlineFDR/reference/LORD.md).
This function also includes the extension to the LORD procedure, called
LORD++ (`version='++'`), proposed by Ramdas et al. (2017). Setting
`version='discard'` implements a modified version of LORD that can
improve the power of the procedure in the presence of conservative nulls
by adaptively \`discarding' these p-values, as proposed by Tian and
Ramdas (2019). LORD++ and the discarding version provably control the
FDR under independence of the p-values; for LORD 3 this has only been
shown empirically. FDR control is not guaranteed for positively
dependent p-values in general (see
[`vignette("theory")`](https://dsrobertson.github.io/onlineFDR/articles/theory.md)).
Setting `version='dep'` provides a modified version of LORD that is
valid for arbitrarily dependent p-values.

The second method is LOND, which stands for (significance) Levels based
On Number of Discoveries and is implemented by the function
[`LOND`](https://dsrobertson.github.io/onlineFDR/reference/LOND.md).
This procedure controls the FDR under independence of the p-values. The
slightly modified version of LOND of Javanmard and Montanari (2015,
Remark 2.3) also provably controls the FDR under positive dependence
(PRDS condition), as shown by Zrnic et al. (2021), and their argument
shows that the original version does too (see
[`LOND`](https://dsrobertson.github.io/onlineFDR/reference/LOND.md)). In
addition, by specifying `dep = TRUE`, this function runs a modified
version of LOND which is valid for arbitrarily dependent p-values.

Another method for online FDR control proposed by Ramdas et al. (2018)
is the
[`SAFFRON`](https://dsrobertson.github.io/onlineFDR/reference/SAFFRON.md)
procedure, which stands for Serial estimate of the Alpha Fraction that
is Futilely Rationed On true Null hypotheses. This provides an adaptive
algorithm for online FDR control. SAFFRON is related to the
Alpha-investing procedure of Foster and Stine (2008), a monotone version
of which is implemented by the function
[`Alpha_investing`](https://dsrobertson.github.io/onlineFDR/reference/Alpha_investing.md).
Both these procedures provably control the FDR under independence of the
p-values, but not for positively dependent p-values in general, even
with the default settings (see
[`vignette("theory")`](https://dsrobertson.github.io/onlineFDR/articles/theory.md)).

Tian and Ramdas (2019) proposed the
[`ADDIS`](https://dsrobertson.github.io/onlineFDR/reference/ADDIS.md)
algorithm, which stands for an ADaptive algorithm that DIScards
conservative nulls. The algorithm compensates for the power loss of
SAFFRON with conservative nulls, by including both adaptivity in the
fraction of null hypotheses (like SAFFRON) and the conservativeness of
nulls (unlike SAFFRON). The ADDIS procedure provably controls the FDR
for independent p-values, but not for positively dependent p-values in
general. Tian and Ramdas (2019) also presented a version for an
asynchronous testing process, consisting of tests that start and finish
at (potentially) random times.

For testing batches of hypotheses, Zrnic et al. (2020) proposed batched
online testing algorithms that control the FDR, where the p-values
across different batches are independent, and within a batch the
p-values are either positively dependent or independent.

Zrnic et al. (2021) generalised LOND, LORD and SAFFRON for asynchronous
online testing, where each hypothesis test can itself be a sequential
process and the tests can overlap in time. Note though that these
algorithms are designed for the control of a modified FDR (mFDR),
although LONDstar, LORDstar and SAFFRONstar with `version='async'` also
control the FDR when the null p-values are independent of each other and
of the non-nulls, and the decision times do not depend on the p-values.
The `'dep'` and `'batch'` versions are only guaranteed to control the
mFDR. They are implemented by the functions
[`LONDstar`](https://dsrobertson.github.io/onlineFDR/reference/LONDstar.md),
[`LORDstar`](https://dsrobertson.github.io/onlineFDR/reference/LORDstar.md)
and
[`SAFFRONstar`](https://dsrobertson.github.io/onlineFDR/reference/SAFFRONstar.md).
Zrnic et al. (2021) presented three explicit versions of these
algorithms:

1\) `version='async'` is for an asynchronous testing process, consisting
of tests that start and finish at (potentially) random times. The
discretised finish times of the tests correspond to the decision times.

2\) `version='dep'` is for online testing under local dependence of the
p-values. More precisely, for any \\t\>0\\ we allow the p-value \\p_t\\
to have arbitrary dependence on the previous \\L_t\\ p-values. The fixed
sequence \\L_t\\ is referred to as \`lags'.

3\) `version='batch'` is for controlling the mFDR in mini-batch testing,
where a mini-batch represents a grouping of tests run asynchronously
which result in dependent p-values. Once a mini-batch of tests is fully
completed, a new one can start, testing hypotheses independent of the
previous batch.

Recently, Xu and Ramdas (2021) proposed the
[`supLORD`](https://dsrobertson.github.io/onlineFDR/reference/supLORD.md)
algorithm, which provably controls the false discovery exceedance (FDX)
for p-values that are conditionally superuniform under the null. supLORD
also controls the supFDR and hence the FDR (even at stopping times).

Finally, Tian and Ramdas (2021) proposed a number of algorithms for
online FWER control. The only previously existing procedure for online
FWER control was Alpha-spending, which is an online analogue of the
Bonferroni procedure. This is implemented by the function
[`Alpha_spending`](https://dsrobertson.github.io/onlineFDR/reference/Alpha_spending.md),
and provides strong FWER control for arbitrarily dependent p-values. A
uniformly more powerful method is
[`online_fallback`](https://dsrobertson.github.io/onlineFDR/reference/online_fallback.md),
which again strongly controls the FWER even under arbitrary dependence
amongst the p-values. The
[`ADDIS_spending`](https://dsrobertson.github.io/onlineFDR/reference/ADDIS_spending.md)
procedure compensates for the power loss of Alpha-spending and online
fallback, by including both adaptivity in the fraction of null
hypotheses and the conservativeness of nulls. This procedure controls
the FWER in the strong sense for independent p-values; this is not
guaranteed for positively dependent p-values in general. Tian and Ramdas
(2021) also presented a version for handling local dependence, which can
be specified by setting `dep=TRUE`. Fischer et al. (2024) proposed an
exhaustive version of ADDIS-spending, implemented by the function
[`ADDIS_exhaustive`](https://dsrobertson.github.io/onlineFDR/reference/ADDIS_exhaustive.md),
which is a uniform improvement of ADDIS-spending and also strongly
controls the FWER for independent p-values.

Further details on all these procedures can be found in Fischer et al.
(2024), Javanmard and Montanari (2015, 2018), Ramdas et al. (2017,
2018), Robertson and Wason (2018), Tian and Ramdas (2019, 2021), Xu and
Ramdas (2021), and Zrnic et al. (2020, 2021).

## References

Aharoni, E. and Rosset, S. (2014). Generalized \\\alpha\\-investing:
definitions, optimality results and application to public databases.
*Journal of the Royal Statistical Society (Series B)*, 76(4):771-794.

Fischer, L., Bofill Roig, M. and Brannath, W. (2024). An exhaustive
ADDIS principle for online FWER control. *Biometrical Journal*,
66(3):2300237.

Foster, D. and Stine, R. (2008). \\\alpha\\-investing: a procedure for
sequential control of expected false discoveries. *Journal of the Royal
Statistical Society (Series B)*, 70(2):429-444.

Javanmard, A. and Montanari, A. (2015). On Online Control of False
Discovery Rate. *arXiv preprint*, <https://arxiv.org/abs/1502.06197>.

Javanmard, A. and Montanari, A. (2018). Online Rules for Control of
False Discovery Rate and False Discovery Exceedance. *Annals of
Statistics*, 46(2):526-554.

Ramdas, A., Yang, F., Wainwright, M.J. and Jordan, M.I. (2017). Online
control of the false discovery rate with decaying memory. *Advances in
Neural Information Processing Systems 30*, 5650-5659.

Ramdas, A., Zrnic, T., Wainwright, M.J. and Jordan, M.I. (2018).
SAFFRON: an adaptive algorithm for online control of the false discovery
rate. *Proceedings of the 35th International Conference on Machine
Learning*, 80:4286-4294.

Robertson, D.S. and Wason, J.M.S. (2018). Online control of the false
discovery rate in biomedical research. *arXiv preprint*,
<https://arxiv.org/abs/1809.07292>.

Robertson, D.S., Wason, J.M.S. and Ramdas, A. (2023). Online multiple
hypothesis testing. *Statistical Science*, 38(4):557-575,
[doi:10.1214/23-STS901](https://doi.org/10.1214/23-STS901) .

Robertson, D.S., Wildenhain, J., Javanmard, A. and Karp, N.A. (2019).
onlineFDR: an R package to control the false discovery rate for growing
data repositories. *Bioinformatics*, 35(20):4196-4199,
[doi:10.1093/bioinformatics/btz191](https://doi.org/10.1093/bioinformatics/btz191)
.

Tian, J. and Ramdas, A. (2019). ADDIS: an adaptive discarding algorithm
for online FDR control with conservative nulls. *Advances in Neural
Information Processing Systems 32*, 9388-9396.

Tian, J. and Ramdas, A. (2021). Online control of the familywise error
rate. *Statistical Methods in Medical Research*, 30(4):976-993.

Xu, Z. and Ramdas, A. (2021). Dynamic Algorithms for Online Multiple
Testing. *Annual Conference on Mathematical and Scientific Machine
Learning*, PMLR, 145:955-986.

Zrnic, T., Jiang, D., Ramdas, A. and Jordan, M.I. (2020). The Power of
Batching in Multiple Hypothesis Testing. *International Conference on
Artificial Intelligence and Statistics*, PMLR, 108:3806-3815.

Zrnic, T., Ramdas, A. and Jordan, M.I. (2021). Asynchronous Online
Testing of Multiple Hypotheses. *Journal of Machine Learning Research*,
22(33):1-39.

## See also

Useful links:

- <https://dsrobertson.github.io/onlineFDR/>

- <https://dsrobertson.github.io/onlineFDR/index.html>

## Author

David S. Robertson (<david.robertson@mrc-bsu.cam.ac.uk>), Lathan Liou,
Adel Javanmard, Aaditya Ramdas, Jinjin Tian, Tijana Zrnic, Andrea
Montanari, Natasha A. Karp and Nikita Mozgunov.
