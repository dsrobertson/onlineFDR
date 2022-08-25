<!-- badges: start -->

[![R-CMD-check](https://github.com/dsrobertson/onlineFDR/workflows/R-CMD-check/badge.svg)](https://github.com/dsrobertson/onlineFDR/actions)
[![codecov](https://codecov.io/gh/dsrobertson/onlineFDR/branch/master/graph/badge.svg)](https://codecov.io/gh/dsrobertson/onlineFDR)
<!-- badges: end -->

# onlineFDR <img src="man/figures/logo.png" align="right" />

`onlineFDR` allows users to control the false discovery rate (FDR) or
familywise error rate (FWER) for online hypothesis testing, where
hypotheses arrive sequentially in a stream. In this framework, a null
hypothesis is rejected based on the evidence against it and the previous
rejection decisions.

## Installation

To install the latest (development) version of the onlineFDR package
from Bioconductor, please run the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# The following initializes usage of Bioc
BiocManager::install()

BiocManager::install("onlineFDR")
```

Alternatively, you can install the package directly from GitHub:

``` r
# install.packages("devtools") # If devtools not installed

devtools::install_github("dsrobertson/onlineFDR")
```

## Documentation

Documentation is hosted at <https://dsrobertson.github.io/onlineFDR/>

To view the vignette for the version of this package installed in your
system, start R and enter:

``` r
browseVignettes("onlineFDR")
```

## References

Aharoni, E. and Rosset, S. (2014). Generalized alpha-investing:
definitions, optimality results and applications to public databases.
*Journal of the Royal Statistical Society (Series B)*, 76(4):771–794.

Foster, D. and Stine R. (2008). alpha-investing: a procedure for
sequential control of expected false discoveries. *Journal of the Royal
Statistical Society (Series B)*, 29(4):429-444.

Javanmard, A., and Montanari, A. (2015). On Online Control of False
Discovery Rate. *arXiv preprint*, <https://arxiv.org/abs/1502.06197>.

Javanmard, A., and Montanari, A. (2018). Online Rules for Control of
False Discovery Rate and False Discovery Exceedance. *Annals of
Statistics*, 46(2):526-554.

Ramdas, A., Yang, F., Wainwright M.J. and Jordan, M.I. (2017). Online
control of the false discovery rate with decaying memory. *Advances in
Neural Information Processing Systems 30*, 5650-5659.

Ramdas, A., Zrnic, T., Wainwright M.J. and Jordan, M.I. (2018). SAFFRON:
an adaptive algorithm for online control of the false discovery rate.
*Proceedings of the 35th International Conference in Machine Learning*,
80:4286-4294.

Robertson, D.S. and Wason, J.M.S. (2018). Online control of the false
discovery rate in biomedical research. *arXiv preprint*,
<https://arxiv.org/abs/1809.07292>.

Robertson, D.S., Wason, J.M.S. and Ramdas, A. (2022). Online multiple
hypothesis testing for reproducible research.*arXiv preprint*,
<https://arxiv.org/abs/2208.11418>.

Robertson, D.S., Wildenhain, J., Javanmard, A. and Karp, N.A. (2019).
onlineFDR: an R package to control the false discovery rate for growing
data repositories. *Bioinformatics*, 35:4196-4199,
<https://doi.org/10.1093/bioinformatics/btz191>.

Tian, J. and Ramdas, A. (2019). ADDIS: an adaptive discarding algorithm
for online FDR control with conservative nulls. *Advances in Neural
Information Processing Systems*, 9388-9396.

Tian, J. and Ramdas, A. (2021). Online control of the familywise error
rate. *Statistical Methods for Medical Research*, 30(4):976–993.

Zrnic, T., Jiang D., Ramdas A. and Jordan M. (2020). The Power of
Batching in Multiple Hypothesis Testing. *International Conference on
Artificial Intelligence and Statistics*, PMLR, 108:3806-3815.

Zrnic, T., Ramdas, A. and Jordan, M.I. (2021). Asynchronous Online
Testing of Multiple Hypotheses. *Journal of Machine Learning Research*, 22:1-33.
