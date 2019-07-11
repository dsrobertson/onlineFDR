[![Build Status](https://travis-ci.org/dsrobertson/onlineFDR.svg?branch=master)](https://travis-ci.org/dsrobertson/onlineFDR)

# onlineFDR

## Introduction
This package allows users to control the false discovery rate for online
hypothesis testing, where hypotheses arrive sequentially in a stream, as
presented by Javanmard and Montanari (2015, 2018), Ramdas et al. (2017, 2018),
Tian and Ramdas (2019) and Zrnic et al. (2018). In this framework, a null
hypothesis is rejected based only on the previous decisions, as the future
p-values and the number of hypotheses to be tested are unknown.

## Installation
To install the latest (development) version of the onlineFDR package in R,
please run the following code:
```{r}
devtools::install_github("dsrobertson/onlineFDR")
```

## Documentation
The documentation is hosted at https://dsrobertson.github.io/onlineFDR/

## References
Javanmard, A., and Montanari, A. (2015). On Online Control of False
Discovery Rate. *arXiv preprint*, https://arxiv.org/abs/1502.06197.

Javanmard, A., and Montanari, A. (2018). Online Rules for Control of False
Discovery Rate and False Discovery Exceedance. *Annals of Statistics*,
46(2):526-554.

Ramdas, A., Yang, F., Wainwright M.J. and Jordan, M.I. (2017). Online control
of the false discovery rate with decaying memory. 
*Advances in Neural Information Processing Systems 30*, 5650-5659.

Ramdas, A., Zrnic, T., Wainwright M.J. and Jordan, M.I. (2018). SAFFRON: an
adaptive algorithm for online control of the false discovery rate. 
*Proceedings of the 35th International Conference in Machine Learning*,
80:4286-4294.

Robertson, D.S. and Wason, J.M.S. (2018). Online control of the false discovery
rate in biomedical research. *arXiv preprint*, https://arxiv.org/abs/1809.07292.

Robertson, D.S., Wildenhain, J., Javanmard, A. and Karp, N.A. (2019). onlineFDR:
an R package to control the false discovery rate for growing data repositories.
*Bioinformatics*, https://doi.org/10.1093/bioinformatics/btz191.

Tian, J. and Ramdas, A. (2019). ADDIS: an adaptive discarding algorithm for 
online FDR control with conservative nulls. *arXiv preprint*, 
https://arxiv.org/abs/1905.11465. 

Zrnic, T., Ramdas, A. and Jordan, M.I. (2018). Asynchronous Online Testing of
Multiple Hypotheses. *arXiv preprint*, https://arxiv.org/abs/1812.05068.
