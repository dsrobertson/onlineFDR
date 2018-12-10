## Introduction
This package allows users to control the false discovery rate for online
hypothesis testing, where hypotheses arrive sequentially in a stream, as
presented by Javanmard and Montanari (2015, 2018) and Ramdas et al. (2017). 
In this framework, a null hypothesis is rejected based only on the previous
decisions, as the future p-values and the number of hypotheses to be tested
are unknown.  

## Installation
To install the onlineFDR package in R, please run the following code:
```{r}
library(devtools)
install_github("dsrobertson/onlineFDR")
```

## References
Javanmard, A., and Montanari, A. (2015). On Online Control of False
Discovery Rate. *arXiv preprint*, https://arxiv.org/abs/1502.06197.

Javanmard, A., and Montanari, A. (2018). Online Rules for Control of False
Discovery Rate and False Discovery Exceedance. *Annals of Statistics*,
46(2):526-554.

Ramdas, A. et al. (2017). Online control of the false discovery rate with
decaying memory. *Advances in Neural Information Processing Systems 30*,
5650-5659.
