#' onlineFDR: A package for online FDR control
#'
#' The onlineFDR package provides methods to control the false discovery rate
#' (FDR) for online hypothesis testing, where hypotheses arrive sequentially
#' in a stream. A null hypothesis is rejected based only on the previous
#' decisions, as the future p-values and the number of hypotheses to be tested
#' are unknown.
#'
#' @section Details:
#'
#' \tabular{ll}{
#' Package: \tab onlineFDR \cr
#' Type: \tab Package\cr
#' Version: \tab 1.3.6\cr
#' Date: \tab 2019-09-18\cr
#' License: \tab GPL-3 \cr
#' }
#'
#' Javanmard and Montanari (2015, 2018) proposed two methods for online FDR
#' control. The first is LORD, which stands for (significance)
#' Levels based On Recent Discovery and is implemented by the function
#' \code{\link{LORD}}. This function also includes the extension to the LORD 
#' procedure, called LORD++, proposed by Ramdas et al. (2017).
#' Setting \code{version='dep'} provides a modified version of LORD that is
#' valid for dependent p-values. 
#'
#' The second method is LOND, which stands for (significance) Levels based On
#' Number of Discoveries and is implemented by the function \code{\link{LOND}}.
#' By specifying \code{dep = TRUE}, thus function runs a modified version of
#' LOND which is valid for arbitrarily dependent p-values.
#' 
#' A related method proposed by Ramdas et al. (2018) is the
#' \code{\link{SAFFRON}} procedure, which stands for Serial estimate of the
#' Alpha Fraction that is Futiley Rationed On true Null hypotheses. This
#' provides an adaptive algorithm for online FDR control. SAFFRON is related to 
#' the Alpha-Investing procedure of Foster and Stine (2008), a monotone version
#' of which is implemented by the function \code{\link{AlphaInvesting}}.
#' 
#' Zrnic et al. (2018) generalised these algorithms for the context of
#' asynchronous online testing, where each hypothesis test can itself be a
#' sequential process and the tests can overlap in time. These algorithms
#' are designed for the control of the modified FDR (mFDR).
#'
#' As an alternative to all these methods, a Bonferroni-like test is implemented
#' by the function \code{\link{bonfInfinite}}. This procedure is also valid
#' for dependent p-values.
#'
#' Further details on all these procedures can be found in Javanmard and
#' Montanari (2015, 2018), Ramdas et al. (2017, 2018) and Zrnic et al. (2018).
#'
#' @author David S. Robertson (\email{david.robertson@@mrc-bsu.cam.ac.uk}),
#' Adel Javanmard, Aaditya Ramdas, Jinjin Tian, Tijana Zrnic, Andrea Montanari
#' and Natasha A. Karp.
#'
#' @references
#' 
#' Foster, D. and Stine R. (2008). \eqn{\alpha}-investing: a procedure for 
#' sequential control of expected false discoveries. \emph{Journal of the Royal
#' Statistical Society (Series B)}, 29(4):429-444.
#' 
#' Javanmard, A. and Montanari, A. (2015) On Online Control of False Discovery
#' Rate. \emph{arXiv preprint}, \url{https://arxiv.org/abs/1502.06197}
#'
#' Javanmard, A. and Montanari, A. (2018) Online Rules for Control of False
#' Discovery Rate and False Discovery Exceedance. \emph{Annals of Statistics},
#' 46(2):526-554.
#'
#' Ramdas, A., Yang, F., Wainwright M.J. and Jordan, M.I. (2017). Online control
#' of the false discovery rate with decaying memory. \emph{Advances in Neural
#' Information Processing Systems 30}, 5650-5659.
#'
#' Ramdas, A., Zrnic, T., Wainwright M.J. and Jordan, M.I. (2018). SAFFRON: an
#' adaptive algorithm for online control of the false discovery rate.
#' \emph{Proceedings of the 35th International Conference in Machine Learning},
#' 80:4286-4294.
#'
#' Robertson, D.S. and Wason, J.M.S. (2018). Online control of the false
#' discovery rate in biomedical research. \emph{arXiv preprint},
#' \url{https://arxiv.org/abs/1809.07292}.
#'
#' Robertson, D.S., Wildenhain, J., Javanmard, A. and Karp, N.A. (2019).
#' onlineFDR: an R package to control the false discovery rate for growing data
#' repositories. \emph{Bioinformatics},
#' \url{https://doi.org/10.1093/bioinformatics/btz191}.
#'
#' Tian, J. and Ramdas, A. (2019). ADDIS: an adaptive discarding algorithm for
#' online FDR control with conservative nulls. \emph{arXiv preprint},
#' \url{https://arxiv.org/abs/1905.11465}.
#'
#' Zrnic, T. et al. (2018). Asynchronous Online Testing of Multiple Hypotheses.
#' \emph{arXiv preprint}, \url{https://arxiv.org/abs/1812.05068}
#'
#' @docType package
#' @name onlineFDR-package
NULL
