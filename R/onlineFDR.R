#' onlineFDR: A package for online FDR control
#'
#' The onlineFDR package provides methods to control the false discovery rate
#' (FDR) or familywise error rate (FWER) for online hypothesis testing, where
#' hypotheses arrive sequentially in a stream. A null hypothesis is rejected
#' based only on the previous decisions, as the future p-values and the number
#' of hypotheses to be tested are unknown.
#'
#' @section Details:
#'
#' \tabular{ll}{
#' Package: \tab onlineFDR \cr
#' Type: \tab Package\cr
#' Version: \tab 1.3.7\cr
#' Date: \tab 2019-10-14\cr
#' License: \tab GPL-3 \cr
#' }
#'
#' Javanmard and Montanari (2015, 2018) proposed two methods for online FDR
#' control. The first is LORD, which stands for (significance) Levels based On
#' Recent Discovery and is implemented by the function \code{\link{LORD}}. This
#' function also includes the extension to the LORD procedure, called LORD++
#' (\code{version='++'}), proposed by Ramdas et al. (2017). Setting
#' \code{version='discard'} implements a modified version of LORD that can
#' improve the power of the procedure in the presence of conservative nulls by
#' adaptively `discarding' these p-values, as proposed by Tian and Ramdas
#' (2019a). All these LORD procedures provably control the FDR under
#' independence of the p-values. However, setting \code{version='dep'} provides
#' a modified version of LORD that is valid for arbitrarily dependent p-values.
#' 
#' The second method is LOND, which stands for (significance) Levels based On
#' Number of Discoveries and is implemented by the function \code{\link{LOND}}.
#' This procedure controls the FDR under independence of the p-values, but the
#' slightly modified version of LOND proposed by Zrnic et al. (2018) also
#' provably controls the FDR under positive dependence (PRDS conditioN). In
#' addition, by specifying \code{dep = TRUE}, thus function runs a modified
#' version of LOND which is valid for arbitrarily dependent p-values.
#' 
#' Another method for online FDR control proposed by Ramdas et al. (2018) is the
#' \code{\link{SAFFRON}} procedure, which stands for Serial estimate of the
#' Alpha Fraction that is Futiley Rationed On true Null hypotheses. This
#' provides an adaptive algorithm for online FDR control. SAFFRON is related to 
#' the Alpha-Investing procedure of Foster and Stine (2008), a monotone version
#' of which is implemented by the function \code{\link{AlphaInvesting}}. Both 
#' these procedure provably control the FDR under independence of the p-values.
#' 
#' Zrnic et al. (2018) generalised these algorithms (i.e. LOND, LORD and
#' SAFFRON) for the context of asynchronous online testing, where each
#' hypothesis test can itself be a sequential process and the tests can overlap
#' in time. These algorithms are designed for the control of the modified FDR
#' (mFDR). They are implemented by the functions \code{\link{LONDstar}},
#' \code{\link{LORDstar}} and \code{\link{SAFFRONstar}}. Zrnic et al. (2018)
#' presented three explicit versions of these algorithms:
#' 
#' 1) \code{version='async'} is for an asynchoronous testing
#' process, consisting of tests that start and finish at (potentially) random 
#' times. The discretised finish times of the test correspond to the decision 
#' times.
#' 
#' 2) \code{version='dep'} is for online testing under local
#' dependence of the p-values. More precisely, for any \eqn{t>0} we allow the
#' p-value \eqn{p_t} to have arbitrary dependence on the previous \eqn{L_t}
#' p-values. The fixed sequence \eqn{L_t} is referred to as `lags'.
#' 
#' 3) \code{version='batch'} is for controlling the mFDR in
#' mini-batch testing, where a mini-batch represents a grouping of tests run
#' asynchronously which result in dependent p-values. Once a mini-batch of tests
#' is fully completed, a new one can start, testing hypotheses independent of
#' the previous batch.
#' 
#' Meanwhile, Tian and Ramdas (2019a) proposed the \code{\link{ADDIS}}
#' algorithm, which stands for an ADaptive algorithm that DIScards conservative
#' nulls. The algorithm compensates for the power loss of SAFFRON with
#' conservative nulls, by including both adapativity in the fraction of null
#' hypotheses (like SAFFRON) and the conservativeness of nulls (unlike SAFFRON).
#' The ADDIS procedure provably controls the FDR for independent p-values. Tian
#' and Ramdas (2019) also presented a version for an asynchoronous testing
#' process, consisting of tests that start and finish at (potentially) random
#' times.
#'
#' Finally, Tian and Ramdas (2019b) proposed a number of algorithms for online
#' FWER control. The only previously existing procedure for online FWER control
#' is Alpha-spending, which is an online analog of the Bonferroni procedure.
#' This is implemented by the function \code{\link{AlphaSpending}}, and provides
#' strong FWER control for arbitrarily dependent p-values. A uniformly more
#' powerful method is \code{\link{AlphaSaving}}, whcich again strongly controls
#' the FWER even under arbitrary dependence amongst the p-values. The
#' \code{\link{ADDIS.spending}} procedure compensates for the power loss of
#' Alpha-spending and Alpha-saving, by including both adapativity in the
#' fraction of null hypotheses and the conservativeness of nulls. This procedure
#' controls the FWER in the strong sense for independent p-values. Tian and
#' Ramdas (2019b) also presented a version for handling local dependence, which
#' can be specified by setting \code{dep=TRUE}.
#'
#' Further details on all these procedures can be found in Javanmard and
#' Montanari (2015, 2018), Ramdas et al. (2017, 2018), Robertson and Wason
#' (2018), Tian and Ramdas (2019a, 2019b) and Zrnic et al. (2018).
#'
#' @author David S. Robertson (\email{david.robertson@@mrc-bsu.cam.ac.uk}),
#' Adel Javanmard, Aaditya Ramdas, Jinjin Tian, Tijana Zrnic, Andrea Montanari
#' and Natasha A. Karp.
#'
#' @references
#' 
#' Aharoni, E. and Rosset, S. (2014). Generalized \eqn{\alpha}-investing:
#' definitions, optimality results and applications to publci databases.
#' \emph{Journal of the Royal Statistical Society (Series B)}, 76(4):771--794.
#' 
#' Foster, D. and Stine R. (2008). \eqn{\alpha}-investing: a procedure for 
#' sequential control of expected false discoveries. \emph{Journal of the Royal
#' Statistical Society (Series B)}, 29(4):429-444.
#' 
#' Javanmard, A. and Montanari, A. (2015) On Online Control of False Discovery
#' Rate. \emph{arXiv preprint}, \url{https://arxiv.org/abs/1502.06197}.
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
#' Tian, J. and Ramdas, A. (2019a). ADDIS: an adaptive discarding algorithm for
#' online FDR control with conservative nulls. \emph{arXiv preprint},
#' \url{https://arxiv.org/abs/1905.11465}.
#'
#' Tian, J. and Ramdas, A. (2019b). Online control of the familywise error rate.
#' \emph{arXiv preprint}, \url{https://arxiv.org/abs/1910.04900}.
#'
#' Zrnic, T. et al. (2018). Asynchronous Online Testing of Multiple Hypotheses.
#' \emph{arXiv preprint}, \url{https://arxiv.org/abs/1812.05068}.
#' 
#' @docType package
#' @name onlineFDR-package
NULL
