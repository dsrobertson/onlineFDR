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
#' Version: \tab 1.1.2\cr
#' Date: \tab 2019-04-12\cr
#' License: \tab GPL-3 \cr
#' }
#'
#' Javanmard and Montanari (2015, 2018) proposed two methods for online FDR
#' control. The first is LORD, which stands for (significance)
#' Levels based On Recent Discovery and is implemented by the function
#' \code{\link{LORD}}. This function also includes the extension to the LORD 
#' procedure, called LORD++, proposed by Ramdas et al. (2017).
#' \code{\link{LORDdep}} provides a modified version of LORD that is valid for
#' dependent p-values.
#'
#' The second method is LOND, which stands for (significance) Levels based On
#' Number of Discoveries and is implemented by the function \code{\link{LOND}}.
#' By specifying \code{dep = TRUE}, thus function runs a modified version of
#' LOND which is valid for dependent p-values.
#'
#' As an alternative to these methods, a Bonferroni-like test is implemented
#' by the function \code{\link{bonfInfinite}}. This procedure is also valid
#' for dependent p-values.
#'
#' Further details on all these procedures can be found in Javanmard and
#' Montanari (2015, 2018) and Ramdas et al. (2017).
#'
#' @author David Robertson (\email{david.robertson@@mrc-bsu.cam.ac.uk}), Adel
#' Javanmard, Andrea Montanari and Natasha Karp.
#'
#' @references
#' Javanmard, A. and Montanari, A. (2015) On Online Control of False Discovery
#' Rate. \emph{arXiv preprint}, \url{https://arxiv.org/abs/1502.06197}
#'
#' Javanmard, A. and Montanari, A. (2018) Online Rules for Control of False
#' Discovery Rate and False Discovery Exceedance. \emph{Annals of Statistics},
#' 46(2):526-554.
#' 
#' Ramdas, A. et al. (2017). Online control of the false discovery rate with
#' decaying memory. \emph{Advances in Neural Information Processing Systems 30},
#' 5650-5659.
#'
#' @docType package
#' @name onlineFDR-package
NULL
