#' ADDIS: Adaptive discarding algorithm for online FDR control
#'
#' Implements the ADDIS algorithm for online FDR control, where ADDIS stands for
#' an ADaptive algorithm that DIScards conservative nulls, as presented by Tian
#' and Ramdas (2019). The algorithm compensates for the power loss of SAFFRON
#' with conservative nulls, by including both adaptivity in the fraction of
#' null hypotheses (like SAFFRON) and the conservativeness of nulls (unlike
#' SAFFRON).
#'
#' The function takes as its input either a vector of p-values, or a dataframe
#' with three columns. The dataframe requires an identifier (`id'), 
#' date (`date') and p-value (`pval'). If the asynchronous
#' version is specified (see below), then the column date should be replaced 
#' by the decision times. 
#' 
#' Given an overall significance level \eqn{\alpha}, ADDIS depends on constants
#' \eqn{w_0}, \eqn{\lambda} and \eqn{\tau}. Here \eqn{w_0} represents the
#' initial `wealth' of the procedure and satisfies \eqn{0 \le w_0 \le \alpha}.
#' \eqn{\tau \in (0,1)} represents the threshold for a
#' hypothesis to be selected for testing: p-values greater than \eqn{\tau} are
#' implicitly `discarded' by the procedure. Finally, \eqn{\lambda \in [0,\tau)}
#' sets the threshold for a p-value to be a candidate for rejection: ADDIS will
#' never reject a p-value larger than \eqn{\lambda}. The algorithm also
#' require a sequence of non-negative non-increasing numbers \eqn{\gamma_i} that
#' sum to 1.
#'
#' The ADDIS procedure provably controls the FDR for independent p-values. Tian
#' and Ramdas (2019) also presented a version for an asynchronous testing
#' process, consisting of tests that start and finish at (potentially) random
#' times. The discretised finish times of the test correspond to the decision
#' times. These decision times are given as the input \code{decision.times}.
#' Note that this asynchronous version controls a modified version of the FDR.
#' 
#' Further details of the ADDIS algorithms can be found in Tian and Ramdas
#' (2019).
#'
#' @param d Either a vector of p-values, or a dataframe with three columns: an
#'   identifier (`id'), p-value (`pval'), and decision times 
#'   (`decision.times').
#'
#' @param alpha Overall significance level of the procedure, the default is
#'   0.05.
#'
#' @param gammai Optional vector of \eqn{\gamma_i}. A default is provided with
#'   \eqn{\gamma_j} proportional to \eqn{1/j^(1.6)}.
#'
#' @param w0 Initial `wealth' of the procedure, defaults to \eqn{\alpha/2}.
#'   
#' @param tau Optional threshold for hypotheses to be selected for testing. Must
#'   be between 0 and 1, defaults to 0.5.
#'
#' @param lambda Optional parameter that sets the threshold for `candidate'
#'   hypotheses. Must be between 0 and tau, defaults to 0.25.
#'
#' @param async Logical. If \code{TRUE} runs the version for an asynchronous
#'   testing process. Defaults to FALSE.
#'   
#' @param random Logical. If \code{TRUE} (the default), then the order of the
#'   p-values in each batch (i.e. those that have exactly the same date) is
#'   randomised. Only needed if async=FALSE.
#'   
#' @param display_progress Logical. If \code{TRUE} prints out a progress bar for the algorithm runtime. 
#'   
#' @param date.format Optional string giving the format that is used for dates.
#'
#' @return \item{out}{A dataframe with the original p-values \code{pval}, the
#'   adjusted testing levels \eqn{\alpha_i} and the indicator function of
#'   discoveries \code{R}. Hypothesis \eqn{i} is rejected if the \eqn{i}-th
#'   p-value is less than or equal to \eqn{\alpha_i}, in which case \code{R[i] =
#'   1}  (otherwise \code{R[i] = 0}).}
#'
#'
#' @references Tian, J. and Ramdas, A. (2019). ADDIS: an adaptive discarding
#'   algorithm for online FDR control with conservative nulls.
#'   \emph{Advances in Neural Information Processing Systems}, 9388-9396.
#'
#'
#' @seealso
#'
#' ADDIS is identical to \code{\link{SAFFRON}} with option \code{discard=TRUE}.
#'
#' ADDIS with option \code{async=TRUE} is identical to \code{\link{SAFFRONstar}}
#' with option \code{discard=TRUE}.
#'
#'
#' @examples
#' sample.df1 <- data.frame(
#' id = c('A15432', 'B90969', 'C18705', 'B49731', 'E99902',
#'     'C38292', 'A30619', 'D46627', 'E29198', 'A41418',
#'     'D51456', 'C88669', 'E03673', 'A63155', 'B66033'),
#' date = as.Date(c(rep('2014-12-01',3),
#'                 rep('2015-09-21',5),
#'                 rep('2016-05-19',2),
#'                 '2016-11-12',
#'                 rep('2017-03-27',4))),
#' pval = c(2.90e-08, 0.06743, 0.01514, 0.08174, 0.00171,
#'         3.60e-05, 0.79149, 0.27201, 0.28295, 7.59e-08,
#'         0.69274, 0.30443, 0.00136, 0.72342, 0.54757))
#' 
#' ADDIS(sample.df1, random=FALSE)
#' 
#' 
#' sample.df2 <- data.frame(
#' id = c('A15432', 'B90969', 'C18705', 'B49731', 'E99902',
#'     'C38292', 'A30619', 'D46627', 'E29198', 'A41418',
#'     'D51456', 'C88669', 'E03673', 'A63155', 'B66033'),
#' pval = c(2.90e-08, 0.06743, 0.01514, 0.08174, 0.00171,
#'         3.60e-05, 0.79149, 0.27201, 0.28295, 7.59e-08,
#'         0.69274, 0.30443, 0.00136, 0.72342, 0.54757),
#' decision.times = seq_len(15) + 1)
#'
#' ADDIS(sample.df2, async = TRUE) # Asynchronous
#' 
#'
#' @export

ADDIS <- function(d, alpha = 0.05, gammai, w0, lambda = 0.25, tau = 0.5,
                  async = FALSE, random = TRUE, display_progress = FALSE, date.format = "%Y-%m-%d") {
    
    d <- checkPval(d)
    
    if (is.data.frame(d)) {
        pval <- d$pval
        if (async) {
            if (!("decision.times" %in% colnames(d))) {
                stop("d needs to have a column of decision.times")
            }
        } else {
            d <- checkdf(d, random, date.format)
            pval <- d$pval
        }
    } else if (is.vector(d)) {
        pval <- d
        if (async) {
            stop("d needs to be a dataframe with a column of decision.times")
        }
    } else {
        stop("d must either be a dataframe or a vector of p-values.")
    }
    
    if (alpha <= 0 || alpha > 1) {
        stop("alpha must be between 0 and 1.")
    }
    
    if (tau <= 0 || tau > 1) {
      stop("tau must be between 0 and 1.")
    }
    
    if (lambda <= 0 || lambda > tau) {
        stop("lambda must be between 0 and tau.")
    }
    
    N <- length(pval)
    
    if (missing(gammai)) {
        gammai <- 0.4374901658/(seq_len(N + 1)^(1.6))
    } else if (any(gammai < 0)) {
        stop("All elements of gammai must be non-negative.")
    } else if (sum(gammai) > 1) {
        stop("The sum of the elements of gammai must not be greater than 1.")
    }
    
    if (missing(w0)) {
        w0 = alpha/2
    } else if (w0 < 0) {
        stop("w0 must be non-negative.")
    } else if (w0 > alpha) {
        stop("w0 must be less than alpha.")
    }
    
    if (!(async)) {
        
            out <- addis_sync_faster(pval,
                                     gammai,
                                     lambda = lambda,
                                     alpha = alpha,
                                     tau = tau,
                                     w0 = w0,
                                     display_progress = display_progress)
            out$R <- as.numeric(out$R)
            if(is.data.frame(d) && !is.null(d$id)) {
                out$id <- d$id
            }
            out
        
    } else {
        
        if (any(is.na(d$decision.times))) {
            stop("Please provide a decision time for each p-value.")
        }
        
        E <- d$decision.times
        
        out <- addis_async_faster(pval, 
                                  E,
                                  gammai,
                                  lambda = lambda,
                                  alpha = alpha,
                                  tau = tau,
                                  w0 = w0,
                                  display_progress = display_progress)
        out$R <- as.numeric(out$R)
        if(is.data.frame(d) && !is.null(d$id)) {
            out$id <- d$id
        }
        out
    }
}
