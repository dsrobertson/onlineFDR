#' SAFFRONstar: Adaptive online mFDR control for asynchronous testing
#'
#' Implements the SAFFRON algorithm for asynchronous online testing, as
#' presented by Zrnic et al. (2018).
#'
#' The function takes as its input either a vector of p-values, or a dataframe
#' with three columns: an identifier (`id'),
#' p-value (`pval'), or a column describing the conflict sets for the hypotheses. 
#' This takes the form of a vector of decision times or lags. Batch sizes can be 
#' specified as a separate argument (see below).
#'
#' Zrnic et al. (2018) present explicit three versions of SAFFRONstar:
#'
#' 1) \code{version='async'} is for an asynchronous testing process, consisting
#' of tests that start and finish at (potentially) random times. The discretised
#' finish times of the test correspond to the decision times. These decision
#' times are given as the input \code{decision.times} for this version of the
#' SAFFRONstar algorithm. For this version of SAFFRONstar, Tian and Ramdas
#' (2019) presented an algorithm that can improve the power of the procedure in
#' the presence of conservative nulls by adaptively `discarding' these p-values.
#' This can be called by setting the option \code{discard=TRUE}.
#'
#' 2) \code{version='dep'} is for online testing under local dependence of the
#' p-values. More precisely, for any \eqn{t>0} we allow the p-value \eqn{p_t} to
#' have arbitrary dependence on the previous \eqn{L_t} p-values. The fixed
#' sequence \eqn{L_t} is referred to as `lags', and is given as the input
#' \code{lags} for this version of the SAFFRONstar algorithm.
#'
#' 3) \code{version='batch'} is for controlling the mFDR in mini-batch testing,
#' where a mini-batch represents a grouping of tests run asynchronously which
#' result in dependent p-values. Once a mini-batch of tests is fully completed,
#' a new one can start, testing hypotheses independent of the previous batch.
#' The batch sizes are given as the input \code{batch.sizes} for this version of
#' the SAFFRONstar algorithm.
#'
#' Given an overall significance level \eqn{\alpha}, SAFFRONstar depends on
#' constants \eqn{w_0} and \eqn{\lambda}, where \eqn{w_0} satisfies \eqn{0 \le
#' w_0 \le (1 - \lambda)\alpha} and represents the intial `wealth' of the
#' procedure, and \eqn{0 < \lambda < 1} represents the threshold for a
#' `candidate' hypothesis. A `candidate' refers to p-values smaller than
#' \eqn{\lambda}, since SAFFRONstar will never reject a p-value larger than
#' \eqn{\lambda}. The algorithms also require a sequence of non-negative
#' non-increasing numbers \eqn{\gamma_i} that sum to 1.
#'
#' Note that these SAFFRONstar algorithms control the \emph{modified} FDR
#' (mFDR). The `async' version also controls the usual FDR if the p-values are
#' assumed to be independent.
#'
#' Further details of the SAFFRONstar algorithms can be found in Zrnic et al.
#' (2018).
#'
#'
#' @param d Either a vector of p-values, or a dataframe with three columns: an
#'   identifier (`id'), 
#'   p-value (`pval'), and either 
#'   `decision.times', or
#'   `lags', depending on which version you're using. See version for more details.
#'
#' @param alpha Overall significance level of the procedure, the default is
#'   0.05.
#'
#' @param gammai Optional vector of \eqn{\gamma_i}. A default is provided with
#'   \eqn{\gamma_j} proportional to \eqn{1/j^(1.6)}.
#'
#' @param version Takes values 'async', 'dep' or 'batch'. This specifies the
#'   version of SAFFRONstar to use. \code{version='async'} requires a 
#' column of decision times (`decision.times'). \code{version='dep'} requires a
#' column of lags (`lags').
#' \code{version='batch'} requires a vector of batch sizes (`batch.sizes').
#'
#' @param w0 Initial `wealth' of the procedure, defaults to \eqn{\alpha/10}.
#'
#' @param lambda Optional threshold for a `candidate' hypothesis, must be
#'   between 0 and 1. Defaults to 0.5.
#'
#' @param batch.sizes A vector of batch sizes, this is required for
#'   \code{version='batch'}.
#'
#' @param discard Logical. If \code{TRUE} then runs the ADDIS algorithm with
#'   adaptive discarding of conservative nulls. The default is \code{FALSE}.
#'
#' @param tau.discard Optional threshold for hypotheses to be selected for
#'   testing. Must be between 0 and 1, defaults to 0.5. This is required if
#'   \code{discard=TRUE}.
#'
#' @return \item{d.out}{A dataframe with the original p-values \code{pval}, the
#'   adjusted testing levels \eqn{\alpha_i} and the indicator function of
#'   discoveries \code{R}. Hypothesis \eqn{i} is rejected if the \eqn{i}-th
#'   p-value is less than or equal to \eqn{\alpha_i}, in which case \code{R[i] =
#'   1}  (otherwise \code{R[i] = 0}).}
#'
#'
#' @references Zrnic, T., Ramdas, A. and Jordan, M.I. (2018). Asynchronous
#'   Online Testing of Multiple Hypotheses. \emph{arXiv preprint},
#'   \url{https://arxiv.org/abs/1812.05068}.
#'
#'   Tian, J. and Ramdas, A. (2019). ADDIS: an adaptive discarding algorithm for
#'   online FDR control with conservative nulls. \emph{arXiv preprint},
#'   \url{https://arxiv.org/abs/1905.11465}.
#'
#'
#' @seealso
#'
#' \code{\link{SAFFRON}} presents versions of SAFFRON for \emph{synchronous}
#' p-values, i.e. where each test can only start when the previous test has
#' finished.
#'
#' If \code{version='async'} and \code{discard=TRUE}, then SAFFRONstar is
#' identical to \code{\link{ADDIS}} with option \code{async=TRUE}.
#'
#'
#' @examples
#' sample.df <- data.frame(
#' id = c('A15432', 'B90969', 'C18705', 'B49731', 'E99902',
#'     'C38292', 'A30619', 'D46627', 'E29198', 'A41418',
#'     'D51456', 'C88669', 'E03673', 'A63155', 'B66033'),
#' pval = c(2.90e-08, 0.06743, 0.01514, 0.08174, 0.00171,
#'         3.60e-05, 0.79149, 0.27201, 0.28295, 7.59e-08,
#'         0.69274, 0.30443, 0.00136, 0.72342, 0.54757),
#' decision.times = seq_len(15) + 1)
#'
#' SAFFRONstar(sample.df, version='async')
#' 
#' sample.df2 <- data.frame(
#' id = c('A15432', 'B90969', 'C18705', 'B49731', 'E99902',
#'     'C38292', 'A30619', 'D46627', 'E29198', 'A41418',
#'     'D51456', 'C88669', 'E03673', 'A63155', 'B66033'),
#' pval = c(2.90e-08, 0.06743, 0.01514, 0.08174, 0.00171,
#'         3.60e-05, 0.79149, 0.27201, 0.28295, 7.59e-08,
#'         0.69274, 0.30443, 0.00136, 0.72342, 0.54757),
#' lags = rep(1,15))
#' 
#' SAFFRONstar(sample.df2, version='dep')
#' 
#' sample.df3 <- data.frame(
#' id = c('A15432', 'B90969', 'C18705', 'B49731', 'E99902',
#'     'C38292', 'A30619', 'D46627', 'E29198', 'A41418',
#'     'D51456', 'C88669', 'E03673', 'A63155', 'B66033'),
#' pval = c(2.90e-08, 0.06743, 0.01514, 0.08174, 0.00171,
#'         3.60e-05, 0.79149, 0.27201, 0.28295, 7.59e-08,
#'         0.69274, 0.30443, 0.00136, 0.72342, 0.54757))
#' 
#' SAFFRONstar(sample.df3, version='batch', batch.sizes = c(4,6,5))
#'
#' @export

SAFFRONstar <- function(d, alpha = 0.05, version, gammai, w0, lambda = 0.5, batch.sizes, 
    discard = FALSE, tau.discard = 0.5) {
    
    d <- checkPval(d)
    
    if (is.data.frame(d)) {
        checkSTARdf(d, version)
        pval <- d$pval
    } else if (is.vector(d)) {
        pval <- d
        if(version == "async") {
            stop("d needs to be a dataframe with a column of decision.times")
        }
        else if(version == "dep") {
            stop("d needs to be a dataframe with a column of lags")
        }
    } else {
        stop("d must either be a dataframe or a vector of p-values.")
    }
    
    if (alpha <= 0 || alpha > 1) {
        stop("alpha must be between 0 and 1.")
    }
    
    if (lambda <= 0 || lambda > 1) {
        stop("lambda must be between 0 and 1.")
    }
    
    if (version == "async" && discard) {
        return(ADDIS(d, alpha, gammai, w0, lambda, tau.discard, async = TRUE))
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
        w0 = (1 - lambda) * alpha/2
    } else if (w0 < 0) {
        stop("w0 must be non-negative.")
    } else if (w0 >= (1 - lambda) * alpha) {
        stop("w0 must be less than (1-lambda)*alpha")
    }
    
    version <- checkStarVersion(d, N, version, batch.sizes)
    
    switch(version, {
        
        ## async = 1
        E <- d$decision.times
        out <- saffronstar_async_faster(pval, E, gammai)
        out$R <- as.numeric(out$R)
        out
        
    }, {
        ## dep = 2
        L <- d$lags
        out <- saffronstar_dep_faster(pval, L, gammai)
        out$R <- as.numeric(out$R)
        out
    }, {
        ## mini-batch = 3
        batch <- batch.sizes
        batchsum <- cumsum(batch)
        
        list_out <- saffronstar_batch_faster(pval, batch, batchsum, gammai)
        
        alphai <- as.vector(t(list_out$alphai))
        R <- as.vector(t(list_out$R))
        x <- alphai != 0
        
        if (length(x) > 0) {
            alphai <- alphai[x]
            R <- as.numeric(R[x])
        }
        
        batch.no <- rep(seq_len(length(batch)), batch)
        out <- data.frame(pval, batch = batch.no, alphai, R)
    })
}
