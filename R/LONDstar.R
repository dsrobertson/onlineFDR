#' LONDstar: Asynchronous online mFDR control based on number of discoveries
#'
#' Implements the LOND algorithm for asynchronous online testing, as presented 
#' by Zrnic et al. (2021).
#'
#' The function takes as its input either a vector of p-values, or a dataframe
#' with three columns: an identifier (`id'),
#' p-value (`pval'), or a column describing the conflict sets for the hypotheses. 
#' This takes the form of a vector of decision times or lags. Batch sizes can be 
#' specified as a separate argument (see below).
#'
#' Zrnic et al. (2021) present explicit three versions of LONDstar:
#'
#' 1) \code{version='async'} is for an asynchronous testing
#' process, consisting of tests that start and finish at (potentially) random 
#' times. The discretised finish times of the test correspond to the decision 
#' times. These decision times are given as the input \code{decision.times}
#' for this version of the LONDstar algorithm.
#' 
#' 2) \code{version='dep'} is for online testing under local
#' dependence of the p-values. More precisely, for any \eqn{t>0} we allow the
#' p-value \eqn{p_t} to have arbitrary dependence on the previous \eqn{L_t}
#' p-values. The fixed sequence \eqn{L_t} is referred to as `lags', and is
#' given as the input \code{lags} for this version of the LONDstar algorithm.
#' 
#' 3) \code{version='batch'} is for controlling the mFDR in
#' mini-batch testing, where a mini-batch represents a grouping of tests run
#' asynchronously which result in dependent p-values. Once a mini-batch of tests
#' is fully completed, a new one can start, testing hypotheses independent of
#' the previous batch. The batch sizes are given as the input \code{batch.sizes}
#' for this version of the LONDstar algorithm.
#'
#' Given an overall significance level \eqn{\alpha}, LONDstar requires a
#' sequence of non-negative non-increasing numbers \eqn{\beta_i} that sum to
#' \eqn{\alpha}.
#'
#' Note that these LONDstar algorithms control the \emph{modified} FDR
#' (mFDR). 
#'
#' Further details of the LONDstar algorithms can be found in 
#' Zrnic et al. (2021).
#'
#'
#' @param d Either a vector of p-values, or a dataframe with three columns: an
#'   identifier (`id'), 
#'   p-value (`pval'), and either 
#'   `decision.times', or
#'   `lags', depending on which version you're using. See version for more details.
#'
#' @param alpha Overall significance level of the procedure, the default
#' is 0.05.
#'
#' @param betai Optional vector of \eqn{\beta_i}. A default is provided as
#' proposed by Javanmard and Montanari (2018), equation 31.
#'
#' @param version Takes values 'async', 'dep' or 'batch'. This 
#' specifies the version of LONDstar to use. \code{version='async'} requires a 
#' column of decision times (`decision.times'). \code{version='dep'} requires a
#' column of lags (`lags').
#' \code{version='batch'} requires a vector of batch sizes (`batch.sizes').
#' 
#' @param batch.sizes A vector of batch sizes, this is required for
#'   \code{version='batch'}.
#'   
#' @param display_progress Logical. If \code{TRUE} prints out a progress bar for the algorithm runtime. 
#'
#' @return
#' \item{out}{A dataframe with the original p-values \code{pval}, the
#' adjusted testing levels \eqn{\alpha_i} and the indicator 
#' function of discoveries \code{R}. Hypothesis \eqn{i} is rejected if the
#' \eqn{i}-th p-value is less than or equal to \eqn{\alpha_i}, in which case
#' \code{R[i] = 1}  (otherwise \code{R[i] = 0}).}
#'
#'
#' @references
#' Javanmard, A. and Montanari, A. (2018) Online Rules for Control of False
#' Discovery Rate and False Discovery Exceedance. \emph{Annals of Statistics},
#' 46(2):526-554.
#' 
#' Zrnic, T., Ramdas, A. and Jordan, M.I. (2021). Asynchronous Online Testing of
#' Multiple Hypotheses. \emph{Journal of Machine Learning Research} (to appear),
#' \url{https://arxiv.org/abs/1812.05068}.
#'
#'
#' @seealso
#'
#' \code{\link{LOND}} presents versions of LOND for \emph{synchronous} p-values,
#' i.e. where each test can only start when the previous test has finished.
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
#' LONDstar(sample.df, version='async')

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
#' LONDstar(sample.df2, version='dep')
#' 
#' sample.df3 <- data.frame(
#' id = c('A15432', 'B90969', 'C18705', 'B49731', 'E99902',
#'     'C38292', 'A30619', 'D46627', 'E29198', 'A41418',
#'     'D51456', 'C88669', 'E03673', 'A63155', 'B66033'),
#' pval = c(2.90e-08, 0.06743, 0.01514, 0.08174, 0.00171,
#'         3.60e-05, 0.79149, 0.27201, 0.28295, 7.59e-08,
#'         0.69274, 0.30443, 0.00136, 0.72342, 0.54757))
#' 
#' LONDstar(sample.df3, version='batch', batch.sizes = c(4,6,5))
#'
#' @export

LONDstar <- function(d, alpha = 0.05, version, betai, batch.sizes, display_progress = FALSE) {
    
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
    
    if (alpha < 0 || alpha > 1) {
        stop("alpha must be between 0 and 1.")
    }

    N <- length(pval)
    
    if (missing(betai)) {
        betai <- 0.07720838 * alpha * log(pmax(seq_len(N), 2))/(seq_len(N) * exp(sqrt(log(seq_len(N)))))
    } else if (any(betai < 0)) {

        stop("All elements of betai must be non-negative.")
    } else if (sum(betai) > alpha + .Machine$double.eps * length(betai)) {
        stop("The sum of the elements of betai must be <= alpha.")
    }
    
    version <- checkStarVersion(d, N, version, batch.sizes)
    
    switch(version, {
        
        ## async = 1
        E <- d$decision.times
        out <- londstar_async_faster(pval, 
                                     E, 
                                     betai,
                                     alpha = alpha,
                                     display_progress = display_progress)
        out$R <- as.numeric(out$R)
        out
        
    }, {
        ## dep = 2
        L <- d$lags
        out <- londstar_dep_faster(pval, 
                                   L,
                                   betai,
                                   alpha = alpha,
                                   display_progress = display_progress)
        out$R <- as.numeric(out$R)
        out
     
    }, {
        ## batch = 3
        batch <- batch.sizes
        batchsum <- cumsum(batch)
        
        list_out <- londstar_batch_faster(pval, 
                                          batch, 
                                          batchsum, 
                                          betai,
                                          alpha = alpha,
                                          display_progress = display_progress)
        
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
