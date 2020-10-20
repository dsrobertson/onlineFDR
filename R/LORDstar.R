#' LORDstar: Asynchronous online mFDR control based on recent discovery
#'
#' Implements LORD algorithms for asynchronous online testing, as presented by
#' Zrnic et al. (2018).
#'
#' The function takes as its input either a vector of p-values, or a dataframe
#' with three columns: an identifier (`id'),
#' p-value (`pval'), or a column describing the conflict sets for the hypotheses. 
#' This takes the form of a vector of decision times or lags. Batch sizes can be 
#' specified as a separate argument (see below).
#'
#' Zrnic et al. (2018) present explicit three versions of LORDstar:
#'
#' 1) \code{version='async'} is for an asynchronous testing process, consisting
#' of tests that start and finish at (potentially) random times. The discretised
#' finish times of the test correspond to the decision times. These decision
#' times are given as the input \code{decision.times} for this version of the
#' LORDstar algorithm.
#'
#' 2) \code{version='dep'} is for online testing under local dependence of the
#' p-values. More precisely, for any \eqn{t>0} we allow the p-value \eqn{p_t} to
#' have arbitrary dependence on the previous \eqn{L_t} p-values. The fixed
#' sequence \eqn{L_t} is referred to as `lags', and is given as the input
#' \code{lags} for this version of the LORDstar algorithm.
#'
#' 3) \code{version='batch'} is for controlling the mFDR in mini-batch testing,
#' where a mini-batch represents a grouping of tests run asynchronously which
#' result in dependent p-values. Once a mini-batch of tests is fully completed,
#' a new one can start, testing hypotheses independent of the previous batch.
#' The batch sizes are given as the input \code{batch.sizes} for this version of
#' the LORDstar algorithm.
#'
#' Given an overall significance level \eqn{\alpha}, LORDstar depends on
#' \eqn{w_0} (where \eqn{0 \le w_0 \le \alpha}), which represents the intial
#' `wealth' of the procedure. The algorithms also require a sequence of
#' non-negative non-increasing numbers \eqn{\gamma_i} that sum to 1.
#'
#' Note that these LORDstar algorithms control the \emph{modified} FDR (mFDR).
#' The `async' version also controls the usual FDR if the p-values are assumed
#' to be independent.
#'
#' Further details of the LORDstar algorithms can be found in Zrnic et al.
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
#' @param gammai Optional vector of \eqn{\gamma_i}. A default is provided as
#'   proposed by Javanmard and Montanari (2018), equation 31.
#'
#' @param version Takes values 'async', 'dep' or 'batch'. This specifies the
#'   version of LORDstar to use. \code{version='async'} requires a 
#' column of decision times (`decision.times'). \code{version='dep'} requires a
#' column of lags (`lags').
#' \code{version='batch'} requires a vector of batch sizes (`batch.sizes').
#'
#' @param w0 Initial `wealth' of the procedure, defaults to \eqn{\alpha/10}.
#'
#' @param batch.sizes A vector of batch sizes, this is required for
#'   \code{version='batch'}.
#'
#'
#' @return \item{d.out}{A dataframe with the original p-values \code{pval}, the
#' adjusted testing levels \eqn{\alpha_i} and the indicator function of
#' discoveries \code{R}. Hypothesis \eqn{i} is rejected if the \eqn{i}-th
#' p-value is less than or equal to \eqn{\alpha_i}, in which case \code{R[i] =
#' 1}  (otherwise \code{R[i] = 0}).}
#'
#'
#' @references Javanmard, A. and Montanari, A. (2018) Online Rules for Control
#' of False Discovery Rate and False Discovery Exceedance. \emph{Annals of
#' Statistics}, 46(2):526-554.
#'
#' Zrnic, T., Ramdas, A. and Jordan, M.I. (2018). Asynchronous Online Testing of
#' Multiple Hypotheses. \emph{arXiv preprint},
#' \url{https://arxiv.org/abs/1812.05068}.
#'
#'
#' @seealso
#'
#' \code{\link{LORD}} presents versions of LORD for \emph{synchronous} p-values,
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
#' LORDstar(sample.df, version='async')
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
#' LORDstar(sample.df2, version='dep')
#' 
#' sample.df3 <- data.frame(
#' #' id = c('A15432', 'B90969', 'C18705', 'B49731', 'E99902',
#'     'C38292', 'A30619', 'D46627', 'E29198', 'A41418',
#'     'D51456', 'C88669', 'E03673', 'A63155', 'B66033'),
#' pval = c(2.90e-08, 0.06743, 0.01514, 0.08174, 0.00171,
#'         3.60e-05, 0.79149, 0.27201, 0.28295, 7.59e-08,
#'         0.69274, 0.30443, 0.00136, 0.72342, 0.54757))
#' 
#' LORDstar(sample.df3, version='batch', batch.sizes = c(4,6,5))
#'
#' @export

LORDstar <- function(d, alpha=0.05, version, gammai, w0, batch.sizes) {

    if(is.data.frame(d)){
        checkSTARdf(d, version)
        pval <- d$pval
    } else if(is.vector(d)){
        pval <- d
    } else {
        stop("d must either be a dataframe or a vector of p-values.")
    }
    
    if(alpha<=0 || alpha>1){
        stop("alpha must be between 0 and 1.")
    }
    
    if(missing(w0)){
        w0 = alpha/10
    } else if(w0 < 0){
        stop("w0 must be non-negative.")
    } else if(w0 > alpha) {
        stop("w0 must not be greater than alpha.")
    }
    
    checkPval(pval)
    N <- length(pval)
    
    if(missing(gammai)){
        gammai <- 0.07720838*log(pmax(seq_len(N),2)) /
            (seq_len(N)*exp(sqrt(log(seq_len(N)))))
    } else if (any(gammai<0)){
        stop("All elements of gammai must be non-negative.")
    } else if(sum(gammai)>1){
        stop("The sum of the elements of gammai must be <= 1.")
    }
    
    version <- checkStarVersion(d, N, version, batch.sizes)
    
    switch(version,
           ## async = 1
           {
               E <- d$decision.times
               alphai <- R <- rep(0, N)
               
               alphai[1] <- gammai[1]*w0
               R[1] <- pval[1] <= alphai[1]
               
               if(N == 1){
                   d.out <- data.frame(pval, alphai, R)
                   return(d.out)
               }
               
               for (i in (seq_len(N-1)+1)){
                   
                   r <- which(R[seq_len(i-1)] == 1 & E[seq_len(i-1)] <= i-1)
                   
                   if(length(r) <= 1){
                       alphai[i] <- gammai[i]*w0 + (alpha - w0)*sum(gammai[i-r])
                       R[i] <- pval[i] <= alphai[i]
                       
                   } else {
                       alphai[i] <- gammai[i]*w0 + (alpha - w0)*gammai[i-r[1]] +
                           alpha*sum(gammai[i-r[-1]])
                       
                       R[i] = pval[i] <= alphai[i]
                   }
               }
               
               d.out <- data.frame(pval, alphai, R)
           },
           ## dep = 2
           {
               L <- d$lags
               
               alphai <- R <- rep(0, N)
               
               alphai[1] <- gammai[1]*w0
               R[1] <- pval[1] <= alphai[1]
               
               if(N == 1){
                   d.out <- data.frame(pval, lag=L, alphai, R)
                   return(d.out)
               }
               
               for (i in (seq_len(N-1)+1)){
                   
                   r <- which(R[seq_len(i-1)] == 1 & seq_len(i-1) <= i-1-L[i])
                   
                   if(length(r) <= 1){
                       alphai[i] <- gammai[i]*w0 + (alpha - w0)*sum(gammai[i-r])
                       R[i] <- pval[i] <= alphai[i]
                       
                   } else {
                       alphai[i] <- gammai[i]*w0 + (alpha - w0)*gammai[i-r[1]] +
                           alpha*sum(gammai[i - r[-1]])
                       
                       R[i] = pval[i] <= alphai[i]
                   }
               }
               
               d.out <- data.frame(pval, lag=L, alphai, R)
           },
           ## mini-batch = 3
           {
               n <- batch.sizes
               ncum <- cumsum(n)
               
               alphai <- matrix(NA, nrow = length(n), ncol = max(n))
               R <- matrix(0, nrow = length(n), ncol = max(n))
               
               for (i in seq_len(n[1])){ 
                   alphai[1,i] <- gammai[i]*w0
                   R[1,i] <- pval[i] <= alphai[1,i]
               }
               
               if(length(n) == 1){
                   d.out <- data.frame(pval, batch = rep(1, n),
                                       alphai = as.vector(t(alphai)),
                                       R = as.vector(t(R)))
                   
                   return(d.out)
                   
               } else {
                   
                   r <- integer(0)
                   
                   for (b in seq_len(length(n)-1)+1){
                       
                       Rcum <- cumsum(rowSums(R))
                       
                       for (i in seq_len(n[b])){ 
                           
                           if(max(Rcum)>0){
                               r <- sapply(seq_len(max(Rcum)),
                                           function(x){match(1,Rcum>=x)})
                           }
                           
                           if(length(r) <= 1){
                               alphai[b,i] <- gammai[ncum[b-1]+i]*w0 + 
                                   (alpha - w0)*sum(gammai[ncum[b-1] + i - ncum[r]])
                               
                               R[b,i] <- pval[ncum[b-1]+i] <= alphai[b,i]
                               
                           } else {
                               alphai[b,i] <- gammai[ncum[b-1]+i]*w0 +
                                   (alpha - w0)*gammai[ncum[b-1]+ i - ncum[r[1]]] +
                                   alpha*sum(gammai[ncum[b-1] + i - ncum[r[-1]]])
                               
                               R[b,i] = pval[ncum[b-1]+i] <= alphai[b,i]
                           }
                       }
                   }
                   
                   alphai <- as.vector(t(alphai))
                   R <- as.vector(t(R))
                   x <- !(is.na(alphai))
                   
                   if(length(x) > 0){
                       alphai <- alphai[x]
                       R <- R[x]
                   }
                   
                   batch.no <- rep(seq_len(length(n)),n)
                   d.out <- data.frame(pval, batch=batch.no, alphai, R)
               }
           })
    
    return(d.out)
}
