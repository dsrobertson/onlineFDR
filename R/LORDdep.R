#' LORD (dep): Online FDR control based on recent discovery for dependent
#' p-values
#'
#' This funcion is deprecated, please use \code{\link{LORD}} instead with
#' \code{version = 'dep'}.
#'
#' LORDdep implements the LORD procedure for online FDR control for dependent
#' p-values, where LORD stands for (significance) Levels based On Recent
#' Discovery, as presented by Javanmard and Montanari (2018).
#'
#' The function takes as its input either a vector of p-values or a dataframe
#' with three columns: an identifier (`id'), date (`date') and p-value (`pval').
#' The case where p-values arrive in batches corresponds to multiple instances
#' of the same date. If no column of dates is provided, then the p-values are
#' treated as being ordered sequentially with no batches.
#'
#' This modified LORD procedure controls FDR for dependent p-values. Given an
#' overall significance level \eqn{\alpha}, we choose a sequence of non-negative
#' numbers \eqn{\xi_i} such that they satisfy a condition given in Javanmard and
#' Montanari (2018), example 3.8.
#'
#' The procedure depends on constants \eqn{w_0} and \eqn{b_0}, where \eqn{w_0
#' \ge 0} represents the intial `wealth' and \eqn{b_0 > 0} represents the
#' `payout' for rejecting a hypothesis. We require \eqn{w_0+b_0 \le \alpha} for
#' FDR control to hold.
#'
#' Further details of the modified LORD procedure can be found in Javanmard and
#' Montanari (2018).
#'
#'
#' @param d Either a vector of p-values, or a dataframe with three columns: an
#'   identifier (`id'), date (`date') and p-value (`pval'). If no column of
#'   dates is provided, then the p-values are treated as being ordered
#'   sequentially with no batches.
#'
#' @param alpha Overall significance level of the FDR procedure, the default is
#'   0.05.
#'
#' @param xi Optional vector of \eqn{\xi_i}. A default is provided to satisfy
#'   the condition given in Javanmard and Montanari (2018), example 3.7.
#'
#' @param w0 Initial `wealth' of the procedure. Defaults to \eqn{\alpha/10}.
#' 
#' @param b0 The `payout' for rejecting a hypothesis. Defaults to \eqn{\alpha -
#'   w_0}.
#'
#' @param random Logical. If \code{TRUE} (the default), then the order of the
#'   p-values in each batch (i.e. those that have exactly the same date) is
#'   randomised.
#'
#' @param date.format Optional string giving the format that is used for dates.
#'
#'
#' @return \item{d.out}{ A dataframe with the original data \code{d} (which
#' will be reordered if there are batches and \code{random = TRUE}), the
#' LORD-adjusted significance thresholds \eqn{\alpha_i} and the indicator
#' function of discoveries \code{R}. Hypothesis \eqn{i} is rejected if the
#' \eqn{i}-th p-value is less than or equal to \eqn{\alpha_i}, in which case
#' \code{R[i] = 1}  (otherwise \code{R[i] = 0}).}
#'
#'
#' @references Javanmard, A. and Montanari, A. (2018) Online Rules for Control
#' of False Discovery Rate and False Discovery Exceedance. \emph{Annals of
#' Statistics}, 46(2):526-554.
#'
#'
#' @export

LORDdep <- function(d, alpha=0.05, xi, w0=alpha/10, b0=alpha - w0, random=TRUE,
                    date.format="%Y-%m-%d") {
    
    .Deprecated("LORD", package="onlineFDR",
    msg = "LORDdep is deprecated. Use LORD instead with version='dep'.")

    if(is.data.frame(d)){
        checkdf(d, random, date.format)
        pval <- d$pval
    } else if(is.vector(d)){
        pval <- d
    } else {
        stop("d must either be a dataframe or a vector of p-values.")
    }
    
    checkPval(pval)
    N <- length(pval)

    if(alpha<=0 || alpha>1){
        stop("alpha must be between 0 and 1.")
    }

    if(missing(xi)){
        xi <- 0.139307*alpha/(b0*seq_len(N)*(log(pmax(seq_len(N),2)))^3)
    } else if (any(xi<0)){
        stop("All elements of xi must be non-negative.")
    } else if(sum(xi)>alpha/b0){
        stop("The sum of the elements of xi must not be greater than alpha/b0.")
    }

    if(w0 < 0){
        stop("w0 must be non-negative.")
    } else if (w0 > b0){
        stop("w0 must be less than b0.")
    } else if(b0 <= 0){
        stop("b0 must be positive.")
    } else if(w0+b0 > alpha & !(isTRUE(all.equal(w0+b0, alpha)))){
        stop("The sum of w0 and b0 must not be greater than alpha.")
    }

    alphai <- rep(0, N)
    R <- W <- rep(0, N+1)

    R[1] <- 1
    W[1] <- w0

    alphai[1] <- phi <- xi[1]*w0
    R[2] <- pval[1] <= alphai[1]
    W[2] <- w0 - phi + R[2]*b0

    if(N == 1){
        R <- R[2]
        d.out <- data.frame(d, alphai, R)
        return(d.out)
    }

    for (i in (seq_len(N-1)+1)){
        tau <- max(which(R[seq_len(i)] == 1))
        alphai[i] <- phi <- xi[i]*W[tau]

        R[i+1] <- pval[i] <= alphai[i]
        W[i+1] <- W[i] - phi + R[i+1]*b0
    }

    R <- R[(seq_len(N)+1)]
    d.out <- data.frame(d, alphai, R)

    return(d.out)
}
