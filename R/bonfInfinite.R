#' Online FDR control based on a Bonferroni-like test
#'
#' This funcion is deprecated, please use \code{\link{AlphaSpending}} instead.
#'
#' Implements online FDR control using a Bonferroni-like test.
#'
#' The function takes as its input either a vector of p-values, or a dataframe
#' with three columns: an identifier (`id'), date (`date') and p-value (`pval').
#' The case where p-values arrive in batches corresponds to multiple instances of
#' the same date. If no column of dates is provided, then the p-values are
#' treated as being ordered sequentially with no batches.
#'
#' The procedure controls FDR for a potentially infinite stream of p-values by
#' using a Bonferroni-like test. Given an overall significance level
#' \eqn{\alpha}, we choose a (potentially infinite) sequence of non-negative
#' numbers \eqn{\alpha_i} such that they sum to \eqn{\alpha}. Hypothesis \eqn{i}
#' is rejected if the \eqn{i}-th p-value is less than or equal to \eqn{\alpha_i}.
#'
#'
#'@param d Either a vector of p-values, or a dataframe with three columns: an
#'  identifier (`id'), date (`date') and p-value (`pval'). If no column of dates
#'  is provided, then the p-values are treated as being ordered sequentially
#'  with no batches.
#'
#'@param alpha Overall significance level of the FDR procedure, the default is
#'  0.05.
#'
#'@param alphai Optional vector of \eqn{\alpha_i}, where hypothesis \eqn{i} is
#'  rejected if the \eqn{i}-th p-value is less than or equal to \eqn{\alpha_i}.
#'  A default is provided as proposed by Javanmard and Montanari (2018),
#'  equation 31.
#'
#'@param random Logical. If \code{TRUE} (the default), then the order of the
#'  p-values in each batch (i.e. those that have exactly the same date) is
#'  randomised.
#'
#'@param date.format Optional string giving the format that is used for dates.
#'
#'
#'@return \item{d.out}{ A dataframe with the original data \code{d} (which
#'  will be reordered if there are batches and \code{random = TRUE}), the
#'  adjusted signifcance thresholds \code{alphai} and the indicator function of
#'  discoveries \code{R}, where \code{R[i] = 1} corresponds to hypothesis
#'  \eqn{i} being rejected (otherwise \code{R[i] = 0}).}
#'
#'
#'@references Javanmard, A. and Montanari, A. (2018) Online Rules for Control of
#'  False Discovery Rate and False Discovery Exceedance. \emph{Annals of
#'  Statistics}, 46(2):526-554.
#'
#'@export

bonfInfinite <- function(d, alpha=0.05, alphai, random=TRUE,
                        date.format="%Y-%m-%d") { # nocov start
    
    .Deprecated("bonfInfinite", package="onlineFDR",
                msg = "bonfInfinite is deprecated, please use AlphaSpending instead.")

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

    if(alpha<0 || alpha>1){
        stop("alpha must be between 0 and 1.")
    }

    if(missing(alphai)){
        alphai <- 0.07720838*alpha*log(pmax(seq_len(N),2)) /
                    ((seq_len(N))*exp(sqrt(log(seq_len(N)))))
    } else if (any(alphai<0)){
        stop("All elements of alphai must be non-negative.")
    } else if(sum(alphai)>alpha){
        stop("The sum of the elements of alphai must not be greater than alpha.")
    }

    R <- as.numeric(pval <= alphai)
    d.out <- data.frame(d, alphai, R)

    return(d.out)
} # nocov end
