#' Alpha-spending for online FWER control
#'
#' Implements online FWER control using a Bonferroni-like test.
#'
#' The function takes as its input either a vector of p-values, or a dataframe
#' with three columns: an identifier (`id'), date (`date') and p-value (`pval').
#' The case where p-values arrive in batches corresponds to multiple instances of
#' the same date. If no column of dates is provided, then the p-values are
#' treated as being ordered sequentially with no batches.
#'
#' Alpha-spending provides strong FWER control for a potentially infinite stream
#' of p-values by using a Bonferroni-like test. Given an overall significance
#' level \eqn{\alpha}, we choose a (potentially infinite) sequence of
#' non-negative numbers \eqn{\gamma_i} such that they sum to 1. Hypothesis
#' \eqn{i} is rejected if the \eqn{i}-th p-value is less than or equal to
#' \eqn{\alpha \gamma_i}.
#' 
#' Note that the procedure controls the generalised familywise error rate
#' (k-FWER) for \eqn{k > 1} if \eqn{\alpha} is replaced by min(\eqn{1,
#' k\alpha}).
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
#'@param gammai Optional vector of \eqn{\gamma_i}, where hypothesis \eqn{i} is
#'  rejected if the \eqn{i}-th p-value is less than or equal to \eqn{\alpha
#'  \gamma_i}. A default is provided as proposed by Javanmard and Montanari
#'  (2018), equation 31.
#'
#'@param random Logical. If \code{TRUE} (the default), then the order of the
#'  p-values in each batch (i.e. those that have exactly the same date) is
#'  randomised.
#'
#'@param date.format Optional string giving the format that is used for dates.
#'
#'
#'@return \item{d.out}{ A dataframe with the original data \code{d} (which will
#'  be reordered if there are batches and \code{random = TRUE}), the adjusted
#'  signifcance thresholds \code{alphai} and the indicator function of
#'  discoveries \code{R}, where \code{R[i] = 1} corresponds to hypothesis
#'  \eqn{i} being rejected (otherwise \code{R[i] = 0}).}
#'
#'
#'@references Javanmard, A. and Montanari, A. (2018) Online Rules for Control of
#'  False Discovery Rate and False Discovery Exceedance. \emph{Annals of
#'  Statistics}, 46(2):526-554.
#'
#' Tian, J. and Ramdas, A. (2019b). Online control of the familywise error rate.
#' \emph{arXiv preprint}, \url{https://arxiv.org/abs/1910.04900}.
#'
#'
#' @examples
#' sample.df <- data.frame(
#' id = c('A15432', 'B90969', 'C18705', 'B49731', 'E99902',
#'     'C38292', 'A30619', 'D46627', 'E29198', 'A41418',
#'     'D51456', 'C88669', 'E03673', 'A63155', 'B66033'),
#' date = as.Date(c(rep('2014-12-01',3),
#'                 rep('2015-09-21',5),
#'                 rep('2016-05-19',2),
#'                 '2016-11-12',
#'                 rep('2017-03-27',4))),
#' pval = c(2.90e-17, 0.06743, 0.01514, 0.08174, 0.00171,
#'         3.60e-05, 0.79149, 0.27201, 0.28295, 7.59e-08,
#'         0.69274, 0.30443, 0.00136, 0.72342, 0.54757))
#'
#' set.seed(1); Alpha_spending(sample.df)
#' 
#' Alpha_spending(sample.df, random=FALSE)
#' 
#' set.seed(1); Alpha_spending(sample.df, alpha=0.1)
#'
#'
#'@export

Alpha_spending <- function(d, alpha = 0.05, gammai, random = TRUE, date.format = "%Y-%m-%d") {
    
    if (is.data.frame(d)) {
        d <- checkdf(d, random, date.format)
        pval <- d$pval
    } else if (is.vector(d)) {
        pval <- d
    } else {
        stop("d must either be a dataframe or a vector of p-values.")
    }
    
    checkPval(pval)
    N <- length(pval)
    
    if (alpha < 0 || alpha > 1) {
        stop("alpha must be between 0 and 1.")
    }
    
    if (missing(gammai)) {
        gammai <- 0.07720838 * log(pmax(seq_len(N), 2))/((seq_len(N)) * exp(sqrt(log(seq_len(N)))))
    } else if (any(gammai < 0)) {
        stop("All elements of gammai must be non-negative.")
    } else if (sum(gammai) > 1) {
        stop("The sum of the elements of gammai must not be greater than 1.")
    }
    
    R <- as.numeric(pval <= alpha * gammai)
    d.out <- data.frame(d, alpha * gammai, R)
    
    return(d.out)
}
TRUE
TRUE
