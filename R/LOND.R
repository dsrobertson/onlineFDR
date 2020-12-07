#' LOND: Online FDR control based on number of discoveries
#'
#' Implements the LOND algorithm for online FDR control, where LOND stands for
#' (significance) Levels based On Number of Discoveries, as presented by
#' Javanmard and Montanari (2015).
#'
#' The function takes as its input either a vector of p-values, or a dataframe
#' with three columns: an identifier (`id'), date (`date') and p-value (`pval').
#' The case where p-values arrive in batches corresponds to multiple instances
#' of the same date. If no column of dates is provided, then the p-values are
#' treated as being ordered sequentially with no batches.
#'
#' The LOND algorithm controls the FDR for independent p-values (see below for
#' the modification for dependent p-values). Given an overall significance level
#' \eqn{\alpha}, we choose a sequence of non-negative numbers \eqn{\beta_i} such
#' that they sum to \eqn{\alpha}. The values of the adjusted significance
#' thresholds \eqn{\alpha_i} are chosen as follows: \deqn{\alpha_i = (D(i-1) +
#' 1)\beta_i} where \eqn{D(n)} denotes the number of discoveries in the first
#' \eqn{n} hypotheses.
#'
#' A slightly modified version of LOND with thresholds \eqn{\alpha_i =
#' max(D(i-1), 1)\beta_i} provably controls the FDR under positive dependence
#' (PRDS condition), see Zrnic et al. (2018).
#'
#' For arbitrarily dependent p-values, LOND controls the FDR if it is modified
#' with \eqn{\beta_i / H(i)} in place of \eqn{\beta_i}, where \eqn{H(j)} is the
#' i-th harmonic number.
#'
#' Further details of the LOND algorithm can be found in Javanmard and Montanari
#' (2015).
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
#' @param betai Optional vector of \eqn{\beta_i}. A default is provided as
#'   proposed by Javanmard and Montanari (2018), equation 31.
#'
#' @param dep Logical. If \code{TRUE}, runs the modified LOND algorithm which
#'   guarantees FDR control for \emph{dependent} p-values. Defaults to
#'   \code{FALSE}.
#'
#' @param random Logical. If \code{TRUE} (the default), then the order of the
#'   p-values in each batch (i.e. those that have exactly the same date) is
#'   randomised.
#'
#' @param date.format Optional string giving the format that is used for dates.
#'
#' @param original Logical. If \code{TRUE}, runs the original LOND algorithm 
#' of Javanmard and Montanari (2015), otherwise runs the modified algorithm 
#' of Zrnic et al. (2018). Defaults to \code{TRUE}.
#'
#'
#' @return \item{d.out}{ A dataframe with the original data \code{d} (which
#' will be reordered if there are batches and \code{random = TRUE}), the
#' LOND-adjusted significance thresholds \eqn{\alpha_i} and the indicator
#' function of discoveries \code{R}. Hypothesis \eqn{i} is rejected if the
#' \eqn{i}-th p-value is less than or equal to \eqn{\alpha_i}, in which case
#' \code{R[i] = 1}  (otherwise \code{R[i] = 0}).}
#'
#'
#' @references Javanmard, A. and Montanari, A. (2015) On Online Control of False
#' Discovery Rate. \emph{arXiv preprint}, \url{https://arxiv.org/abs/1502.06197}.
#'
#' Javanmard, A. and Montanari, A. (2018) Online Rules for Control of False
#' Discovery Rate and False Discovery Exceedance. \emph{Annals of Statistics},
#' 46(2):526-554.
#'
#' Zrnic, T., Ramdas, A. and Jordan, M.I. (2018). Asynchronous Online Testing of
#' Multiple Hypotheses. \emph{arXiv preprint}, 
#' \url{https://arxiv.org/abs/1812.05068}.
#' 
#'
#' @seealso
#'
#' \code{\link{LONDstar}} presents versions of LORD for \emph{synchronous}
#' p-values, i.e. where each test can only start when the previous test has
#' finished.
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
#' pval = c(2.90e-08, 0.06743, 0.01514, 0.08174, 0.00171,
#'         3.60e-05, 0.79149, 0.27201, 0.28295, 7.59e-08,
#'         0.69274, 0.30443, 0.00136, 0.72342, 0.54757))
#'
#' set.seed(1); LOND(sample.df)
#' 
#' LOND(sample.df, random=FALSE)
#' 
#' set.seed(1); LOND(sample.df, alpha=0.1)
#'
#' @export

LOND <- function(d, alpha = 0.05, betai, dep = FALSE, random = TRUE, date.format = "%Y-%m-%d", 
    original = TRUE) {
    
    checkPval(d)
    
    if (is.data.frame(d)) {
        d <- checkdf(d, random, date.format)
        pval <- d$pval
    } else if (is.vector(d)) {
        pval <- d
    } else {
        stop("d must either be a dataframe or a vector of p-values.")
    }
    
    N <- length(pval)
    
    if (alpha < 0 || alpha > 1) {
        stop("alpha must be between 0 and 1.")
    }
    
    if (missing(betai)) {
        betai <- 0.07720838 * alpha * log(pmax(seq_len(N), 2))/(seq_len(N) * exp(sqrt(log(seq_len(N)))))
    } else if (any(betai < 0)) {
        stop("All elements of betai must be non-negative.")
    } else if (sum(betai) > alpha) {
        stop("The sum of the elements of betai must not be greater than alpha.")
    }
    
    if (dep) {
        den <- cumsum(1/seq_len(N))
        betai <- betai/den
    }
    
    ### Start LOND procedure
    
    out <- lond_faster(pval, betai, original = original)
    out$R <- as.numeric(out$R)
    out
}
TRUE
TRUE
