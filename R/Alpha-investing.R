#' Alpha-investing for online FDR control
#'
#' Implements a variant of the Alpha-investing algorithm of Foster and Stine
#' (2008) that guarantees FDR control, as proposed by Ramdas et al. (2018). This
#' procedure uses SAFFRON's update rule with the constant \eqn{\lambda} replaced
#' by a sequence \eqn{\lambda_i = \alpha_i}. This is also equivalent to using
#' the ADDIS algorithm with \eqn{\tau = 1} and \eqn{\lambda_i = \alpha_i}.
#'
#' The function takes as its input either a vector of p-values or a dataframe
#' with three columns: an identifier (`id'), date (`date') and p-value (`pval').
#' The case where p-values arrive in batches corresponds to multiple instances
#' of the same date. If no column of dates is provided, then the p-values are
#' treated as being ordered sequentially with no batches.
#'
#' The Alpha-investing procedure provably controls FDR for independent p-values.
#' Given an overall significance level \eqn{\alpha}, we choose a sequence of
#' non-negative non-increasing numbers \eqn{\gamma_i} that sum to 1.
#' Alpha-investing depends on a constant \eqn{w_0}, which satisfies \eqn{0 \le
#' w_0 \le \alpha} and represents the initial `wealth' of the procedure.
#' 
#' Further details of the Alpha-investing procedure and its modification can be
#' found in Foster and Stine (2008) and Ramdas et al. (2018).
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
#' @param gammai Optional vector of \eqn{\gamma_i}. A default is provided with
#'   \eqn{\gamma_j} proportional to \eqn{1/j^(1.6)}.
#'
#' @param w0 Initial `wealth' of the procedure, defaults to \eqn{\alpha/2}. Must
#'   be between 0 and \eqn{\alpha}.
#'
#' @param random Logical. If \code{TRUE} (the default), then the order of the
#'   p-values in each batch (i.e. those that have exactly the same date) is
#'   randomised.
#'
#' @param date.format Optional string giving the format that is used for dates.
#'
#'
#' @return \item{out}{ A dataframe with the original data \code{d} (which will
#'   be reordered if there are batches and \code{random = TRUE}), the
#'   LORD-adjusted significance thresholds \eqn{\alpha_i} and the indicator
#'   function of discoveries \code{R}. Hypothesis \eqn{i} is rejected if the
#'   \eqn{i}-th p-value is less than or equal to \eqn{\alpha_i}, in which case
#'   \code{R[i] = 1}  (otherwise \code{R[i] = 0}).}
#'
#'
#' @references Foster, D. and Stine R. (2008). \eqn{\alpha}-investing: a
#'   procedure for sequential control of expected false discoveries.
#'   \emph{Journal of the Royal Statistical Society (Series B)}, 29(4):429-444.
#'
#'   Ramdas, A., Zrnic, T., Wainwright M.J. and Jordan, M.I. (2018). SAFFRON: an
#'   adaptive algorithm for online control of the false discovery rate.
#'   \emph{Proceedings of the 35th International Conference in Machine
#'   Learning}, 80:4286-4294.
#'
#' @seealso
#'
#' \code{\link{SAFFRON}} uses the update rule of Alpha-investing but with
#' constant \eqn{\lambda}.
#'
#'
#' @examples
#' sample.df <- data.frame(
#' id = c('A15432', 'B90969', 'C18705', 'B49731', 'E99902',
#'     'C38292', 'A30619', 'D46627', 'E29198', 'A41418',
#'     'D51456', 'C88669', 'E03673', 'A63155', 'B66033'),
#' date = as.Date(c(rep('2014-12-01',3),
#'                rep('2015-09-21',5),
#'                 rep('2016-05-19',2),
#'                 '2016-11-12',
#'                rep('2017-03-27',4))),
#' pval = c(2.90e-08, 0.06743, 0.01514, 0.08174, 0.00171,
#'         3.60e-05, 0.79149, 0.27201, 0.28295, 7.59e-08,
#'         0.69274, 0.30443, 0.00136, 0.72342, 0.54757))
#'
#' Alpha_investing(sample.df, random=FALSE)
#' 
#' set.seed(1); Alpha_investing(sample.df)
#' 
#' set.seed(1); Alpha_investing(sample.df, alpha=0.1, w0=0.025)
#' 
#' @export

Alpha_investing <- function(d, alpha = 0.05, gammai, w0, random = TRUE, date.format = "%Y-%m-%d") {
    
    d <- checkPval(d)
    
    if (is.data.frame(d)) {
        d <- checkdf(d, random, date.format)
        pval <- d$pval
    } else if (is.vector(d)) {
        pval <- d
    } else {
        stop("d must either be a dataframe or a vector of p-values.")
    }
    
    N <- length(pval)
    
    if (alpha <= 0 || alpha > 1) {
        stop("alpha must be between 0 and 1.")
    }
    
    if (missing(gammai)) {
        gammai <- 0.4374901658/(seq_len(N)^(1.6))
    } else if (any(gammai < 0)) {
        stop("All elements of gammai must be non-negative.")
    } else if (sum(gammai) > 1) {
        stop("The sum of the elements of gammai must not be greater than 1.")
    }
    
    if (missing(w0)) {
        w0 = alpha/2
    } else if (w0 < 0) {
        stop("w0 must be non-negative.")
    } else if (w0 >= alpha) {
        stop("w0 must be less than alpha.")
    }
    
    ### Start algorithm
    out <- alphainvesting_faster(pval)
    out$R <- as.numeric(out$R)
    out
}
TRUE
TRUE
