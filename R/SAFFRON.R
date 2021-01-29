#' SAFFRON: Adaptive online FDR control
#'
#' Implements the SAFFRON procedure for online FDR control, where SAFFRON stands
#' for Serial estimate of the Alpha Fraction that is Futilely Rationed On true
#' Null hypotheses, as presented by Ramdas et al. (2018). The algorithm is based
#' on an estimate of the proportion of true null hypotheses. More precisely,
#' SAFFRON sets the adjusted test levels based on an estimate of the amount of
#' alpha-wealth that is allocated to testing the true null hypotheses.
#'
#' The function takes as its input either a vector of p-values or a dataframe
#' with three columns: an identifier (`id'), date (`date') and p-value (`pval').
#' The case where p-values arrive in batches corresponds to multiple instances
#' of the same date. If no column of dates is provided, then the p-values are
#' treated as being ordered sequentially with no batches.
#'
#' SAFFRON procedure provably controls FDR for independent p-values. Given an
#' overall significance level \eqn{\alpha}, we choose a sequence of non-negative
#' non-increasing numbers \eqn{\gamma_i} that sum to 1.
#'
#' SAFFRON depends on constants \eqn{w_0} and \eqn{\lambda}, where \eqn{w_0}
#' satisfies \eqn{0 \le w_0 \le \alpha} and represents the initial `wealth' of
#' the procedure, and \eqn{0 < \lambda < 1} represents the threshold for a
#' `candidate' hypothesis. A `candidate' refers to p-values smaller than
#' \eqn{\lambda}, since SAFFRON will never reject a p-value larger than
#' \eqn{\lambda}.
#'
#' Note that FDR control also holds for the SAFFRON procedure if only the
#' p-values corresponding to true nulls are mutually independent, and
#' independent from the non-null p-values.
#'
#' The SAFFRON procedure can lose power in the presence of conservative nulls,
#' which can be compensated for by adaptively `discarding' these p-values. This
#' option is called by setting \code{discard=TRUE}, which is the same algorithm
#' as ADDIS.
#'
#' Further details of the SAFFRON procedure can be found in Ramdas et al.
#' (2018).
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
#' @param lambda Optional threshold for a `candidate' hypothesis, must be
#'   between 0 and 1. Defaults to 0.5.
#'
#' @param discard Logical. If \code{TRUE} then runs the ADDIS algorithm with
#'   adaptive discarding of conservative nulls. The default is \code{FALSE}.
#'
#' @param tau.discard Optional threshold for hypotheses to be selected for
#'   testing. Must be between 0 and 1, defaults to 0.5. This is required if
#'   \code{discard=TRUE}.
#'
#' @param random Logical. If \code{TRUE} (the default), then the order of the
#'   p-values in each batch (i.e. those that have exactly the same date) is
#'   randomised.
#'
#' @param date.format Optional string giving the format that is used for dates.
#'
#'
#' @return \item{out}{ A dataframe with the original data \code{d} (which
#'   will be reordered if there are batches and \code{random = TRUE}), the
#'   LORD-adjusted significance thresholds \eqn{\alpha_i} and the indicator
#'   function of discoveries \code{R}. Hypothesis \eqn{i} is rejected if the
#'   \eqn{i}-th p-value is less than or equal to \eqn{\alpha_i}, in which case
#'   \code{R[i] = 1}  (otherwise \code{R[i] = 0}).}
#'
#'
#' @references Ramdas, A., Zrnic, T., Wainwright M.J. and Jordan, M.I. (2018).
#'   SAFFRON: an adaptive algorithm for online control of the false discovery
#'   rate. \emph{Proceedings of the 35th International Conference in Machine
#'   Learning}, 80:4286-4294.
#'
#' @seealso
#'
#' \code{\link{SAFFRONstar}} presents versions of SAFFRON for
#' \emph{asynchronous} testing, i.e. where each hypothesis test can itself be a
#' sequential process and the tests can overlap in time.
#'
#' If option \code{discard=TRUE}, SAFFRON is the same as \code{\link{ADDIS}}.
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
#' SAFFRON(sample.df, random=FALSE)
#' 
#' set.seed(1); SAFFRON(sample.df)
#' 
#' set.seed(1); SAFFRON(sample.df, alpha=0.1, w0=0.025)
#'
#' SAFFRON(sample.df, discard=TRUE, random=FALSE)
#'
#' @export

SAFFRON <- function(d, alpha = 0.05, gammai, w0 = 0.025, lambda = 0.5, random = TRUE, date.format = "%Y-%m-%d", 
    discard = FALSE, tau.discard = 0.5) {
    
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
    
    if (lambda <= 0 || lambda > 1) {
        stop("lambda must be between 0 and 1.")
    }
    
    if (discard == TRUE) {
        return(ADDIS(pval, alpha, gammai, w0, lambda, tau.discard))
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
    
    ### Start SAFFRON algorithm
    out <- saffron_faster(pval, gammai)
    out$R <- as.numeric(out$R)
    out

}
