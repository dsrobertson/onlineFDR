#' LORD: Online FDR control based on recent discovery
#'
#' Implements the LORD procedure for online FDR control, where LORD stands for
#' (significance) Levels based On Recent Discovery, as presented by Javanmard
#' and Montanari (2018) and Ramdas et al. (2017).
#'
#' The function takes as its input either a vector of p-values or a dataframe
#' with three columns: an identifier (`id'), date (`date') and p-value (`pval').
#' The case where p-values arrive in batches corresponds to multiple instances
#' of the same date. If no column of dates is provided, then the p-values are
#' treated as being ordered in sequence, arriving one at a time.
#'
#' LORD++ and the discarding version of LORD provably control the FDR for
#' independent p-values; for LORD 3 this has only been shown empirically (see
#' below for dependent p-values). Given an overall significance level
#' \eqn{\alpha}, we
#' choose a sequence of non-negative non-increasing numbers \eqn{\gamma_i} that
#' sum to 1.
#'
#' Javanmard and Montanari (2018) presented versions of LORD which differ in the
#' way the adjusted significance thresholds \eqn{\alpha_i} are calculated. The
#' significance thresholds for LORD 2 are based on all previous discovery times.
#' LORD 2 has been superseded by the algorithm given in Ramdas et al. (2017),
#' LORD++ (\code{version='++'}), which is the default version. The significance
#' thresholds for LORD 3 (\code{version=3}) are based on the time of the last
#' discovery as well as the 'wealth' accumulated at that time. Finally, Tian and
#' Ramdas (2019) presented a version of LORD (\code{version='discard'}) that can
#' improve the power of the procedure in the presence of conservative nulls by
#' adaptively `discarding' these p-values.
#'
#' LORD depends on constants \eqn{w_0} and (for versions 3 and 'dep') \eqn{b_0},
#' where \eqn{0 \le w_0 \le \alpha} represents the initial `wealth' of the
#' procedure and \eqn{b_0 > 0} represents the `payout' for rejecting a
#' hypothesis. We require \eqn{w_0+b_0 \le \alpha} for FDR control to hold.
#' Version 'discard' also depends on a constant \eqn{\tau}, where \eqn{\tau \in
#' (0,1]} represents the threshold for a hypothesis to be selected for testing:
#' p-values greater than \eqn{\tau} are implicitly `discarded' by the procedure.
#'
#' Note that FDR control also holds for the LORD procedure if only the p-values
#' corresponding to true nulls are mutually independent, and independent from
#' the non-null p-values.
#'
#' FDR control is not guaranteed for positively dependent p-values in general,
#' even with the default settings. For example, the FDR of LORD++ can exceed
#' \eqn{\alpha} when long runs of positively correlated null p-values follow
#' strong signals; see \code{vignette("theory")}.
#'
#' For dependent p-values, a modified LORD procedure was proposed in Javanmard
#' and Montanari (2018), which is called by setting \code{version='dep'}. Given
#' an overall significance level \eqn{\alpha}, we choose a sequence of
#' non-negative numbers \eqn{\xi_i} (supplied through the argument
#' \code{gammai}) such that \eqn{\sum_i \xi_i (1 + \log i) \le \alpha/b_0} if
#' \eqn{w_0 \le b_0}, which is the condition in Javanmard and Montanari (2018),
#' example 3.7. If \eqn{w_0 > b_0}, the condition is instead \eqn{\sum_i (w_0 +
#' b_0 \log i) \xi_i \le \alpha}, which follows from their Theorem 3.6.
#'
#' Further details of the LORD procedures can be found in Javanmard and
#' Montanari (2018), Ramdas et al. (2017) and Tian and Ramdas (2019).
#'
#'
#' @param d Either a vector of p-values, or a dataframe with three columns: an
#'   identifier (`id'), date (`date') and p-value (`pval'). If no column of
#'   dates is provided, then the p-values are treated as being ordered
#'   in sequence, arriving one at a time.
#'
#' @param alpha Overall significance level of the FDR procedure, the default is
#'   0.05.
#'
#' @param gammai Optional vector of \eqn{\gamma_i}. A default is provided as
#'   proposed by Javanmard and Montanari (2018), equation 31 for all versions of
#'   LORD except 'dep'. The latter is provided a default to satisfy a condition
#'   given in Javanmard and Montanari (2018), example 3.7.
#'
#' @param version Takes values '++', 3, 'discard', or 'dep'. This specifies the
#'   version of LORD to use, and defaults to '++'.
#'
#' @param w0 Initial `wealth' of the procedure, defaults to \eqn{\alpha/10}
#'   (\eqn{\min(\alpha/10, \tau\alpha/2)} for \code{version='discard'}).
#'   Version 'discard' requires \eqn{w_0 \le \tau\alpha}.
#'
#' @param b0 The `payout' for rejecting a hypothesis in versions 3 and 'dep' of
#'   LORD. Defaults to \eqn{\alpha - w_0}.
#'
#' @param tau.discard Optional threshold for hypotheses to be selected for
#'   testing. Must be in (0, 1], defaults to 0.5. This is only used if
#'   \code{version='discard'}.
#'
#' @param random Logical. If \code{TRUE} (the default), then the order of the
#'   p-values in each batch (i.e. those that have exactly the same date) is
#'   randomised.
#'   
#' @param display_progress Logical. If \code{TRUE} prints out a progress bar for the algorithm runtime. 
#'
#' @param date.format Optional string giving the format that is used for dates.
#'
#'
#' @return \item{d.out}{ A dataframe with the p-values \code{pval} (reordered if
#'   there are batches and \code{random = TRUE}), the LORD-adjusted
#'   significance thresholds \eqn{\alpha_i}, the indicator function of
#'   discoveries \code{R} and, if \code{d} is a dataframe with a column
#'   \code{id}, the identifiers \code{id}. Hypothesis \eqn{i} is rejected if
#'   the \eqn{i}-th p-value is less than or equal to \eqn{\alpha_i}, in which
#'   case \code{R[i] = 1} (otherwise \code{R[i] = 0}).}
#'
#'
#' @references Javanmard, A. and Montanari, A. (2018). Online Rules for Control
#' of False Discovery Rate and False Discovery Exceedance. \emph{Annals of
#' Statistics}, 46(2):526-554.
#'
#' Ramdas, A., Yang, F., Wainwright, M.J. and Jordan, M.I. (2017). Online control
#' of the false discovery rate with decaying memory. \emph{Advances in Neural
#' Information Processing Systems 30}, 5650-5659.
#'
#' Tian, J. and Ramdas, A. (2019). ADDIS: an adaptive discarding algorithm for
#' online FDR control with conservative nulls.
#' \emph{Advances in Neural Information Processing Systems 32}, 9388-9396.
#'
#'
#' @seealso
#'
#' \code{\link{LORDstar}} presents versions of LORD for \emph{asynchronous}
#' testing, i.e. where each hypothesis test can itself be a sequential process
#' and the tests can overlap in time.
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
#' pval = c(2.90e-08, 0.06743, 0.01514, 0.08174, 0.00171,
#'         3.60e-05, 0.79149, 0.27201, 0.28295, 7.59e-08,
#'         0.69274, 0.30443, 0.00136, 0.72342, 0.54757))
#'
#' LORD(sample.df, random=FALSE)
#' 
#' set.seed(1); LORD(sample.df, version='dep')
#' 
#' set.seed(1); LORD(sample.df, version='discard')
#' 
#' set.seed(1); LORD(sample.df, alpha=0.1, w0=0.05)
#' 
#'
#' @export

LORD <- function(d, alpha = 0.05, gammai, version = "++", w0, b0, tau.discard = 0.5, 
    random = TRUE, display_progress = FALSE, date.format = "%Y-%m-%d") {
    
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
    
    if (version %in% c(1, 2)) {
        stop("LORD 1 and LORD 2 have been superseded by LORD++, please use ",
             "version '++' instead.")
    } else if (!(version %in% c("++", 3, "3", "discard", "dep"))) {
        stop("version must be '++', 3, 'discard' or 'dep'.")
    }
    
    if (version == "discard") {
        if (missing(tau.discard)) {
            tau.discard = 0.5
        } else if (tau.discard <= 0 || tau.discard > 1) {
            stop("tau.discard must be between 0 and 1.")
        }
    }

    if (missing(w0)) {
        if (version == "discard") {
            w0 = min(alpha/10, tau.discard * alpha/2)
        } else {
            w0 = alpha/10
        }
    } else if (w0 < 0) {
        stop("w0 must be non-negative.")
    }

    if (version %in% c(3, "3", "dep")) {
        if (missing(b0)) {
            if (w0 >= alpha) {
                stop("w0 must be less than alpha.")
            }
            b0 = alpha - w0
        } else if (b0 <= 0) {
            stop("b0 must be positive.")
        } else if (w0 + b0 > alpha & !(isTRUE(all.equal(w0 + b0, alpha)))) {
            stop("The sum of w0 and b0 must not be greater than alpha.")
        }
    } else {
        if (missing(b0)) {
            b0 = alpha - w0
        }
        if (w0 > alpha) {
            stop("w0 must not be greater than alpha.")
        }
        if (version == "discard" &&
            w0 - tau.discard * alpha > sqrt(.Machine$double.eps)) {
            stop("w0 must not be greater than tau.discard * alpha.")
        }
    }

    if (version != "dep") {
        if (missing(gammai)) {
            gammai <- 0.07908197 * log(pmax(seq_len(N + 1), 2))/
                (seq_len(N + 1) *
                exp(sqrt(log(seq_len(N + 1)))))
        } else if (any(gammai < 0)) {
            stop("All elements of gammai must be non-negative.")
        } else if (sumExceeds(gammai, 1)) {
            stop("The sum of the elements of gammai must be <= 1.")
        } else if (!isNonIncreasing(gammai)) {
            stop("The elements of gammai must be non-increasing.")
        } else if (length(gammai) < N) {
            stop("gammai must have length at least equal to the number of ",
                 "p-values.")
        }
    } else if (w0 <= b0) {
        if (missing(gammai)) {
            gammai <- 0.1393066 /(seq_len(N) * (log(pmax(seq_len(N), 2)))^3)
        } else if (any(gammai < 0)) {
            stop("All elements of gammai must be non-negative.")
        } else if (sumExceeds(gammai * (1 + log(seq_along(gammai))),
                              alpha/b0)) {
            stop("The sum of the elements of gammai*(1 + log(i)) must be <= ",
                 "alpha/b0.")
        } else if (length(gammai) < N) {
            stop("gammai must have length at least equal to the number of ",
                 "p-values.")
        }
    } else {
        if (missing(gammai)) {
            # Scaled so that sum_i (w0 + b0*log(i))*gammai[i] = alpha over all
            # i >= 1, so the sequence does not depend on N. 5.068667247 and
            # 2.109742802 are the sums over i >= 1 of 1/(i*log(max(i,2))^3) and
            # log(i)/(i*log(max(i,2))^3) (rounded up)
            gammai <- alpha/(5.068667247*w0 + 2.109742802*b0) /
                (seq_len(N) * (log(pmax(seq_len(N), 2)))^3)
        } else if (any(gammai < 0)) {
            stop("All elements of gammai must be non-negative.")
        } else if (sumExceeds((w0 + b0*log(seq_along(gammai)))*gammai, alpha)) {
            stop("The sum of the elements of (w0 + b0*log(seq_len(N)))*gammai must be <= alpha.")
        } else if (length(gammai) < N) {
            stop("gammai must have length at least equal to the number of ",
                 "p-values.")
        }
    }
    
    if (version == "++") {
        version <- 1
    } else if (version == "discard") {
        version <- 2
    } else if (version == "3" || version == 3) {
        version <- 3
    } else if (version == "dep") {
        version <- 4
    }
    
    out <- lord_faster(pval, 
                       gammai,
                       version, 
                       alpha = alpha,
                       w0 = w0,
                       b0 = b0,
                       taudiscard = tau.discard,
                       display_progress = display_progress)
    out$R <- as.numeric(out$R)
    if(is.data.frame(d) && !is.null(d$id)) {
        out$id <- d$id
    }
    out
}
