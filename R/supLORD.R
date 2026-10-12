#' supLORD: Online control of the false discovery exceedance (FDX) and the 
#' FDR at stopping times
#'
#' Implements the supLORD procedure, which controls both FDX and FDR, including
#' the FDR at stopping times, as presented by Xu and Ramdas (2021).
#'
#' The function takes as its input either a vector of p-values or a dataframe
#' with three columns: an identifier (`id'), date (`date') and p-value (`pval').
#' The case where p-values arrive in batches corresponds to multiple instances
#' of the same date. If no column of dates is provided, then the p-values are
#' treated as being ordered in sequence, arriving one at a time.
#'
#' The supLORD procedure provably controls the FDX for p-values that are 
#' conditionally superuniform under the null. supLORD also controls the supFDR 
#' and hence the FDR (even at stopping times). Given an overall significance
#' level \eqn{\alpha}, we choose a sequence of non-negative non-increasing
#' numbers \eqn{\gamma_i} that sum to 1.
#'
#' supLORD requires the user to specify r, a threshold of rejections after which
#' the error control begins to apply, eps, the upper bound on the false 
#' discovery proportion (FDP), and delta, an upper bound on the probability
#' that the FDP exceeds eps at any time step after making r rejections. As
#' well, the user
#' should specify the variables eta, which controls the pace at which wealth is 
#' spent (as a function of the algorithm's current wealth), and rho, which 
#' controls the length of time before the spending sequence exhausts 
#' the wealth earned from a rejection.
#'
#' The test levels follow the dynamic alpha schedule of Xu and Ramdas (2021),
#' equations 16 and 17, with the default boost sequence of their Appendix A,
#' equation 21: the initial wealth and each of the first r - 1 rejections earn
#' \eqn{\beta_0 = b - a/r}, and the r-th and later rejections earn
#' \eqn{b = \epsilon a \log(1 + L/a)/L}, where \eqn{L = \log(1/\delta)}. The
#' offset \eqn{a} is the canonical solution of their equation 20.
#'
#' Further details of the supLORD procedure can be found in Xu and Ramdas (2021).
#'
#' @param d Either a vector of p-values, or a dataframe with three columns: an
#'   identifier (`id'), date (`date') and p-value (`pval'). If no column of
#'   dates is provided, then the p-values are treated as being ordered
#'   in sequence, arriving one at a time.
#'   
#' @param delta Upper bound on the probability that the FDP exceeds eps (at any
#' time step after making r rejections). Must be between 0 and 1, defaults to
#' 0.05.
#'
#' @param eps The upper bound on the FDP. Must be between 0 and 1.
#'
#' @param r The threshold of rejections after which the error control 
#' begins to apply. Must be a positive integer.
#' 
#' @param eta Controls the pace at which wealth is spent as a function of the 
#' algorithm's current wealth. Must be a positive real number.
#' 
#' @param rho Controls the length of time before the spending sequence exhausts 
#' the wealth earned from a rejection. Must be a positive integer, no greater
#' than the length of a user-supplied \code{gammai}.
#' 
#' @param gammai Optional vector of \eqn{\gamma_i}. A default is provided as
#'   proposed by Javanmard and Montanari (2018), equation 31.
#'
#' @param random Logical. If \code{TRUE} (the default), then the order of the
#'   p-values in each batch (i.e. those that have exactly the same date) is
#'   randomised.
#'   
#' @param display_progress Logical. If \code{TRUE} prints out a progress bar for
#'  the algorithm runtime. 
#'
#' @param date.format Optional string giving the format that is used for dates.
#'
#' @return \item{d.out}{ A dataframe with the original data \code{d} (which
#'   will be reordered if there are batches and \code{random = TRUE}), the
#'   supLORD-adjusted significance thresholds \eqn{\alpha_i} and the indicator
#'   function of discoveries \code{R}. Hypothesis \eqn{i} is rejected if the
#'   \eqn{i}-th p-value is less than or equal to \eqn{\alpha_i}, in which case
#'   \code{R[i] = 1}  (otherwise \code{R[i] = 0}).}
#'
#'
#' @references Javanmard, A. and Montanari, A. (2018). Online Rules for Control
#' of False Discovery Rate and False Discovery Exceedance. \emph{Annals of
#' Statistics},
#' 46(2):526-554.
#'
#' Xu, Z. and Ramdas, A. (2021). Dynamic Algorithms for Online 
#' Multiple Testing. \emph{Annual Conference on Mathematical and Scientific 
#' Machine Learning}, PMLR, 145:955-986.
#'
#'
#' @examples
#' 
#' set.seed(1)
#' N <- 1000
#' B <- rbinom(N, 1, 0.5)
#' Z <- rnorm(N, mean = 3*B)
#' pval <- pnorm(-Z)
#'
#' out <- supLORD(pval, eps=0.15, r=30, eta=0.05, rho=30, random=FALSE)
#' head(out)
#' sum(out$R)
#'
#' @export

supLORD <- function(d, delta = 0.05, eps, r, eta, rho, gammai, random = TRUE, display_progress = FALSE, date.format = "%Y-%m-%d") {
    
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

    if (!missing(gammai)) {
        if (any(gammai < 0)) {
            stop("All elements of gammai must be non-negative.")
        } else if (sumExceeds(gammai, 1)) {
            stop("The sum of the elements of gammai must be <= 1.")
        }
    }

    if (delta <=0 | delta >= 1){
        stop("delta must be between 0 and 1.")
    } else if (eps <=0 | eps >= 1) {
        stop("eps must be between 0 and 1.")
    } else if (r %% 1 != 0 | r < 1) {
        stop("r must be a positive integer")
    } else if (eta <= 0){
        stop("eta must be a positive real number.")
    } else if (rho %% 1 != 0 | rho < 1){
        stop("rho must be a positive integer.")
    }

    if (missing(gammai)) {
        M <- max(N, rho)
        gammai <- 0.07908197 * log(pmax(seq_len(M), 2))/(seq_len(M) *
        exp(sqrt(log(seq_len(M)))))
    } else if (length(gammai) < N) {
        stop("gammai must have length at least equal to the number of ",
             "p-values.")
    } else if (length(gammai) < rho) {
        stop("rho must not be greater than the length of gammai.")
    }

    R <- alphai <- rep(0, N)

    # Offset a: the canonical solution of Xu and Ramdas (2021), equation 20,
    # found on the log scale. The bracket is valid because the left-hand side
    # minus the right-hand side is decreasing in a, is at least 1 at the lower
    # end and is negative at the upper end.
    L <- log(1/delta)
    f <- function(la) {
        log(exp(la) + L) - la - L/(exp(la) + L) - L/(eps*r)
    }
    la <- stats::uniroot(f, c(log(L) - 2 - L/(eps*r), log(sqrt(eps*r*L/2)) + 1),
                         tol = 1e-12)$root
    a <- exp(la)

    # Boosts, equation 21: beta0 for the initial wealth and the first r - 1
    # rejections, beta1 for the r-th and later rejections
    logbar <- L/(a*log(1 + L/a))
    beta0 <- (eps*r/logbar - a)/r
    beta1 <- eps/logbar

    # Spending sequence gamma-bar of equation 16, computed on the log scale so
    # that large exponents do not underflow
    lg <- log(gammai[seq_len(rho)])
    lgmax <- max(lg)
    gamma_bar_active <- function(dj, e, lnorm) {
        ifelse(dj <= rho,
               exp(e*(log(gammai[pmin(dj, rho)]) - lgmax) - lnorm), 0)
    }
    lognorm <- function(e) {
        log(sum(exp(e*(lg - lgmax))))
    }

    # Spending sequence for the initial wealth, gamma-bar_i(W(0), W(0))
    if (eta > 1 && is.finite(lgmax)) {
        gamma_bar0 <- gamma_bar_active(seq_len(N), eta, lognorm(eta))
    } else if (eta > 1) {
        gamma_bar0 <- rep(0, N)
    } else {
        gamma_bar0 <- gammai[seq_len(N)]
    }

    # For each rejection: its time, boost, exponent eta*W/W(0) and, if this is
    # greater than 1, the log normalising constant of gamma-bar, all of which
    # are fixed once the rejection is made
    tj <- bj <- ej <- lnj <- rep(0, N)
    m <- 0
    W <- beta0

    if (display_progress) {
        pb <- progress::progress_bar$new(
            format = "  Computing [:bar] :percent eta: :eta",
            total = N, clear = FALSE, width = 60)
    }

    for (i in seq_len(N)) {

        if (display_progress) {
            pb$tick()
        }

        # Dynamic alpha schedule, equation 17
        alphai[i] <- gamma_bar0[i]*beta0
        if (m > 0) {
            j <- seq_len(m)
            dj <- i - tj[j]
            gamma_bar <- gammai[dj]
            act <- ej[j] > 1
            if (any(act)) {
                gamma_bar[act] <- if (is.finite(lgmax)) {
                    gamma_bar_active(dj[act], ej[j][act], lnj[j][act])
                } else 0
            }
            alphai[i] <- alphai[i] + sum(bj[j]*gamma_bar)
        }

        R[i] <- (pval[i] <= alphai[i])

        if (R[i] == 1) {
            m <- m + 1
            if (m <= r - 1) {
                bj[m] <- beta0
            } else {
                bj[m] <- beta1
            }
            # wealth, equation 10
            W <- W + bj[m] - alphai[i]
            tj[m] <- i
            ej[m] <- eta*W/beta0
            if (ej[m] > 1 && is.finite(lgmax)) {
                lnj[m] <- lognorm(ej[m])
            }
        } else {
            W <- W - alphai[i]
        }
    }

    out <- data.frame(d, alphai, R)
    out
}
