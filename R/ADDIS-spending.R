#' ADDIS-spending: Adaptive discarding algorithm for online FWER control
#'
#' Implements the ADDIS algorithm for online FWER control, where ADDIS stands
#' for an ADaptive algorithm that DIScards conservative nulls, as presented by
#' Tian and Ramdas (2019b). The procedure compensates for the power loss of
#' Alpha-spending, by including both adaptivity in the fraction of null
#' hypotheses and the conservativeness of nulls.
#'
#' The function takes as its input either a vector of p-values, or a dataframe
#' with three columns: an identifier (`id'),
#' p-value (`pval'), and lags, if the dependent version is specified (see below). 
#' Given an overall significance level \eqn{\alpha}, ADDIS depends on constants \eqn{\lambda} and
#' \eqn{\tau}, where \eqn{\lambda < \tau}. Here \eqn{\tau \in (0,1)} represents
#' the threshold for a hypothesis to be selected for testing: p-values greater
#' than \eqn{\tau} are implicitly `discarded' by the procedure, while
#' \eqn{\lambda \in (0,1)} sets the threshold for a p-value to be a candidate
#' for rejection: ADDIS-spending will never reject a p-value larger than
#' \eqn{\lambda}. The algorithms also require a sequence of non-negative
#' non-increasing numbers \eqn{\gamma_i} that sum to 1.
#'
#' The ADDIS-spending procedure provably controls the FWER in the strong sense
#' for independent p-values. Note that the procedure also controls the
#' generalised familywise error rate (k-FWER) for \eqn{k > 1} if \eqn{\alpha} is
#' replaced by min(\eqn{1, k\alpha}).
#'
#' Tian and Ramdas (2019b) also presented a version for handling local
#' dependence. More precisely, for any \eqn{t>0} we allow the p-value \eqn{p_t}
#' to have arbitrary dependence on the previous \eqn{L_t} p-values. The fixed
#' sequence \eqn{L_t} is referred to as `lags', and is given as the input
#' \code{lags} for this version of the ADDIS-spending algorithm.
#'
#' Further details of the ADDIS-spending algorithms can be found in Tian and
#' Ramdas (2019b).
#'
#' @param d Either a vector of p-values, or a dataframe with three columns: an
#'   identifier (`id'), 
#'   p-value (`pval'), and lags 
#'   (`lags').
#'
#' @param alpha Overall significance level of the procedure, the default is
#'   0.05.
#'
#' @param gammai Optional vector of \eqn{\gamma_i}. A default is provided with
#'   \eqn{\gamma_j} proportional to \eqn{1/j^(1.6)}.
#'
#' @param lambda Optional parameter that sets the threshold for `candidate'
#'   hypotheses. Must be between 0 and 1, defaults to 0.25.
#'
#' @param tau Optional threshold for hypotheses to be selected for testing. Must
#'   be between 0 and 1, defaults to 0.5.
#'
#' @param dep Logical. If \code{TRUE} runs the version for locally dependent 
#' p-values
#'
#'
#' @return \item{d.out}{A dataframe with the original p-values \code{pval}, the
#'   adjusted testing levels \eqn{\alpha_i} and the indicator function of
#'   discoveries \code{R}. Hypothesis \eqn{i} is rejected if the \eqn{i}-th
#'   p-value is less than or equal to \eqn{\alpha_i}, in which case \code{R[i] =
#'   1}  (otherwise \code{R[i] = 0}).}
#'
#'
#' @references Tian, J. and Ramdas, A. (2019b). Online control of the familywise error rate.
#' \emph{arXiv preprint}, \url{https://arxiv.org/abs/1910.04900}.
#'
#'
#' @seealso
#'
#' \code{\link{ADDIS}} provides online control of the FDR.
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
#' lags = rep(1,15))
#' 
#' ADDIS_spending(sample.df) #independent
#' 
#' ADDIS_spending(sample.df, dep = TRUE) #Locally dependent
#'
#' @export

ADDIS_spending <- function(d, alpha = 0.05, gammai, lambda = 0.25, tau = 0.5, dep = FALSE) {
    
    d <- checkPval(d)
    
    if (is.data.frame(d)) {
        pval <- d$pval
    } else if (is.vector(d)) {
        pval <- d
    } else {
        stop("d must either be a dataframe or a vector of p-values.")
    }
    
    if (alpha <= 0 || alpha > 1) {
        stop("alpha must be between 0 and 1.")
    }
    
    if (lambda <= 0 || lambda > 1) {
        stop("lambda must be between 0 and 1.")
    }
    
    if (tau <= 0 || tau > 1) {
        stop("tau must be between 0 and 1.")
    }
    
    if (lambda >= tau) {
        stop("lambda must be less than tau.")
    }
    
    N <- length(pval)
    
    if (missing(gammai)) {
        gammai <- 0.4374901658/(seq_len(N)^(1.6))
    } else if (any(gammai < 0)) {
        stop("All elements of gammai must be non-negative.")
    } else if (sum(gammai) > 1) {
        stop("The sum of the elements of gammai must not be greater than 1.")
    }
    
    
    if (!(dep)) {
        
        alphai <- R <- rep(0, N)
        
        alphai[1] <- alpha * (tau - lambda) * gammai[1]
        R[1] <- (pval[1] <= alphai[1])
        
        if (N == 1) {
            d.out <- data.frame(pval, alphai, R)
            return(d.out)
        }
        
        select.sum <- (pval[1] <= tau)
        cand.sum <- (pval[1] <= lambda)
        
        for (i in (seq_len(N - 1) + 1)) {
            
            alphai[i] <- alpha * (tau - lambda) * gammai[1 + select.sum - cand.sum]
            R[i] <- (pval[i] <= alphai[i])
            
            select.sum <- select.sum + (pval[i] <= tau)
            cand.sum <- cand.sum + (pval[i] <= lambda)
            
        }
    } else {
        
        checkStarVersion(d, N, "dep")
        
        L <- d$lags
        
        R <- select <- cand <- alphai <- rep(0, N)
        
        alphai[1] <- alpha * (tau - lambda) * gammai[1]
        R[1] <- pval[1] <= alphai[1]
        
        if (N == 1) {
            d.out <- data.frame(pval, lag = L, alphai, R)
            return(d.out)
        }
        
        select[1] <- (pval[1] <= tau)
        cand[1] <- (pval[1] <= lambda)
        
        for (i in (seq_len(N - 1) + 1)) {
            
            alphai[i] <- alpha * (tau - lambda) * gammai[1 + min(L[i], i - 1) + sum(select[seq_len(max(0, 
                i - L[i] - 1))]) - sum(cand[seq_len(max(0, i - L[i] - 1))])]
            
            R[i] <- (pval[i] <= alphai[i])
            
            select[i] <- (pval[i] <= tau)
            cand[i] <- (pval[i] <= lambda)
        }
    }
    
    d.out <- data.frame(pval, alphai, R)
    
    return(d.out)
}
TRUE
TRUE
