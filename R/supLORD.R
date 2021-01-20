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
#' treated as being ordered sequentially with no batches.
#'
#' The supLORD procedure provably controls the FDX for p-values that are 
#' conditionally superuniform under the null. supLORD also controls the supFDR 
#' and hence the FDR (even at stopping times). Given an overall significance
#' level \eqn{\alpha}, we choose a sequence of non-negative non-increasing
#' numbers \eqn{\gamma_i} that sum to 1.
#'
#' supLORD requires the user to specify r, a threshold of rejections after which
#' the error control begins to apply, eps, the upper bound on the false 
#' discovery proportion (FDP), and delta, the probability at which the FDP
#' exceeds eps at any time step after making r rejections. As well, the user
#' should specify the variables eta, which controls the pace at which wealth is 
#' spent (as a function of the algorithm's current wealth), and rho, which 
#' controls the length of time before the spending sequence exhausts 
#' the wealth earned from a rejection.
#'
#' Further details of the supLORD procedure can be found in Xu and Ramdas (2021).
#'
#' @param d Either a vector of p-values, or a dataframe with three columns: an
#'   identifier (`id'), date (`date') and p-value (`pval'). If no column of
#'   dates is provided, then the p-values are treated as being ordered
#'   sequentially with no batches.
#'   
#' @param delta The probability at which the FDP exceeds eps (at any time step
#' after making r rejections). Must be between 0 and 1, defaults to 0.05.
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
#' the wealth earned from a rejection. Must be a positive integer.
#' 
#' @param gammai Optional vector of \eqn{\gamma_i}. A default is provided as
#'   proposed by Javanmard and Montanari (2018).
#'
#' @param random Logical. If \code{TRUE} (the default), then the order of the
#'   p-values in each batch (i.e. those that have exactly the same date) is
#'   randomised.
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
#' @references Xu, Z. and Ramdas, A. (2021). Dynamic Algorithms for Online 
#' Multiple Testing. \emph{arXiv preprint}, 
#' \url{https://arxiv.org/abs/2010.13953}.
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

supLORD <- function(d, delta = 0.05, eps, r, eta, rho, gammai, random = TRUE,
                    date.format = "%Y-%m-%d") {
    
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
    
    if (missing(gammai)) {
        gammai <- 0.07720838 * log(pmax(seq_len(N), 2))/(seq_len(N) * 
        exp(sqrt(log(seq_len(N)))))
    } else if (any(gammai < 0)) {
        stop("All elements of gammai must be non-negative.")
    } else if (sum(gammai) > 1) {
        stop("The sum of the elements of gammai must be <= 1.")
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
    
    R <- betai <- alphai <- W <- rep(0, N)
    
    obj <- function(a) {
      (log(1 + log(1/delta)/a) - log(1/delta)/(a+log(1/delta)) - 
         log(1/delta)/(eps*r))^2
    }
    
    a <- optimize(obj, c(0,r/10), tol = 1e-10)$minimum
    beta0 <- ((eps*r/(log(1/delta)/(a*log(1+log(1/delta)/a))))-a)/r
    betai[1] <- beta0

    if (eta > 1){
        gamma_bar0 <- (gammai[1]^eta)/sum(gammai[seq_len(rho)]^eta)
    } else {
        gamma_bar0 <- gammai[1]
    }
    
    alphai[1] <- gamma_bar0*beta0
    R[1] <- (pval[1] <= alphai[1])
    
    if(N == 1){
        out <- data.frame(d, alphai, R)
        return(out)
    }
    
    W[1] <- beta0 + betai[1]*R[1] - alphai[1]
    beta1 <- eps/(log(1/delta)/(a*log(1+log(1/delta)/a)))
    
    for (i in (seq_len(N-1)+1)){
        
        tau <- which(R[seq_len(i-1)] == 1)
        
        if (sum(R) <= r-1){
            betai[i] = beta0
        } else {
            betai[i] = beta1
        }
        
        if (eta > 1 & i <= rho){
            gamma_bar0 <- (gammai[i]^eta)/sum(gammai[seq_len(rho)]^eta)
        } else if(eta > 1 & i > rho){
            gamma_bar0 <- 0
        } else {
            gamma_bar0 <- gammai[i]
        }
    
        if (sum(R) <= 1){
            
            alphai[i] <- gamma_bar0*beta0
            R[i] <- (pval[i] <= alphai[i])

        } else {
            
            m <- length(tau)
            gamma_bar <- rep(0, m)
            
            for(j in seq_len(m)){
                
                cond = eta*W[tau[j]]/beta0
                
                if (cond > 1 & (i - tau[j]) <= rho){
                    gamma_bar[j] <- (gammai[i-tau[j]]^max(cond, 1)) / 
                        sum(gammai[seq_len(rho)]^max(cond, 1))
                } else if (cond > 1 & (i - tau[j]) > rho){
                    gamma_bar[j] <- 0
                } else {
                    gamma_bar[j] <- gammai[i-tau[j]]
                }
            }

            alphai[i] <- beta0*gamma_bar0 + sum(betai[tau]*gamma_bar)
            R[i] = (pval[i] <= alphai[i])
        }
        W[i] <- W[i-1] + betai[i]*R[i] - alphai[i]
    }
    
    out <- data.frame(d, alphai, R)
    return(out)
}