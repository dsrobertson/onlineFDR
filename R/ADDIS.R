#' ADDIS: Adaptive discarding algorithm for online FDR control
#'
#' Implements the ADDIS algorithm for online FDR control, where ADDIS stands for
#' an ADaptive algorithm that DIScards conservative nulls, as presented by Tian
#' and Ramdas (2019). The algorithm compensates for the power loss of SAFFRON
#' with conservative nulls, by including both adaptivity in the fraction of
#' null hypotheses (like SAFFRON) and the conservativeness of nulls (unlike
#' SAFFRON).
#'
#' The function takes as its input a vector of p-values. Given an overall
#' significance level \eqn{\alpha}, ADDIS depends on constants \eqn{w_0}
#' \eqn{\lambda} and \eqn{\tau}. \eqn{w_0} represents the intial `wealth' of the
#' procedure and satisfies \eqn{0 \le w_0 \le \tau \lambda \alpha}. \eqn{\tau
#' \in (0,1)} represents the threshold for a hypothesis to be selected for
#' testing: p-values greater than \eqn{\tau} are implicitly `discarded' by the
#' procedure. Finally, \eqn{\lambda \in (0,1)} sets the threshold for a p-value
#' to be a candidate for rejection: ADDIS will never reject a p-value larger
#' than \eqn{\tau \lambda}. The algorithm also require a sequence of
#' non-negative non-increasing numbers \eqn{\gamma_i} that sum to 1.
#'
#' The ADDIS procedure provably controls the FDR for independent p-values. Tian
#' and Ramdas (2019) also presented a version for an asynchoronous testing
#' process, consisting of tests that start and finish at (potentially) random
#' times. The discretised finish times of the test correspond to the decision
#' times. These decision times are given as the input \code{decision.times}.
#' 
#' Further details of the ADDIS algorithms can be found in Tian and Ramdas
#' (2019).
#'
#' @param pval A vector of p-values.
#'
#' @param alpha Overall significance level of the procedure, the default is
#'   0.05.
#'
#' @param gammai Optional vector of \eqn{\gamma_i}. A default is provided with
#'   \eqn{\gamma_j} proportional to \eqn{1/j^(1.6)}.
#'
#' @param w0 Initial `wealth' of the procedure, defaults to \eqn{\tau \lambda
#'   \alpha/2}.
#'
#' @param lambda Optional parameter that sets the threshold for `candidate'
#'   hypotheses. Must be between 0 and 1, defaults to 0.5.
#'
#' @param tau Optional threshold for hypotheses to be selected for testing. Must
#'   be between 0 and 1, defaults to 0.5.
#'
#' @param async Logical. If \code{TRUE} runs the version for an asynchronous
#'   testing process
#'
#' @param decision.times A vector of decision times for the hypothesis tests,
#'   this is required if \code{async=TRUE}.
#'
#'
#' @return \item{d.out}{A dataframe with the original p-values \code{pval}, the
#'   adjusted testing levels \eqn{\alpha_i} and the indicator function of
#'   discoveries \code{R}. Hypothesis \eqn{i} is rejected if the \eqn{i}-th
#'   p-value is less than or equal to \eqn{\alpha_i}, in which case \code{R[i] =
#'   1}  (otherwise \code{R[i] = 0}).}
#'
#'
#' @references Tian, J. and Ramdas, A. (2019). ADDIS: an adaptive discarding
#'   algorithm for online FDR control with conservative nulls. \emph{arXiv
#'   preprint}, \url{https://arxiv.org/abs/1905.11465}.
#'
#'
#' @seealso
#'
#' ADDIS is identical to \code{\link{SAFFRON}} with option \code{discard=TRUE}.
#'
#' ADDIS with option \code{async=TRUE} is identical to \code{\link{SAFFRONstar}}
#' with option \code{discard=TRUE}.
#'
#'
#' @examples
#' pval = c(2.90e-08, 0.06743, 3.51e-04, 0.0154, 0.04723,
#'         3.60e-05, 0.79149, 0.27201, 0.28295, 7.59e-06,
#'         0.69274, 0.30443, 0.00136, 0.82342, 0.54757)
#'
#' ADDIS(pval)
#'
#' ADDIS(pval, async=TRUE, decision.times=seq_len(15)) # Same as above
#' 
#' ADDIS(pval, async=TRUE, decision.times=seq_len(15)+1) # Asynchronous
#'
#' @export

ADDIS <- function(pval, alpha=0.05, gammai, w0, lambda=0.5, tau=0.5,
                  async=FALSE, decision.times) {
    
    if(alpha<=0 || alpha>1){
        stop("alpha must be between 0 and 1.")
    }
    
    if(lambda<=0 || lambda>1){
        stop("lambda must be between 0 and 1.")
    }
    
    if(tau<=0 || tau>1){
        stop("tau must be between 0 and 1.")
    }
    
    checkPval(pval)
    N <- length(pval)
    
    if(missing(gammai)){
        gammai <- 0.4374901658/(seq_len(N+1)^(1.6))
    } else if (any(gammai<0)){
        stop("All elements of gammai must be non-negative.")
    } else if(sum(gammai)>1){
        stop("The sum of the elements of gammai must not be greater than 1.")
    }
    
    if(missing(w0)){
        w0 = tau*lambda*alpha/2
    } else if(w0 < 0){
        stop("w0 must be non-negative.")
    } else if(w0 > tau*lambda*alpha){
        stop("w0 must be less than tau*lambda*alpha")
    }
    
    if(!(async)){
        
        alphai <- R <- cand <- Cj.plus <- rep(0, N)
        
        alphai[1] <- w0*gammai[1]
        R[1] <- (pval[1] <= alphai[1])
        
        if(N == 1){
            d.out <- data.frame(pval, alphai, R)
            return(d.out)
        }
        
        cand.sum <- 0
        selected <- (pval <= tau)
        S <- cumsum(selected)
        
        for (i in (seq_len(N-1)+1)){
            
            kappai <- which(R[seq_len(i-1)] == 1)
            K <- length(kappai)
            
            cand[i-1] <- (pval[i-1] <= tau*lambda)
            cand.sum <- cand.sum + cand[i-1]
            
            if (K > 1) {
                
                kappai.star <- sapply(kappai,
                                      function(x){sum(selected[seq_len(x)])})
                
                Kseq <- seq_len(K-1)
                
                Cj.plus[Kseq] <- Cj.plus[Kseq] + cand[i-1]
                Cj.plus.sum <- sum(gammai[S[i-1]-
                                          kappai.star[Kseq]-Cj.plus[Kseq]+1])
                
                Cj.plus[K] <- sum(cand[seq(from=kappai[K]+1,
                                           to=max(i-1, kappai+1))])
                
                Cj.plus.sum <- Cj.plus.sum + 
                    gammai[S[i-1]-kappai.star[K]-Cj.plus[K]+1] - 
                    gammai[S[i-1]-kappai.star[1]-Cj.plus[1]+1]
                
                alphai.hat <- w0*gammai[S[i-1]-cand.sum+1] + 
                    (tau*(1-lambda)*alpha - w0)*
                    gammai[S[i-1]-kappai.star[1]-Cj.plus[1]+1] + 
                    tau*(1-lambda)*alpha*Cj.plus.sum
                
                alphai[i] <- min(tau*lambda, alphai.hat)
                R[i] <- (pval[i] <= alphai[i])
                
            } else if(K == 1){
                
                kappai.star <- sum(selected[seq_len(kappai)])
                
                Cj.plus[1] <- sum(cand[seq(from=kappai+1,
                                           to=max(i-1, kappai+1))])
                
                alphai.hat <- w0*gammai[S[i-1]-cand.sum+1] + 
                    (tau*(1-lambda)*alpha - w0)*
                    gammai[S[i-1]-kappai.star-Cj.plus[1]+1]
                
                alphai[i] <- min(tau*lambda, alphai.hat)
                R[i] <- (pval[i] <= alphai[i])
                
            } else {
                
                alphai.hat <- w0*gammai[S[i-1]-cand.sum+1]
                alphai[i] <- min(tau*lambda, alphai.hat)
                R[i] <- (pval[i] <= alphai[i])
            }
        }
    } else {
        
        if(length(decision.times) != length(pval)){
            stop("Please provide a decision time for each p-value.")
        }
        
        E <- decision.times
        
        alphai <- R <- S <- cand <- Cj.plus <- rep(0, N)
        
        selected <- (pval <= tau)
        
        alphai[1] <- w0*gammai[1]
        R[1] <- (pval[1] <= alphai[1])
        
        if(N == 1){
            d.out <- data.frame(pval, alphai, R)
            return(d.out)
        }
        
        for (i in (seq_len(N-1)+1)){
            
            kappai <- which(R[seq_len(i-1)]==1 & E[seq_len(i-1)] <= i-1)
            K <- length(kappai)
            
            cand[i-1] <- (pval[i-1] <= tau*lambda)
            cand.sum <- sum(cand[seq_len(i-1)] & E[seq_len(i-1)] <= i-1)
            
            S[i-1] <- sum((selected[seq_len(i-1)] & E[seq_len(i-1)] <= i-1) +
                              (E[seq_len(i-1)] >= i))
            
            if (K > 1) {
                
                kappai.star <- sapply(kappai,
                                      function(x){sum(selected[seq_len(x)])})
                
                Kseq <- seq_len(K)
                
                Cj.plus[Kseq] <- sapply(Kseq,
                                 function(x){sum(cand[seq(from=kappai[x]+1,
                                              to=max(i-1,kappai[x]+1))] &
                                             E[seq(from=kappai[x]+1,
                                               to=max(i-1,kappai[x]+1))]<=i-1)})
                
                Cj.plus.sum <- sum(gammai[S[i-1]-
                                          kappai.star[Kseq]-Cj.plus[Kseq]+1]) -
                    gammai[S[i-1]-kappai.star[1]-Cj.plus[1]+1]
                
                alphai.tilde <- w0*gammai[S[i-1]-cand.sum+1] + 
                    (tau*(1-lambda)*alpha - w0)*
                    gammai[S[i-1]-kappai.star[1]-Cj.plus[1]+1] + 
                    tau*(1-lambda)*alpha*Cj.plus.sum
                
                alphai[i] <- min(tau*lambda, alphai.tilde)
                R[i] <- (pval[i] <= alphai[i])
                
            } else if(K == 1){
                
                kappai.star <- sum(selected[seq_len(kappai)])
                
                Cj.plus[1] <- sum(cand[seq(from=kappai+1,
                                           to=max(i-1,kappai+1))] &
                                      E[seq(from=kappai+1,
                                            to=max(i-1,kappai+1))] <= i-1)
                
                alphai.tilde <- w0*gammai[S[i-1]-cand.sum+1] + 
                    (tau*(1-lambda)*alpha - w0)*
                    gammai[S[i-1]-kappai.star-Cj.plus[1]+1]
                
                alphai[i] <- min(tau*lambda, alphai.tilde)
                R[i] <- (pval[i] <= alphai[i])
                
            } else {
                
                alphai.tilde <- w0*gammai[S[i-1]-cand.sum+1]
                alphai[i] <- min(tau*lambda, alphai.tilde)
                R[i] <- (pval[i] <= alphai[i])
            }
        }
    }

    d.out <- data.frame(pval, alphai, R)
    
    return(d.out)
}
