#' setBound
#'
#' Calculates a default sequence of non-negative numbers \eqn{\gamma_i} that sum
#' to 1, given an upper bound \eqn{N} on the number of hypotheses to be tested.
#' 
#' @param alg A string that takes the value of one of the following: LOND, LORD,
#'   LORDdep, SAFFRON, ADDIS, LONDstar, LORDstar, SAFFRONstar, or
#'   Alpha_investing
#'
#'
#' @param alpha Overall significance level of the FDR procedure, the default is
#'   0.05. The bounds for LOND and LORDdep depend on alpha.
#'
#' @param N An upper bound on the number of hypotheses to be tested
#'
#' @return \item{bound}{ A vector giving the values of a default sequence
#' \eqn{\gamma_i} of nonnegative numbers.}
#'   
#' @noRd

setBound <- function(alg, alpha = 0.05, N) {
  
  if(!(alg %in% c('LOND','LORD','LORDdep','SAFFRON','ADDIS','LONDstar',
                  'LORDstar','SAFFRONstar','Alpha_investing', 'Alpha_spending', 'ADDIS_spending', 'online_fallback'))){
    stop('alg must be one of LOND, LORD, LORDdep, SAFFRON, ADDIS, LONDstar, LORDstar, SAFFRONstar, Alpha_investing, Alpha_investing, ADDIS_spending or online_fallback')
  }
  
  if (alpha <= 0 || alpha > 1) {
    stop("alpha must be between 0 and 1.")
  }
  
  bound <- switch(alg,
                  LOND = rep(alpha/N, N),
                  
                  LORD = (1/sum(log(pmax(seq_len(N),2))/((seq_len(N)) * 
                      exp(sqrt(log(seq_len(N)))))))*log(pmax(seq_len(N),2))/
                      ((seq_len(N))*exp(sqrt(log(seq_len(N))))),
                  
                  LORDdep = rep(1/N, N),
                  
                  SAFFRON = (1/(seq_len(N))^1.6)/sum(1/(seq_len(N))^1.6),
                  
                  ADDIS = (1/(seq_len(N))^1.6)/sum(1/(seq_len(N))^1.6),
                  
                  LONDstar = (alpha/sum(log(pmax(seq_len(N),2))/((seq_len(N)) * 
                      exp(sqrt(log(seq_len(N)))))))*log(pmax(seq_len(N),2))/
                      ((seq_len(N))*exp(sqrt(log(seq_len(N))))),
                  
                  LORDstar = (1/sum(log(pmax(seq_len(N),2))/((seq_len(N)) * 
                      exp(sqrt(log(seq_len(N)))))))*log(pmax(seq_len(N),2))/
                      ((seq_len(N))*exp(sqrt(log(seq_len(N))))),
                  
                  SAFFRONstar = (1/(seq_len(N))^1.6)/sum(1/(seq_len(N))^1.6),
                  
                  Alpha_investing = (1/(seq_len(N))^1.6)/sum(1/(seq_len(N))^1.6),
                  
                  Alpha_spending = rep(1/N, N),
                  
                  online_fallback = (1/sum(log(pmax(seq_len(N),2))/((seq_len(N)) * 
                          exp(sqrt(log(seq_len(N)))))))*log(pmax(seq_len(N),2))/
                          ((seq_len(N))*exp(sqrt(log(seq_len(N))))),
                  
                  ADDIS_spending = (1/(seq_len(N))^1.6)/sum(1/(seq_len(N))^1.6)
  )
  
  bound
  
}
