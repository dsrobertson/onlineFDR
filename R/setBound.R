#' setBound
#'
#' Calculates a default sequence of non-negative numbers \eqn{\gamma_i} that sum
#' to 1 (for LOND and LONDstar, a sequence \eqn{\beta_i} that sums to
#' \eqn{\alpha}), given an upper bound \eqn{N} on the number of hypotheses to be
#' tested.
#' For LORDdep the sequence \eqn{\xi_i} instead satisfies \eqn{\sum_{i=1}^N
#' \xi_i (1 + \log i) = 1}, the condition in Javanmard and Montanari (2018),
#' example 3.7, for \code{LORD(version='dep')}.
#' 
#' @param alg A string that takes the value of one of the following: LOND, LORD,
#'   LORDdep, SAFFRON, ADDIS, LONDstar, LORDstar, SAFFRONstar, Alpha_investing,
#'   Alpha_spending, ADDIS_spending or online_fallback
#'
#' @param alpha Overall significance level of the FDR procedure, the default is
#'   0.05. The bounds for LOND and LONDstar depend on alpha.
#'
#' @param N An upper bound on the number of hypotheses to be tested
#'
#' @return \item{bound}{ A vector giving the values of a default sequence
#' \eqn{\gamma_i} of non-negative numbers.}
#'   
#' @export

setBound <- function(alg, alpha = 0.05, N) {
  
  if(!(alg %in% c('LOND','LORD','LORDdep','SAFFRON','ADDIS','LONDstar',
                  'LORDstar','SAFFRONstar','Alpha_investing', 'Alpha_spending', 'ADDIS_spending', 'online_fallback'))){
    stop('alg must be one of LOND, LORD, LORDdep, SAFFRON, ADDIS, LONDstar, ',
         'LORDstar, SAFFRONstar, Alpha_investing, Alpha_spending, ',
         'ADDIS_spending or online_fallback')
  }
  
  if (alpha <= 0 || alpha > 1) {
    stop("alpha must be between 0 and 1.")
  }
  
  bound <- switch(alg,
                  LOND = rep(alpha/N, N),
                  
                  LORD = (1/sum(log(pmax(seq_len(N),2))/((seq_len(N)) * 
                      exp(sqrt(log(seq_len(N)))))))*log(pmax(seq_len(N),2))/
                      ((seq_len(N))*exp(sqrt(log(seq_len(N))))),
                  
                  LORDdep = (1/sum((1 + log(seq_len(N)))/(seq_len(N) *
                      log(pmax(seq_len(N), 2))^3)))/(seq_len(N) *
                      log(pmax(seq_len(N), 2))^3),
                  
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
