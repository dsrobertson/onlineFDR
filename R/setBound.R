#' setBound
#'
#' Sets an upper bound to the number of hypotheses to be tested
#' 
#' Note that when specifying the number of hypotheses to be tested, \code{N} should not be less than the number of hypotheses within your target dataset. 
#' 
#' @param alg A string that takes the value of one of the following: LOND, LORD, LORDdep, SAFFRON, ADDIS, LONDSTAR, LORDSTAR, SAFFRONSTAR, or alphainvesting
#' 
#' @param alpha Overall significance level of the FDR procedure, the default is 0.05.
#' 
#' @param N The number of hypotheses to be tested
#' 
#' @return \item{bound}{A vector of nonnegative numbers that define the update rule to each subsequent significance threshold}
#' 
setBound <- function(alg, alpha = 0.05, N) {
  bound <- switch(alg,
                  LOND = 0.07720838 * alpha * log(pmax(seq_len(N), 2))/(seq_len(N) * exp(sqrt(log(seq_len(N))))),
                  LORD = 0.07720838 * log(pmax(seq_len(N + 1), 2))/(seq_len(N + 1) * 
                                                                      exp(sqrt(log(seq_len(N + 1))))),
                  LORDdep = 0.139307 * alpha/(b0 * seq_len(N) * (log(pmax(seq_len(N), 2)))^3),
                  SAFFRON = 0.4374901658/(seq_len(N)^(1.6)),
                  ADDIS = 0.4374901658/(seq_len(N + 1)^(1.6)),
                  LONDSTAR = 0.07720838 * alpha * log(pmax(seq_len(N), 2))/(seq_len(N) * exp(sqrt(log(seq_len(N))))),
                  LORDSTAR = 0.07720838 * log(pmax(seq_len(N), 2))/(seq_len(N) * exp(sqrt(log(seq_len(N))))),
                  SAFFRONSTAR = 0.4374901658/(seq_len(N + 1)^(1.6)),
                  alphainvesting = 0.4374901658/(seq_len(N)^(1.6))
  )
}
