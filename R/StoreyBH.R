#' StoreyBH: Offline FDR control using the St-BH procedure
#'
#' Implements the Storey-BH algorithm for offline FDR control, as presented by
#' Storey (2002).
#'
#' The function takes as its input either a vector of p-values, or a dataframe
#' with a column of p-values (`pval').
#'
#' @param d Either a vector of p-values, or a dataframe with the column: p-value
#'   (`pval').
#'
#' @param alpha Overall significance level of the FDR procedure, the default is
#'   0.05.
#'
#' @param lambda Threshold for Storey-BH, must be between 0 and 1. Defaults to
#'   0.5.
#'
#' @return \item{ordered_d}{ A dataframe with the original data \code{d} and the
#'   indicator function of discoveries \code{R}. Hypothesis \eqn{i} is rejected
#'   if the \eqn{i}-th p-value is less than or equal to \eqn{(r/n)\alpha}, where
#'   \eqn{r} is the rank of the \eqn{i}-th p-value within an ordered set and
#'   \eqn{n} is the total number of hypotheses. If hypothesis \eqn{i} is
#'   rejected, \code{R[i] = 1} (otherwise \code{R[i] = 0}).}
#'
#' @references Storey, J.D. (2002). A direct approach to false discovery rates.
#'   \emph{J. R. Statist. Soc. B}: 64, Part 3, 479-498.
#'
#' @examples
#'
#' pvals <- c(2.90e-08, 0.06743, 0.01514, 0.08174, 0.00171,
#'         3.60e-05, 0.79149, 0.27201, 0.28295, 7.59e-08,
#'         0.69274, 0.30443, 0.00136, 0.72342, 0.54757)
#'
#' StoreyBH(pvals)
#' 
#' @export

StoreyBH <- function(d, alpha = 0.05, lambda = 0.5){
  
  d <- checkPval(d)
  
  if (is.vector(d)) {
    d <- as.data.frame(d)
    colnames(d) <- "pval"
  } else {
    stop("d must either be a dataframe or a vector of p-values.")
  }
  
  if (alpha < 0 || alpha > 1) {
    stop("alpha must be between 0 and 1.")
  }
  
  if (lambda <= 0 || lambda > 1) {
    stop("lambda must be between 0 and 1.")
  }
  
  pvals <- .subset2(d, "pval")
  o <- order(pvals, decreasing = TRUE)
  ro <- order(o)
  
  #calculate pi0
  n <- length(pvals)
  candsum <- sum(pvals > lambda)
  pi0 <- (candsum + 1)/((1 - lambda) * n)
  
  jvec <- n:1L
  out_R <- pmin(1, cummin(n/jvec * pi0 * pvals[o]))[ro] <= alpha
  
  out <- d
  out$R <- as.numeric(out_R)
  out
}
