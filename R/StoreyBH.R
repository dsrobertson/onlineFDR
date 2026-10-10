#' StoreyBH: Offline FDR control using the St-BH procedure
#'
#' Implements the Storey-BH algorithm for offline FDR control, as presented by
#' Storey (2002), with the null proportion estimate of Storey, Taylor and
#' Siegmund (2004).
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
#' @param lambda Threshold for Storey-BH, must be in (0, 1). Defaults to
#'   0.5.
#'
#' @return \item{out}{ A dataframe with the original data \code{d} and the
#'   indicator function of discoveries \code{R}. Hypothesis \eqn{i} is rejected
#'   if the \eqn{i}-th p-value is less than or equal to \eqn{(k/n)\alpha /
#'   \hat{\pi}_0}, where \eqn{k} is the largest \eqn{r} such that the
#'   \eqn{r}-th smallest p-value is at most \eqn{(r/n)\alpha / \hat{\pi}_0},
#'   \eqn{n} is the total number of hypotheses and \eqn{\hat{\pi}_0 = (1 +
#'   n_\lambda) / (n(1 - \lambda))} is Storey's estimate of the proportion of
#'   true null hypotheses, with \eqn{n_\lambda} the number of p-values greater
#'   than \eqn{\lambda}. If hypothesis \eqn{i} is
#'   rejected, \code{R[i] = 1} (otherwise \code{R[i] = 0}).}
#'
#' @references Storey, J.D. (2002). A direct approach to false discovery rates.
#'   \emph{Journal of the Royal Statistical Society (Series B)}, 64(3):479-498.
#'
#'   Storey, J.D., Taylor, J.E. and Siegmund, D. (2004). Strong control,
#'   conservative point estimation and simultaneous conservative consistency of
#'   false discovery rates: a unified approach. \emph{Journal of the Royal
#'   Statistical Society (Series B)}, 66(1):187-205.
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
  } else if (!is.data.frame(d)) {
    stop("d must either be a dataframe or a vector of p-values.")
  }
  
  if (alpha < 0 || alpha > 1) {
    stop("alpha must be between 0 and 1.")
  }
  
  if (lambda <= 0 || lambda >= 1) {
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
  out_R <- cummin(n/jvec * pi0 * pvals[o])[ro] <= alpha
  
  out <- d
  out$R <- as.numeric(out_R)
  out
}
