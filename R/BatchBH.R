#' BatchBH: Online batch FDR control using the BH procedure
#'
#' Implements the BatchBH algorithm for online FDR control, as presented by 
#' Zrnic et al. (2020).
#'
#' The function takes as its input a dataframe with three columns: identifiers
#' (`id'), batch numbers (`batch') and p-values (`pval').
#'
#' The BatchBH algorithm controls the FDR when the p-values in a batch 
#' are independent, and independent across batches. Given an overall
#' significance level \eqn{\alpha}, we choose a sequence of non-negative numbers
#' \eqn{\gamma_i} such that they sum to 1. The algorithm runs the
#' Benjamini-Hochberg procedure on each batch, where the values of the adjusted
#' significance thresholds \eqn{\alpha_{t+1}} depend on the number of previous 
#' discoveries.
#'
#' Further details of the BatchBH algorithm can be found in Zrnic et al. (2020).
#'
#' @param d A dataframe with three columns: identifiers (`id'),
#'   batch numbers (`batch') and p-values (`pval').
#'
#' @param alpha Overall significance level of the FDR procedure, the default is
#'   0.05.
#'
#' @param gammai Optional vector of \eqn{\gamma_i}. A default is provided with
#'   \eqn{\gamma_j} proportional to \eqn{1/j^(1.6)}.
#'
#' @return \item{out}{ A dataframe with the original data \code{d} and the
#'   indicator function of discoveries \code{R}. Hypothesis \eqn{i} is rejected
#'   if the \eqn{i}-th p-value within the \eqn{t}-th batch is less than or equal
#'   to \eqn{(r/n)\alpha_t}, where \eqn{r} is the rank of the \eqn{i}-th p-value
#'   within an ordered set and \eqn{n} is the total number of hypotheses within
#'   the \eqn{t}-th batch. If hypothesis \eqn{i} is rejected, \code{R[i] = 1}
#'   (otherwise \code{R[i] = 0}).}
#'
#' @references Zrnic, T., Jiang D., Ramdas A. and Jordan M. (2020). The Power of
#'   Batching in Multiple Hypothesis Testing. \emph{International Conference on
#'   Artificial Intelligence and Statistics}, 3806-3815.
#'
#' @examples
#'
#' sample.df <- data.frame(
#' id = c('A15432', 'B90969', 'C18705', 'B49731', 'E99902',
#'     'C38292', 'A30619', 'D46627', 'E29198', 'A41418',
#'     'D51456', 'C88669', 'E03673', 'A63155', 'B66033'),
#' pval = c(2.90e-08, 0.06743, 0.01514, 0.08174, 0.00171,
#'         3.60e-05, 0.79149, 0.27201, 0.28295, 7.59e-08,
#'         0.69274, 0.30443, 0.00136, 0.72342, 0.54757),
#' batch = c(rep(1,5), rep(2,6), rep(3,4)))
#'
#' BatchBH(sample.df)
#'
#' @export

BatchBH <- function(d, alpha = 0.05, gammai){
  
  d <- checkPval(d)
  
  if (!is.data.frame(d)) {
    stop("d must be a dataframe")
  } else if (!("batch" %in% colnames(d))) {
    stop("d needs to have a column called batch")
  }
  
  if (alpha < 0 || alpha > 1) {
    stop("alpha must be between 0 and 1.")
  }
  
  if(!is.numeric(d$batch)) {
    stop("Check that your batch labels are numeric values.")
  }
  
  #check that batches were labeled correctly
  
  if(is.unsorted(d$batch)) {
    d <- d[order(d$batch),]
    warning("Batches were re-ordered in increasing numeric value.")
  }
  
  n_batch <- length(unique(d$batch))
  if (missing(gammai)) {
    gammai <- 0.4374901658/(seq_len(n_batch)^(1.6))
  } else if (any(gammai < 0)) {
    stop("All elements of gammai must be non-negative.")
  } else if (sum(gammai) > 1) {
    stop("The sum of the elements of gammai must not be greater than 1.")
  }
  
  ### Start Batch BH procedure
  R <- NULL
  Rplus <- Rsum <- Rrsum <- alphai <- rep(0, n_batch)
  alphai[1] <- gammai[1] * alpha
  
  nt <- as.vector(table(d$batch))
  batch_indices <- c(0, cumsum(nt))
  
  for(i in seq_len(n_batch)){
    idx_b <- batch_indices[i]+1
    idx_e <- batch_indices[i+1]
    batch_pval <- .subset2(d, "pval")[idx_b:idx_e]
    
    k <- nt[i]:1L
    #sort pvals and then return the original indices of the sorted pvals
    o <- order(batch_pval, decreasing = TRUE)
    #sort the indices and then return the indices of the sorted indices
    #effectively reverses the order
    ro <- order(o)
    out_R <- pmin(1, cummin(nt[i]/k * batch_pval[o]))[ro] <= alphai[i]
    
    R <- c(R, out_R)
    
    Rsum[i] <- sum(out_R)
    
    #calculate Rsplus
    aug_rej <- rep(0,nt[i])
    
    for (j in seq_len(nt[i])) {
      
      #run BH procedure with hallucinated p-value
      hallucinated_pval <- batch_pval
      hallucinated_pval[j] <- 0
      oh <- order(hallucinated_pval, decreasing = TRUE)
      roh <- order(oh)
      hallucinated_R <- pmin(1, cummin(nt[i]/k * hallucinated_pval[oh]))[roh] <= alphai[i]
      aug_rej[j] <- sum(hallucinated_R)
    }
    Rplus[i] = max(aug_rej)
    
    #update alphai
    if(i < n_batch) {
      gammasum <- sum(gammai[seq_len(i+1)]) * alpha
      
      Rrsum[1:i] = sum(Rsum)-Rsum[1:i]
      
      alphai[i+1] <- (gammasum - sum(alphai[1:i]*(Rplus[1:i]/(Rplus[1:i] + Rrsum[1:i])))) * 
        ((nt[i+1] + sum(Rsum))/nt[i+1])
    }
  }
  out <- d
  out$R <- as.numeric(R)
  out$alphai <- rep(alphai, nt)
  out
}
