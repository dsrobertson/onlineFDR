#' BatchPRDS: Online batch FDR control under Positive Dependence
#'
#' Implements the BatchPRDS algorithm for online FDR control, where PRDS stands
#' for positive regression dependency on a subset, as presented by Zrnic et. al.
#' (2020).
#'
#' The function takes as its input a dataframe with three columns: identifiers
#' (`id'), batch numbers (`batch') and p-values (`pval').
#'
#' The BatchPRDS algorithm controls the FDR when the p-values in one batch are
#' positively dependent, and independent across batches. Given an overall
#' significance level \eqn{\alpha}, we choose a sequence of non-negative numbers
#' \eqn{\gamma_i} such that they sum to 1. The algorithm runs the
#' Benjamini-Hochberg procedure on each batch, where the values of the adjusted
#' significance thresholds \eqn{\alpha_{t+1}} depend on the number of previous 
#' discoveries.
#'
#' Further details of the BatchPRDS algorithm can be found in Zrnic et. al.
#' (2020).
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
#'   Artificial Intelligence and Statistics}: 3806-3815
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
#' BatchPRDS(sample.df)
#'
#' @export

BatchPRDS <- function(d, alpha = 0.05, gammai){
  
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
  n_batch <- length(unique(d$batch))
  if(max(d$batch, na.rm = TRUE) > n_batch) {
    d <- d[order(d$batch),]
    d$batch2 <- rep(1:length(table(d$batch)), table(d$batch))
    d$batch2 <- as.numeric(d$batch2)
    warning("Your batches were not labelled in a fully sequential manner. A new batch id was created labelled sequentially starting from 1.")
  } else {
    d$batch2 <- rep(1:length(table(d$batch)), table(d$batch))
    d$batch2 <- as.numeric(d$batch2)
  }
  
  if (missing(gammai)) {
    gammai <- 0.4374901658/(seq_len(n_batch)^(1.6))
  } else if (any(gammai < 0)) {
    stop("All elements of gammai must be non-negative.")
  } else if (sum(gammai) > 1) {
    stop("The sum of the elements of gammai must not be greater than 1.")
  }
  
  ### Start Batch PRDS procedure
  
  R <- rep(0, nrow(d))
  alphai <- rep(0, n_batch)
  alphai[1] <- gammai[1] * alpha
  
  for(i in seq_len(n_batch)){
    batch_pval <- d[which(d$batch2 == i),]$pval
    # batch_pval <- subset_rcpp(d$batch2, d$pval, i)
    n <- length(batch_pval)
    
    ordered_pval <- sort(batch_pval)
    batchR <- ordered_pval <= ((1:n)/n)*alphai[i]

    max_entry <- suppressWarnings(which.max(batchR))
    if(is.finite(max_entry)) {
      batchR[1:max_entry] <- 1
    }
    out_R <- batchR[order(batch_pval)]
    idx <- which(d$batch2 == i)
    R[idx] <- out_R
    
    #update alphai
    if(i < n_batch) {
      ntplus <- length(which(d$batch2 == i+1))
      alphai[i+1] <- alpha * (gammai[i+1]/ntplus) * (ntplus + sum(R, na.rm = T))
    }
  }
  out <- d
  out$R = R
  out$alphai = rep(alphai, table(d$batch2))
  #deduplicate columns, so if batch2 is the same, remove it
  if(identical(out$batch, out$batch2)) {
    out$batch2 <- NULL
  }
  return(out)
}

