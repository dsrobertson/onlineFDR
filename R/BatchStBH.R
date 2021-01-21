#' BatchStBH: Online batch FDR control using the St-BH procedure
#'
#' Implements the BatchSt-BH algorithm for online FDR control, as presented by
#' Zrnic et al. (2020). This algorithm makes one modification to the original
#' Storey-BH algorithm (Storey 2002), by adding 1 to the numerator of
#' the null proportion estimate for more stable results.
#'
#' The function takes as its input a dataframe with three columns: identifiers
#' (`id'), batch numbers (`batch') and p-values (`pval').
#'
#' The BatchSt-BH algorithm controls the FDR when the p-values in a batch are
#' independent, and independent across batches. Given an overall significance
#' level \eqn{\alpha}, we choose a sequence of non-negative numbers
#' \eqn{\gamma_i} such that they sum to 1. The algorithm runs the
#' Benjamini-Hochberg procedure on each batch, where the values of the adjusted
#' significance thresholds \eqn{\alpha_{t+1}} depend on the number of previous
#' discoveries.
#'
#' Further details of the BatchSt-BH algorithm can be found in Zrnic et al.
#' (2020).
#'
#' @param d A dataframe with three columns: identifiers (`id'), batch numbers
#'   (`batch') and p-values (`pval').
#'
#' @param alpha Overall significance level of the FDR procedure, the default is
#'   0.05.
#'
#' @param gammai Optional vector of \eqn{\gamma_i}. A default is provided with
#'   \eqn{\gamma_j} proportional to \eqn{1/j^(1.6)}.
#'
#' @param lambda Threshold for Storey-BH, must be between 0 and 1. Defaults to
#'   0.5.
#'
#' @return \item{out}{ A dataframe with the original data \code{d} and the
#'   indicator function of discoveries \code{R}. Hypothesis \eqn{i} is rejected
#'   if the \eqn{i}-th p-value within the \eqn{t}-th batch is less than or equal
#'   to \eqn{(r/n)\alpha_t}, where \eqn{r} is the rank of the \eqn{i}-th p-value
#'   within an ordered set and \eqn{n} is the total number of hypotheses within
#'   the \eqn{t}-th batch. If hypothesis \eqn{i} is rejected, \code{R[i] = 1}
#'   (otherwise \code{R[i] = 0}).}
#'
#' @references Storey, J.D. (2002). A direct approach to false discovery rates.
#'   \emph{J. R. Statist. Soc. B}: 64, Part 3, 479-498.
#'
#'   Zrnic, T., Jiang D., Ramdas A. and Jordan M. (2020). The Power of Batching
#'   in Multiple Hypothesis Testing. \emph{International Conference on
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
#' BatchStBH(sample.df)
#'
#' @export

BatchStBH <- function(d, alpha = 0.05, gammai, lambda = 0.5){
  
  if (!is.data.frame(d)) {
    stop("d must be a dataframe")
  } else if (!("batch" %in% colnames(d))) {
    stop("d needs to have a column called batch")
  }
  
  if (alpha < 0 || alpha > 1) {
    stop("alpha must be between 0 and 1.")
  }
  
  if (lambda <= 0 || lambda > 1) {
    stop("lambda must be between 0 and 1.")
  }
  
  if(!is.numeric(d$batch)) {
    stop("Check that your batches are numeric")
  }
  
  #check that batches were labeled correctly
  n_batch <- length(unique(d$batch))
  
  if(max(d$batch, na.rm = TRUE) > n_batch) {
    stop("Check that your batches labelled in ascending order starting from 1")
  }
  
  if (missing(gammai)) {
    gammai <- 0.4374901658/(seq_len(n_batch)^(1.6))
  } else if (any(gammai < 0)) {
    stop("All elements of gammai must be non-negative.")
  } else if (sum(gammai) > 1) {
    stop("The sum of the elements of gammai must not be greater than 1.")
  }
  
  ### Start Batch St-BH procedure
  
  all_batches <- list()
  
  Rplus <- Rsum <- Rrsum <- alphai <- k <- rep(0, n_batch)
  alphai[1] <- gammai[1] * alpha
  
  for(i in seq_len(n_batch)){
    batch_data <- d[d$batch == i,]
    batch_data$ind <- seq.int(nrow(batch_data))
    n <- length(batch_data$pval)
    ordered_batch_data <- batch_data[order(batch_data$pval),]
    candsum <- sum(ordered_batch_data$pval > lambda, na.rm = T)
    
    #calculate pi0
    pi0 <- (candsum + 1)/((1 - lambda)*n)
    
    ordered_batch_data$R <- pi0*ordered_batch_data$pval <= ((1:n)/n)*alphai[i]
    max_entry <- suppressWarnings(max(which(ordered_batch_data$R)))
    if(is.finite(max_entry)) {
      ordered_batch_data$R[1:max_entry] <- 1
    }
    ordered_batch_data <- ordered_batch_data[order(ordered_batch_data$ind),]
    ordered_batch_data$ind <- NULL
    all_batches[[i]] <- ordered_batch_data
    out <- do.call(rbind, all_batches)
    
    Rsum[i] <- sum(ordered_batch_data$R)
    
    ## k
    if(max(ordered_batch_data$pval) > lambda) {
      k[i] <- 1
    }
    
    #calculate Rsplus
    #aug_rej is the number of rejections if we hallucinate jth pval to be 0
    aug_rej <- rep(0,n)
    
    hallucinated_data <- ordered_batch_data
    
    for (j in seq_len(n)) {
      ## Rsplus
      hallucinated_pvals <- c(0, hallucinated_data$pval[-j])
      #calulcate pi0
      hallucinated_pi0 <- (sum(hallucinated_pvals > lambda) + 1)/((1 - lambda)*n)
      #run St-BH procedure with hallucinated p-value
      hallucinated_data$R <- hallucinated_pi0*hallucinated_pvals <= ((1:n)/n)*alphai[i]
      max_entry <- suppressWarnings(max(which(hallucinated_data$R)))
      if(is.finite(max_entry)) {
        hallucinated_data$R[1:max_entry] <- 1
      }
      
      aug_rej[j] <- sum(hallucinated_data$R, na.rm = T)
      
    }
    
    Rplus[i] = max(aug_rej)
    
    #update alphai
    if(i < n_batch) {
      gammasum <- sum(gammai[seq_len(i+1)]) * alpha
      
      for (r in seq_len(i)) {
        Rrsum[r] = sum(Rsum[-r])
      }
      
      ntplus <- nrow(d[d$batch == i+1,])
      alphai[i+1] <- (gammasum - sum(k[1:i]*alphai[1:i]*(Rplus[1:i]/(Rplus[1:i] + Rrsum[1:i])))) * 
        ((ntplus + sum(Rsum))/ntplus)
      
    }
  }
  out$alphai <- rep(alphai, table(d$batch))
  return(out)
}
