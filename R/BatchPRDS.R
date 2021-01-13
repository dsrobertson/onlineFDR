#' BatchPRDS: Online batch FDR control under Positive Dependence
#' 
#' Implements the BatchPRDS algorithm for online FDR control, where PRDS stands for
#' positive regression dependency on a subset, as presented by
#' Zrnic et. al. (2020).
#'
#' The function takes as its input a dataframe
#' with three columns: an identifier (`id'), 
#' and a p-value (`pval').
#'
#' The BatchPRDS algorithm controls the FDR when the p-values in one batch
#' are positively dependent, and independent across batches.
#' Given an overall significance level \eqn{\alpha}, 
#' we choose a sequence of non-negative numbers \eqn{\gammi_i} such
#' that they sum to \eqn{\alpha}. The algorithm runs the Benjamini-Hochberg procedure
#' on each batch, where the values of the adjusted significance
#' thresholds \eqn{\alpha_{t+1}} are chosen as follows: \deqn{\alpha_{t+1} = \alpha \frac{\gammai_{t+1}}{n_{t+1}}(n_{t+1} + \sum_{s=1}^{t}R_s
#' where \eqn{R_s} denotes the number of discoveries in all the batches up until \eqn{batch_t}.
#' 
#' Further details of the LOND algorithm can be found in Zrnic et. al. (2020).
#' 
#' @param d A dataframe with three columns: an identifier (`id') 
#' and a p-value (`pval').
#'
#' @param alpha Overall significance level of the FDR procedure, the default is
#'   0.05.
#'   
#' @param gammai Optional vector of \eqn{\gamma_i}. A default is provided with
#' \eqn{\gamma_j} proportional to \eqn{1/j^(1.6)}.
#' 
#' @return \item{out}{ A dataframe with the original data \code{d} and the indicator
#' function of discoveries \code{R}. Hypothesis \eqn{i} is rejected if the
#' \eqn{i}-th p-value within the \eqn{t}-th batch is less than or equal to 
#' \eqn{(r/n)*\alpha_t}, where \eqn{r} is the rank of the \eqn{i}-th p-value 
#' within an ordered set and \eqn{n} is the total number of hypotheses 
#' within the \eqn{t}-th batch. As such, \code{R[i] = 1}  (otherwise \code{R[i] = 0}).}
#' 
#' @references Zrnic, T., Jiang D., Ramdas A. and Jordan M. (2020). The Power of Batching in Multiple Hypothesis Testing. \emph{International Conference on Artificial Intelligence and Statistics}: 3806-3815 
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
  
  if (missing(gammai)) {
    gammai <- 0.4374901658/(seq_len(length(d$pval))^(1.6))
  } else if (any(gammai < 0)) {
    stop("All elements of gammai must be non-negative.")
  } else if (sum(gammai) > 1) {
    stop("The sum of the elements of gammai must not be greater than 1.")
  }
  
  ### Start Batch PRDS procedure
  alphai <- c()
  alphai[1] <- gammai[1] * alpha
  n_batch <- length(unique(d$batch))
  all_batches <- list()
  
  for(i in seq_along(1:n_batch)){
    batch_data <- d[d$batch == i,]
    batch_data$ind <- seq.int(nrow(batch_data))
    n <- length(batch_data$pval)
    ordered_batch_data <- batch_data[order(batch_data$pval),]
    ordered_batch_data$R <- ordered_batch_data$pval <= ((1:n)/n)*alphai[i]
    max_entry <- suppressWarnings(max(which(ordered_batch_data$R)))
    if(is.finite(max_entry)) {
      ordered_batch_data$R[1:max_entry] <- 1
    }
    ordered_batch_data <- ordered_batch_data[order(ordered_batch_data$ind),]
    ordered_batch_data$ind <- NULL
    all_batches[[i]] <- ordered_batch_data
    out <- do.call(rbind, all_batches)
    
    #update alphai
    if(i < n_batch) {
      ntplus <- nrow(d[d$batch == i+1,])
      alphai[i+1] <- alpha * (gammai[i+1]/ntplus) * (ntplus + sum(out$R, na.rm = T))
    }
  }
  return(out)
}