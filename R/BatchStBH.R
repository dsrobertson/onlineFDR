#' BatchStBH: Online batch FDR control using the St-BH procedure
#'
#' Implements the BatchSt-BH algorithm for online FDR control, as presented by
#' Zrnic et al. (2020). Within each batch it uses the Storey-BH procedure
#' (Storey 2002) with the null proportion estimate of Storey, Taylor and
#' Siegmund (2004), which adds 1 to the numerator for more stable results.
#'
#' The function takes as its input a dataframe with three columns: identifiers
#' (`id'), batch numbers (`batch') and p-values (`pval').
#'
#' The BatchSt-BH algorithm controls the FDR when the p-values in a batch are
#' independent, and independent across batches. Given an overall significance
#' level \eqn{\alpha}, we choose a sequence of non-negative numbers
#' \eqn{\gamma_i} such that they sum to 1. The algorithm runs the
#' Storey Benjamini-Hochberg procedure on each batch, where the values of the adjusted
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
#'   \eqn{\gamma_j} proportional to \eqn{1/j^{1.6}}.
#'
#' @param lambda Threshold for Storey-BH, must be between 0 and 1. Defaults to
#'   0.5.
#'   
#' @param display_progress Logical. If \code{TRUE} prints out a progress bar for the algorithm runtime. 
#'
#' @return \item{out}{ A dataframe with the original data \code{d}, the
#'   indicator function of discoveries \code{R} and the test level
#'   \code{alphai} (\eqn{\alpha_t}) of each hypothesis's batch. Hypothesis
#'   \eqn{i} in the \eqn{t}-th batch is rejected if its p-value is less than or
#'   equal to \eqn{(k/n)\alpha_t / \hat{\pi}_0}, where \eqn{n} is the number of
#'   hypotheses in the batch, \eqn{k} is the largest \eqn{r} such that the
#'   \eqn{r}-th smallest p-value in the batch is at most \eqn{(r/n)\alpha_t /
#'   \hat{\pi}_0}, and \eqn{\hat{\pi}_0 = (1 + n_\lambda)/(n(1 - \lambda))},
#'   with \eqn{n_\lambda} the number of p-values in the batch greater than
#'   \eqn{\lambda}. If hypothesis \eqn{i} is rejected, \code{R[i] = 1}
#'   (otherwise \code{R[i] = 0}).}
#'
#' @references Storey, J.D. (2002). A direct approach to false discovery rates.
#'   \emph{Journal of the Royal Statistical Society (Series B)}, 64(3):479-498.
#'
#'   Storey, J.D., Taylor, J.E. and Siegmund, D. (2004). Strong control,
#'   conservative point estimation and simultaneous conservative consistency of
#'   false discovery rates: a unified approach. \emph{Journal of the Royal
#'   Statistical Society (Series B)}, 66(1):187-205.
#'
#'   Zrnic, T., Jiang, D., Ramdas, A. and Jordan, M.I. (2020). The Power of
#'   Batching in Multiple Hypothesis Testing. \emph{International Conference on
#'   Artificial Intelligence and Statistics}, PMLR, 108:3806-3815.
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

BatchStBH <- function(d, alpha = 0.05, gammai, lambda = 0.5, display_progress = FALSE){
  
  d <- checkPval(d)
  
  if (!is.data.frame(d)) {
    stop("d must be a dataframe")
  } else if (!("batch" %in% colnames(d))) {
    stop("d needs to have a column called batch")
  }
  
  if (alpha < 0 || alpha > 1) {
    stop("alpha must be between 0 and 1.")
  }
  
  if (lambda <= 0 || lambda >= 1) {
    stop("lambda must be between 0 and 1.")
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
  } else if (sumExceeds(gammai, 1)) {
    stop("The sum of the elements of gammai must not be greater than 1.")
  } else if (length(gammai) < n_batch) {
    stop("gammai must have length at least equal to the number of batches.")
  }
  
  ### Start Batch St-BH procedure
  
  R <- NULL
  Rplus <- Rsum <- Rrsum <- alphai <- k <- rep(0, n_batch)
  alphai[1] <- gammai[1] * alpha
  
  nt <- as.vector(table(d$batch))
  batch_indices <- c(0, cumsum(nt))
  
  if(display_progress) {
    pb <- progress::progress_bar$new(format = "  Computing [:bar] :percent eta: :eta",
                                     total = n_batch, clear = FALSE, width = 60)
    
    for(i in seq_len(n_batch)) {
      pb$tick()
      idx_b <- batch_indices[i]+1
      idx_e <- batch_indices[i+1]
      batch_pval <- .subset2(d, "pval")[idx_b:idx_e]
      
      jvec <- nt[i]:1L
      #sort pvals and then return the original indices of the sorted pvals
      o <- order(batch_pval, decreasing = TRUE)
      #sort the indices and then return the indices of the sorted indices
      #effectively reverses the order
      ro <- order(o)
      
      #calculate pi0
      n <- length(batch_pval)
      candsum <- sum(batch_pval > lambda)
      pi0 <- (candsum + 1)/((1 - lambda) * n)
      
      out_R <- cummin(nt[i]/jvec * pi0 * batch_pval[o])[ro] <= alphai[i]
      
      R <- c(R, out_R)
      
      Rsum[i] <- sum(out_R)
      
      ## k
      if(max(batch_pval) > lambda) {
        k[i] <- 1
      }
      
      #calculate Rsplus
      #aug_rej is the number of rejections if we hallucinate jth pval to be 0
      aug_rej <- rep(0,nt[i])
      
      for (j in seq_len(nt[i])) {
        
        #run St-BH procedure with hallucinated p-value
        hallucinated_pval <- batch_pval
        hallucinated_pval[j] <- 0
        oh <- order(hallucinated_pval, decreasing = TRUE)
        roh <- order(oh)
        
        #calculate pi0
        hallucinated_pi0 <- (sum(hallucinated_pval > lambda) + 1)/((1 - lambda)*n)
        hallucinated_R <- cummin(nt[i]/jvec * hallucinated_pi0 *
                                    hallucinated_pval[oh])[roh] <= alphai[i]
        
        aug_rej[j] <- sum(hallucinated_R, na.rm = TRUE)
      }
      
      Rplus[i] = max(aug_rej)
      
      #update alphai
      if(i < n_batch) {
        gammasum <- sum(gammai[seq_len(i+1)]) * alpha
        
        Rrsum[1:i] = sum(Rsum)-Rsum[1:i]
        
        alphai[i+1] <- (gammasum - sum(k[1:i]*alphai[1:i]*(Rplus[1:i]/(Rplus[1:i] + Rrsum[1:i])))) * 
          ((nt[i+1] + sum(Rsum))/nt[i+1])
      }
    }
    out <- d
    out$R <- as.numeric(R)
    out$alphai <- rep(alphai, nt)
    out
  } else {
    for(i in seq_len(n_batch)) {
      idx_b <- batch_indices[i]+1
      idx_e <- batch_indices[i+1]
      batch_pval <- .subset2(d, "pval")[idx_b:idx_e]
      
      jvec <- nt[i]:1L
      #sort pvals and then return the original indices of the sorted pvals
      o <- order(batch_pval, decreasing = TRUE)
      #sort the indices and then return the indices of the sorted indices
      #effectively reverses the order
      ro <- order(o)
      
      #calculate pi0
      n <- length(batch_pval)
      candsum <- sum(batch_pval > lambda)
      pi0 <- (candsum + 1)/((1 - lambda) * n)
      
      out_R <- cummin(nt[i]/jvec * pi0 * batch_pval[o])[ro] <= alphai[i]
      
      R <- c(R, out_R)
      
      Rsum[i] <- sum(out_R)
      
      ## k
      if(max(batch_pval) > lambda) {
        k[i] <- 1
      }
      
      #calculate Rsplus
      #aug_rej is the number of rejections if we hallucinate jth pval to be 0
      aug_rej <- rep(0,nt[i])
      
      for (j in seq_len(nt[i])) {
        
        #run St-BH procedure with hallucinated p-value
        hallucinated_pval <- batch_pval
        hallucinated_pval[j] <- 0
        oh <- order(hallucinated_pval, decreasing = TRUE)
        roh <- order(oh)
        
        #calculate pi0
        hallucinated_pi0 <- (sum(hallucinated_pval > lambda) + 1)/((1 - lambda)*n)
        hallucinated_R <- cummin(nt[i]/jvec * hallucinated_pi0 *
                                    hallucinated_pval[oh])[roh] <= alphai[i]
        
        aug_rej[j] <- sum(hallucinated_R, na.rm = TRUE)
      }
      
      Rplus[i] = max(aug_rej)
      
      #update alphai
      if(i < n_batch) {
        gammasum <- sum(gammai[seq_len(i+1)]) * alpha
        
        Rrsum[1:i] = sum(Rsum)-Rsum[1:i]
        
        alphai[i+1] <- (gammasum - sum(k[1:i]*alphai[1:i]*(Rplus[1:i]/(Rplus[1:i] + Rrsum[1:i])))) * 
          ((nt[i+1] + sum(Rsum))/nt[i+1])
      }
    }
    out <- d
    out$R <- as.numeric(R)
    out$alphai <- rep(alphai, nt)
    out
  }
}
  