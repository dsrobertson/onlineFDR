#' Online FDR control based on number of discoveries
#'
#' Implements the LOND algorithm for online FDR control, where LOND stands for
#' (significance) Levels based On Number of Discoveries, as presented by
#' Javanmard and Montanari (2015).
#'
#' The function takes as its input a dataframe with three columns: an identifier
#' (`id'), date (`date') and p-value (`pval'). The case where p-values arrive in
#' batches corresponds to multiple instances of the same date. If no column of
#' dates is provided, then the p-values are treated as being ordered
#' sequentially with no batches.
#'
#' The LOND algorithm controls FDR for independent p-values. Given an overall
#' significance level \eqn{\alpha}, we choose a sequence of
#' non-negative numbers \eqn{\beta_i} such that they sum to \eqn{\alpha}. The
#' values of the adjusted test levels \eqn{\alpha_i} are chosen as follows:
#' \deqn{\alpha_i = (D(i-1) + 1)\beta_i}
#' where \eqn{D(n)} denotes the number of discoveries in the first \eqn{n}
#' hypotheses.
#'
#' For dependent p-values, LOND controls FDR if it is modified with
#' \eqn{\beta_i / H(i)} in place of \eqn{\beta_i}, where \eqn{H(j)} is the i-th
#' harmonic number.
#'
#' Further details of the LOND algorithm can be found in Javanmard and
#' Montanari (2015).
#'
#'
#' @param d Dataframe with three columns: an identifier (`id'), date (`date')
#' and p-value (`pval'). If no column of dates is provided, then the p-values
#' are treated as being ordered sequentially with no batches.
#'
#' @param alpha Overall significance level of the FDR procedure, the default
#' is 0.05.
#'
#' @param beta Optional vector of \eqn{\beta_i}. A default is provided as
#' proposed by Javanmard and Montanari (2017), equation 31.
#'
#' @param dep Logical. If \code{TRUE}, runs the modified LOND algorithm which
#' guarantees FDR control for \emph{dependent} p-values.
#' Defaults to \code{FALSE}.
#'
#' @param random Logical. If \code{TRUE} (the default), then the order of the
#' p-values in each batch (i.e. those that have exactly the same date) is
#' randomised.
#'
#' @param date.format Optional string giving the format that is used for dates.
#'
#'
#' @return
#' \item{d.out}{ A dataframe with the original dataframe \code{d} (which will
#' be reordered if there are batches and \code{random = TRUE}), the
#' LOND-adjusted test levels \eqn{\alpha_i} and the indicator function of
#' discoveries \code{R}. Hypothesis \eqn{i} is rejected if the \eqn{i}-th
#' p-value is less than or equal to \eqn{\alpha_i}, in which case
#' \code{R[i] = 1}  (otherwise \code{R[i] = 0}).}
#'
#'
#' @references
#' Javanmard, A. and Montanari, A. (2015) On Online Control of False Discovery
#' Rate. \emph{arXiv preprint}, \url{https://arxiv.org/abs/1502.06197}
#'
#' Javanmard, A. and Montanari, A. (2017) Online Rules for Control of False
#' Discovery Rate and False Discovery Exceedance. \emph{Accepted for publication
#' in Annals of Statistics}, available at
#' \url{https://arxiv.org/abs/1603.09000}.
#'
#'
#' @examples
#' sample.df <- data.frame(
#' id = c('A15432', 'B90969', 'C18705', 'B49731', 'E99902',
#'     'C38292', 'A30619', 'D46627', 'E29198', 'A41418',
#'     'D51456', 'C88669', 'E03673', 'A63155', 'B66033'),
#' date = as.Date(c(rep("2014-12-01",3),
#'                 rep("2015-09-21",5),
#'                 rep("2016-05-19",2),
#'                 "2016-11-12",
#'                 rep("2017-03-27",4))),
#' pval = c(2.90e-17, 0.06743, 0.01514, 0.08174, 0.00171,
#'         3.60E-05, 0.79149, 0.27201, 0.28295, 7.59E-08,
#'         0.69274, 0.30443, 0.00136, 0.72342, 0.54757))
#'
#' LOND(sample.df)
#' LOND(sample.df, random=FALSE)
#' LOND(sample.df, alpha=0.1)
#'

LOND <- function(d, alpha=0.05, beta, dep=FALSE, random=TRUE,
                date.format="%Y-%m-%d") {

    if(length(d$id) == 0){
        stop("The dataframe d is missing a column 'id' of identifiers.")
    } else if(length(d$pval) == 0){
        stop("The dataframe d is missing a column 'pval' of p-values.")
    }

    if(length(d$date) == 0){
        warning("No column of dates is provided, so p-values are treated as
        being ordered sequentially with no batches.")
        random = FALSE
    } else if(any(is.na(as.Date(d$date, date.format)))){
        stop("One or more dates are not in the correct format.")
    } else {
        d <- d[order(as.Date(d$date, format = date.format)),]
    }

    if(alpha<0 || alpha>1 ){
        stop("alpha must be between 0 and 1.")
    }

    if(anyNA(d$pval)){
        warning("Missing p-values were ignored.")
        d <- na.omit(d)
    }

    if(!(is.numeric(d$pval))){
        stop("The column of p-values contains at least one non-numeric
        element.")
    } else if(any(d$pval>1 | d$pval<0)){
        stop("All p-values must be between 0 and 1.")
    }

    N <- length(d$pval)

    if(missing(beta)){
        beta <- 0.07720838*alpha*log(pmax(1:N,2))/((1:N)*exp(sqrt(log(1:N))))
    } else if (any(beta<0)){
        stop("All elements of beta must be non-negative.")
    } else if(sum(beta)>alpha){
        stop("The sum of the elements of beta must not be greater than alpha.")
    }

    if(dep) {
        for(i in 1:N){
            beta[i] <- beta[i]/sum(1/(1:i))
        }
    }

    if(random){
        Nbatch <- length(unique(d$date))
        set.seed(1)

        for(i in 1:Nbatch){
            d.temp <- d[d$date == unique(d$date)[i],]
            d.temp <- d.temp[sample.int(length(d.temp$date)),]
            d[d$date == unique(d$date)[i],] <- d.temp
        }
    }

    pval <- d$pval

    D <- R <- alphai <- rep(0, N)

    alphai[1] <- beta[1]
    R[1] <- pval[1] <= alphai[1]
    D[1] <- R[1]

    for (i in 2:N){
        alphai[i] <- beta[i]*(D[i-1]+1)
        R[i] <- pval[i] <= alphai[i]
        D[i] <- sum(R[1:i])
    }

    d.out <- data.frame(d, alphai, R)

    return(d.out)
}
