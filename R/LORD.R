#' Online FDR control based on recent discovery
#'
#' Implements the LORD procedure for online FDR control where LORD stands for
#' (significance) Levels based On Recent Discovery, as presented by
#' Javanmard and Montanari (2018).
#'
#' The function takes as its input a dataframe with three columns: an identifier
#' (`id'), date (`date') and p-value (`pval'). The case where p-values arrive in
#' batches corresponds to multiple instances of the same date. If no column of
#' dates is provided, then the p-values are treated as being ordered
#' sequentially with no batches.
#'
#' The LORD procedure controls FDR for independent p-values. Given an overall
#' significance level \eqn{\alpha}, we choose a sequence of
#' non-negative numbers \eqn{\gamma_i} such that they sum to 1, and
#' \eqn{\gamma_i \geq \gamma_j} for \eqn{i \leq j}.
#'
#' Javanmard and Montanari (2018) present three versions of LORD which
#' differ in the way the adjusted test levels \eqn{\alpha_i} are calculated. The
#' test levels for LORD 1 are based on the time of the last discovery
#' (i.e. hypothesis rejection), LORD 2 are based on all previous discovery
#' times, and LORD 3 are based on the time of the last discovery as well as
#' the 'wealth' accumulated at that time.
#'
#' LORD depends on constants \eqn{w_0} and \eqn{b_0}, where \eqn{w_0 \ge 0}
#' represents the intial `wealth' of the procedure and \eqn{b_0 > 0} represents
#' the `payout' for rejecting a hypothesis. We require \eqn{w_0+b_0 \le \alpha}
#' for FDR control to hold.
#'
#' Note that FDR control also holds for the LORD procedure if only the p-values
#' corresponding to true nulls are mutually independent, and independent from
#' the non-null p-values.
#'
#' Further details of the LORD procedure can be found in Javanmard and
#' Montanari (2018).
#'
#'
#' @param d Dataframe with three columns: an identifier (`id'), date (`date')
#' and p-value (`pval'). If no column of dates is provided, then the p-values
#' are treated as being ordered sequentially with no batches.
#'
#' @param alpha Overall significance level of the FDR procedure, the default
#' is 0.05.
#'
#' @param gammai Optional vector of \eqn{\gamma_i}. A default is provided as
#' proposed by Javanmard and Montanari (2018), equation 31.
#'
#' @param version An integer from 1 to 3 giving the version of LORD to use.
#' Defaults to 3.
#'
#' @param w0 Initial `wealth' of the procedure. Defaults to \eqn{\alpha/10}.
#' @param b0 The 'payout' for rejecting a hypothesis. Defaults to
#' \eqn{\alpha - w_0}.
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
#' Javanmard, A. and Montanari, A. (2018) Online Rules for Control of False
#' Discovery Rate and False Discovery Exceedance. \emph{Annals of Statistics},
#' 46(2):526-554.
#'
#'
#' @seealso
#'
#' \code{\link{LORDdep}} uses a modified version of the LORD algorithm that is
#' valid for \emph{dependent} p-values.
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
#'         3.60e-05, 0.79149, 0.27201, 0.28295, 7.59e-08,
#'         0.69274, 0.30443, 0.00136, 0.72342, 0.54757))
#'
#' LORD(sample.df, random=FALSE)
#' set.seed(1); LORD(sample.df, version=2)
#' set.seed(1); LORD(sample.df, alpha=0.1, w0=0.05)
#'
#'
#' @export

LORD <- function(d, alpha=0.05, gammai, version=3, w0=alpha/10, b0=alpha-w0,
                random=TRUE, date.format="%Y-%m-%d") {

    if(!(is.data.frame(d))){
        stop("d must be a dataframe.")
    }

    if(length(d$id) == 0){
        stop("The dataframe d is missing a column 'id' of identifiers.")
    } else if(length(d$pval) == 0){
        stop("The dataframe d is missing a column 'pval' of p-values.")
    }

    if(length(d$date) == 0){
        warning("No column of dates is provided, so p-values are treated
        as being ordered sequentially with no batches.")
        random = FALSE
    } else if(any(is.na(as.Date(d$date, date.format)))){
        stop("One or more dates are not in the correct format.")
    } else {
        d <- d[order(as.Date(d$date, format = date.format)),]
    }

    if(alpha<=0 || alpha>1){
        stop("alpha must be between 0 and 1.")
    }

    if(!(version %in% c(1,2,3))){
        stop("version must be 1, 2 or 3.")
    }

    if(anyNA(d$pval)){
        warning("Missing p-values were ignored.")
        d <- stats::na.omit(d)
    }

    if(!(is.numeric(d$pval))){
        stop("The column of p-values contains at least one non-numeric
        element.")
    } else if(any(d$pval>1 | d$pval<0)){
        stop("All p-values must be between 0 and 1.")
    }

    N <- length(d$pval)

    if(missing(gammai)){
        gammai <- 0.07720838*log(pmax(seq_len(N),2)) /
                    (seq_len(N)*exp(sqrt(log(seq_len(N)))))
    } else if (any(gammai<0)){
        stop("All elements of gammai must be non-negative.")
    } else if(sum(gammai)>1){
        stop("The sum of the elements of gammai must not be greater than 1.")
    }

    if(w0 < 0){
        stop("w0 must be non-negative.")
    } else if(b0 <= 0){
        stop("b0 must be positive.")
    } else if(w0+b0 > alpha & !(isTRUE(all.equal(w0+b0, alpha)))){
        stop("The sum of w0 and b0 must not be greater than alpha.")
    }

    if(random){
        d <- randBatch(d)
    }

    alphai <- rep(0, N)
    pval <- d$pval

    switch(version,
        ## 1
        {
        R <- rep(0, N)

        for (i in seq_len(N)){
            tau <- max(0, which(R[seq_len(i-1)] == 1))
            alphai[i] <- gammai[i]*w0*(tau == 0) + gammai[i-tau]*b0*(tau > 0)
            R[i] <- pval[i] <= alphai[i]
            }
        },
        ## 2
        {
        R <- rep(0, N)
        alphai[1] <- gammai[1]*w0
        R[1] <- pval[1] <= alphai[1]

        if(N == 1){
            d.out <- data.frame(d, alphai, R)
            return(d.out)
        }

        for (i in (seq_len(N-1)+1)){
            alphai[i] <- gammai[i]*w0 +
                            sum(gammai[i-which(R[seq_len(i-1)] == 1)])*b0

            R[i] <- pval[i] <= alphai[i]
            }
        },
        ## 3
        {
        R <- W <- rep(0, N+1)
        R[1] <- 1
        W[1] <- w0

        alphai[1] <- phi <- gammai[1]*w0
        R[2] <- pval[1] <= alphai[1]
        W[2] <- w0 - phi + R[2]*b0

        if(N == 1){
            R <- R[2]
            d.out <- data.frame(d, alphai, R)
            return(d.out)
        }

        for (i in (seq_len(N-1)+1)){
            tau <- max(which(R[seq_len(i)] == 1))
            alphai[i] <- phi <- gammai[i-tau+1]*W[tau]

            R[i+1] <- pval[i] <= alphai[i]
            W[i+1] <- W[i] - phi + R[i+1]*b0
            }
        R <- R[(seq_len(N)+1)]
        })

    d.out <- data.frame(d, alphai, R)

    return(d.out)
}
