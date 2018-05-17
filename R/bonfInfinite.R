#' Online FDR control based on a Bonferroni-like test
#'
#' Implements online FDR control using a Bonferroni-like test.
#'
#' The function takes as its input a dataframe with three columns: an identifier
#' (`id'), date (`date') and p-value (`pval'). The case where p-values arrive in
#' batches corresponds to multiple instances of the same date. If no column of
#' dates is provided, then the p-values are treated as being ordered
#' sequentially with no batches.
#'
#' The procedure controls FDR for a potentially infinite stream of p-values by
#' using a Bonferroni-like test. Given an overall significance level
#' \eqn{\alpha}, we choose a (potentially infinite) sequence of non-negative
#' numbers \eqn{\alpha_i} such that they sum to \eqn{\alpha}.
#' Hypothesis \eqn{i} is rejected if the \eqn{i}-th p-value is less than or
#' equal to \eqn{\alpha_i}.
#'
#'
#' @param d Dataframe with three columns: an identifier (`id'), date (`date')
#' and p-value (`pval'). If no column of dates is provided, then the p-values
#' are treated as being ordered sequentially with no batches.
#'
#' @param alpha Overall significance level of the FDR procedure, the default
#' is 0.05.
#'
#' @param alphai Optional vector of \eqn{\alpha_i}, where hypothesis \eqn{i} is
#' rejected if the \eqn{i}-th p-value is less than or equal to \eqn{\alpha_i}.
#' A default is provided as proposed by Javanmard and Montanari (2017),
#' equation 31.
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
#' be reordered if there are batches and \code{random = TRUE}),
#' the test levels \code{alphai} and the indicator function of
#' discoveries \code{R}, where \code{R[i] = 1} corresponds to
#' hypothesis \eqn{i} being rejected (otherwise \code{R[i] = 0}).}
#'
#'
#' @references
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
#' bonfInfinite(sample.df)
#' bonfInfinite(sample.df, random=FALSE)
#' bonfInfinite(sample.df, alpha=0.1)
#'

bonfInfinite <- function(d, alpha=0.05, alphai, random=TRUE,
                        date.format="%Y-%m-%d") {

    if(!(is.data.frame(d))){
        stop("d must be a dataframe.")
    }

    if(length(d$id) == 0){
        stop("The dataframe d is missing a column 'id' of identifiers.")
    } else if(length(d$pval) == 0){
        stop("The dataframe d is missing a column 'pval' of p-values.")
    }

    if(length(d$date) == 0){
        warning("No column of dates is provided, so p-values are treated as
        being ordered sequentially with no batches.")
        random <- FALSE
    } else if(any(is.na(as.Date(d$date, date.format)))){
        stop("One or more dates are not in the correct format.")
    } else {
        d <- d[order(as.Date(d$date, format = date.format)),]
    }

    if(alpha<0 || alpha>1){
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

    if(missing(alphai)){
        alphai <- 0.07720838*alpha*log(pmax(seq_len(N),2)) /
                    ((seq_len(N))*exp(sqrt(log(seq_len(N)))))
    } else if (any(alphai<0)){
        stop("All elements of alphai must be non-negative.")
    } else if(sum(alphai)>alpha){
        stop("The sum of the elements of alphai must not be greater than
        alpha.")
    }

    if(random){

        if(exists(".Random.seed", where = .GlobalEnv)){
            old.seed <- .Random.seed
            on.exit({ .Random.seed <<- old.seed })
        } else {
            on.exit({set.seed(NULL)})
        }

        set.seed(1)

        lst <- lapply(split(d, d$date),
                function(x){x[sample.int(length(x$date)),]})

        d <- do.call('rbind', lst)
        rownames(d) <- NULL
    }

    pval <- d$pval
    R <- as.numeric(pval <= alphai)
    d.out <- data.frame(d, alphai, R)

    return(d.out)
}
