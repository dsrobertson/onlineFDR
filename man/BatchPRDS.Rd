% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BatchPRDS.R
\name{BatchPRDS}
\alias{BatchPRDS}
\title{BatchPRDS: Online batch FDR control under Positive Dependence}
\usage{
BatchPRDS(d, alpha = 0.05, gammai)
}
\arguments{
\item{d}{A dataframe with three columns: an identifier (`id') 
and a p-value (`pval').}

\item{alpha}{Overall significance level of the FDR procedure, the default is
0.05.}

\item{gammai}{Optional vector of \eqn{\gamma_i}. A default is provided with
\eqn{\gamma_j} proportional to \eqn{1/j^(1.6)}.}
}
\value{
\item{out}{ A dataframe with the original data \code{d} and the indicator
function of discoveries \code{R}. Hypothesis \eqn{i} is rejected if the
\eqn{i}-th p-value within the \eqn{t}-th batch is less than or equal to 
\eqn{(r/n)*\alpha_t}, where \eqn{r} is the rank of the \eqn{i}-th p-value 
within an ordered set and \eqn{n} is the total number of hypotheses 
within the \eqn{t}-th batch. As such, \code{R[i] = 1}  (otherwise \code{R[i] = 0}).}
}
\description{
Implements the BatchPRDS algorithm for online FDR control, where PRDS stands for
positive regression dependency on a subset, as presented by
Zrnic et. al. (2020).
}
\details{

}
\examples{

sample.df <- data.frame(
id = c('A15432', 'B90969', 'C18705', 'B49731', 'E99902',
    'C38292', 'A30619', 'D46627', 'E29198', 'A41418',
    'D51456', 'C88669', 'E03673', 'A63155', 'B66033'),
pval = c(2.90e-08, 0.06743, 0.01514, 0.08174, 0.00171,
        3.60e-05, 0.79149, 0.27201, 0.28295, 7.59e-08,
        0.69274, 0.30443, 0.00136, 0.72342, 0.54757),
batch = c(rep(1,5), rep(2,6), rep(3,4)))

BatchPRDS(sample.df)

}
\references{
Zrnic, T., Jiang D., Ramdas A. and Jordan M. (2020). The Power of Batching in Multiple Hypothesis Testing. \emph{International Conference on Artificial Intelligence and Statistics}: 3806-3815
}