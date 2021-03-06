% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setBound.R
\name{setBound}
\alias{setBound}
\title{setBound}
\usage{
setBound(alg, alpha = 0.05, N)
}
\arguments{
\item{alg}{A string that takes the value of one of the following: LOND, LORD,
LORDdep, SAFFRON, ADDIS, LONDstar, LORDstar, SAFFRONstar, or
Alpha_investing}

\item{alpha}{Overall significance level of the FDR procedure, the default is
0.05. The bounds for LOND and LORDdep depend on alpha.}

\item{N}{An upper bound on the number of hypotheses to be tested}
}
\value{
\item{bound}{ A vector giving the values of a default sequence
\eqn{\gamma_i} of nonnegative numbers.}
}
\description{
Calculates a default sequence of non-negative numbers \eqn{\gamma_i} that sum
to 1, given an upper bound \eqn{N} on the number of hypotheses to be tested.
}
