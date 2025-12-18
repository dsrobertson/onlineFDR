#' ADDIS_exhaustive: Exhaustive ADDIS procedure for online FDR control
#'
#' Implements an exhaustive variant of the ADDIS algorithm for online FDR
#' control by adapting code from the Fischer, L.: Exhaustive ADDIS procedures 
#' for online FWER control.
#'
#' @author Lasse Fischer
#'
#' @references Fischer, L.: Exhaustive ADDIS procedures for online FWER control.
#'   arXiv:2308.13827 <https://arxiv.org/abs/2308.13827>
#'
#' @param d Either a vector of p-values, or a dataframe with at least a
#'   `pval` column (and optionally `id`).
#' @param alpha Overall significance level of the procedure, default 0.05.
#' @param tau Optional threshold for hypotheses to be selected for testing.
#'   Must be between 0 and 1, defaults to 0.5.
#' @param lambda Optional parameter that sets the threshold for `candidate'
#'   hypotheses. Must be between 0 and tau, defaults to 0.25.
#' @param gamma Optional vector of initial weights. If `NULL` (the default),
#'   a decreasing sequence proportional to j^(-1.6) is used, as in ADDIS().
#'
#' @return A dataframe with the original p-values `pval`, the per-hypothesis
#'   testing levels `alphai`, and the indicator of discoveries `R`.
#'
#' @export
ADDIS_exhaustive <- function(d, alpha = 0.05, tau = 0.5, lambda = 0.25, gamma = NULL) {
	d <- checkPval(d)

	if (is.data.frame(d)) {
		pval <- d$pval
	} else if (is.vector(d)) {
		pval <- d
	} else {
		stop("d must either be a dataframe or a vector of p-values.")
	}

	if (alpha <= 0 || alpha > 1) {
		stop("alpha must be between 0 and 1.")
	}

	if (tau <= 0 || tau > 1) {
		stop("tau must be between 0 and 1.")
	}

	if (lambda <= 0 || lambda > tau) {
		stop("lambda must be between 0 and tau.")
	}

	res <- .e_addis_spending(pval = pval, alpha = alpha, tau = tau, lambda = lambda, gamma = gamma)

	out <- data.frame(pval = pval,
										alphai = res$alphai,
										R = as.numeric(res$R))

	if (is.data.frame(d) && !is.null(d$id)) {
		out$id <- d$id
	}

	out
}

# Core exhaustive ADDIS implementation adapted from E_ADDIS_Spending in
.e_addis_spending <- function(pval, alpha, tau, lambda, gamma = NULL) {
# n: In original code `n` is always `length(pval)``
    n <- length(pval)

	if (length(tau) == 1) {
		tau <- rep(tau, n)
	}
	if (length(lambda) == 1) {
		lambda <- rep(lambda, n)
	}

	if (length(tau) != n || length(lambda) != n || length(pval) != n) {
		stop("mismatching length between tau, lambda and pval")
	}

	# Default gamma sequence as in ADDIS: proportional to j^(-1.6)
	if (is.null(gamma)) {
		gamma <- 0.4374901658/(seq_len(n + 1)^(1.6))
	} else {
		if (any(gamma < 0)) {
			stop("All elements of gamma must be non-negative.")
		}
		if (length(gamma) < n + 1) {
			stop("gamma must have length at least n + 1.")
		}
	}

	t <- 1
	alpha_k <- alpha
	alphai <- numeric(n)

	for (i in seq_len(n)) {
		alphai[i] <- alpha * gamma[t] * (tau[i] - lambda[i]) / (1 - alpha_k)
		if (pval[i] <= tau[i] && pval[i] > lambda[i]) {
			alpha_k <- alpha_k - (alphai[i] * (1 - alpha_k) / (tau[i] - lambda[i]))
			t <- t + 1
		}
	}

	R <- as.integer(pval <= alphai)
	list(alphai = alphai, R = R)
}
