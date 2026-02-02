#' ADDIS-exhaustive: Exhaustive ADDIS-spending procedure for online FWER control
#'
#' Implements an exhaustive variant of the ADDIS-spending algorithm for online
#' FWER control, as presented by Fischer et al. (2023). The procedure is a 
#' uniform improvement of ADDIS-spending, and no other FWER controlling
#' procedure can enlarge the event of rejecting any hypothesis.
#' 
#' The function takes as its input either a vector of p-values, or a dataframe
#' with two columns: an identifier (`id') and p-value (`pval'). Given an overall
#' significance level \eqn{\alpha}, ADDIS-exhaustive depends on constants
#' \eqn{\lambda} and \eqn{\tau}, where \eqn{\lambda < \tau}. Here \eqn{\tau \in
#' (0,1)} represents the threshold for a hypothesis to be selected for testing:
#' p-values greater than \eqn{\tau} are implicitly `discarded' by the procedure,
#' while \eqn{\lambda \in (0,1)} sets the threshold for a p-value to be a
#' candidate for rejection: ADDIS-exhaustive will never reject a p-value larger
#' than \eqn{\lambda}. The algorithms also require a sequence of non-negative
#' non-increasing numbers \eqn{\gamma_i} that sum to 1.
#'
#' The ADDIS-exhaustive procedure provably controls the FWER in the strong sense
#' for independent p-values.
#' 
#'
#' @author Lasse Fischer
#'
#' @references Fischer, L., Bofill Roig, M. and Brannath W. (2024). An
#' exhaustive ADDIS principle for online FWER control.
#' \emph{Biometrical Journal} 66(3) 2300237.
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
