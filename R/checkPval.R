checkPval <- function(pval) {
    
    if (anyNA(pval)) {
        warning("Missing p-values were ignored.")
        pval <- stats::na.omit(pval)
    }
    
    if (!(is.numeric(pval))) {
        stop("The vector of p-values contain at least one non-numeric element.")
    } else if (any(pval > 1 | pval < 0)) {
        stop("All p-values must be between 0 and 1.")
    }
}
TRUE
TRUE
