checkPval <- function(d) {
    
    if (sum(is.na(d$pval)) > 0) {
        warning("Your data contains missing values. Missing values were omitted.")
        d <- na.omit(d)
    }
    
    if (!(is.numeric(d$pval))) {
        stop("The vector of p-values contain at least one non-numeric element.")
    } else if (any(d$pval > 1 | d$pval < 0)) {
        stop("All p-values must be between 0 and 1.")
    }
    
    return(d)
}
TRUE
TRUE
