checkPval <- function(d) {
    
    if (is.data.frame(d)) {
        if (any(is.na(d$pval))) {
            warning("Your data contains missing p-values. Missing p-values were omitted.")
            d <- na.omit(d)
        }
        
        if (!(is.numeric(d$pval))) {
            stop("The vector of p-values contain at least one non-numeric element.")
        } else if (any(d$pval > 1 | d$pval < 0)) {
            stop("All p-values must be between 0 and 1.")
        } 
        
    } else if (is.vector(d)) {
        
        if (any(is.na(d))) {
            warning("Your data contains missing p-values. Missing p-values were omitted.")
            d <- d[!is.na(d)]
        }
        if (!(is.numeric(d))) {
            stop("The vector of p-values contain at least one non-numeric element.")
        } else if (any(d > 1 | d < 0)) {
            stop("All p-values must be between 0 and 1.")
        }
    }
    return(d)
}
TRUE
TRUE
