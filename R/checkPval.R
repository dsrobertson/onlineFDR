checkPval <- function(d) {
    
    if (is.data.frame(d)) {
        if (!("pval" %in% colnames(d))) {
            stop("The dataframe d is missing a column 'pval' of p-values.")
        }

        if (any(is.na(d$pval))) {
            # decision times and lags refer to positions in the sequence, so
            # rows cannot be dropped without changing their meaning
            if (any(c("decision.times", "lags") %in% colnames(d))) {
                stop("Missing p-values cannot be omitted when d contains ",
                     "decision times or lags.")
            }
            warning("Your data contains missing p-values. Missing p-values were omitted.")
            d <- d[!is.na(d$pval), , drop = FALSE]
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

    if (length(if (is.data.frame(d)) d$pval else d) == 0) {
        stop("d must contain at least one p-value.")
    }
    return(d)
}
