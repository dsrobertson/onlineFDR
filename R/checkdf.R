checkdf <- function(d, random, date.format) {
    
    if (length(d$pval) == 0) {
        stop("The dataframe d is missing a column 'pval' of p-values.")
    }
    
    if (length(d$date) == 0) {
        # warning('No column of dates is provided, so p-values are treated as being
        # ordered sequentially with no batches.')
        random = FALSE
    } else if (any(is.na(as.Date(d$date, date.format)))) {
        stop("One or more dates are not in the correct format.")
    } else {
        d <- d[order(as.Date(d$date, format = date.format)), ]
    }
    
    if (random) {
        d <- randBatch(d)
    }
    
    return(d)
}
TRUE
TRUE
