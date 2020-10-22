checkSTARdf <- function(d, version) {
    switch(as.character(version), async = {
        if (!("decision.times" %in% colnames(d))) {
            stop("d does not contain the column 'decision.times'.")
        }
        if (c("lags") %in% colnames(d)) {
            warning("d should not contain the columns 'lags'.")
        }
    }, dep = {
        if (!("lags" %in% colnames(d))) {
            stop("d does not contain the column 'lags'.")
        }
        if (c("decision.times") %in% colnames(d)) {
            warning("d should not contain the columns 'decision.times'.")
        }
    })
}

TRUE
TRUE
