checkStarVersion <- function(d, N, version, batch.sizes) {
    
    if (!(version %in% c("async", "dep", "batch"))) {
        stop("version must be 'async', 'dep' or 'batch'.")
    }
    
    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) abs(x - round(x)) < 
        tol
    
    if (version == "async") {
        
        if (any(is.na(d$decision.times))) {
            stop("Please provide a decision time for all p-values.")
        }
        
        if (any(!is.wholenumber(d$decision.times))) {
            stop("All decision times should be integers.")
        }
        
        if (any(d$decision.times <= 0)) {
            stop("All decision times should be positive integers.")
        }
        
        version <- 1
        
    } else if (version == "dep") {
        
        if (any(is.na(d$lags))) {
            stop("Please provide lags for all p-values.")
        }
        
        if (N > 1 && any(d$lags[seq_len(N - 1) + 1] > d$lags[seq_len(N - 1)] + 1)) {
            stop("The sequence of lags must satisfy L_(t+1) <= L_t + 1")
        }
        
        if (any(!is.wholenumber(d$lags))) {
            stop("All lags should be integers.")
        }
        
        if (any(d$lags < 0)) {
            stop("All lags should be positive integers.")
        }
        
        version <- 2
        
    } else if (version == "batch") {
        
        if (sum(batch.sizes) != N) {
            stop("The sum of the batch sizes must equal the number of p-values observed.")
        }
        
        if (any(!is.wholenumber(batch.sizes))) {
            stop("All batch sizes should be integers.")
        }
        
        if (any(batch.sizes <= 0)) {
            stop("All batch sizes should be positive integers.")
        }
        
        version <- 3
    }
    
    return(version)
}
TRUE
TRUE
