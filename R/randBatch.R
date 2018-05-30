randBatch <- function(d, seed) {

    if(!(is.null(seed))){
        set.seed(seed)
    }

    lst <- lapply(split(seq_len(nrow(d)), d$date), function(x) {
        x[sample.int(length(x))]
    })

    d <- d[unlist(lst, use.names = FALSE),]
    rownames(d) <- NULL

    return(d)

}
