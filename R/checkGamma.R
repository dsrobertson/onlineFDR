sumExceeds <- function(x, bound) {

    # Allows for rounding error in sum(), e.g. for the output of setBound()
    sum(x) > bound * (1 + .Machine$double.eps * length(x))
}

isNonIncreasing <- function(x) {

    n <- length(x)
    if (n < 2) {
        return(TRUE)
    }
    all(diff(x) <= 4 * .Machine$double.eps * abs(x[-n]))
}
