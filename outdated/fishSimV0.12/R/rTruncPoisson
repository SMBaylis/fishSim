#' rTruncPoisson(): draw from a zero-truncated Poisson distribution.
#'
#' Returns a vector of draws from a zero-truncated Poisson distribution, with specified
#' lambda.
#'
#' @param n The number of values to draw
#' @param T The value of lambda for an equivalent non-truncated Poisson
#' @export

rTruncPoisson <- function(n = 1, T = 0.5) {  ## sample from a zero-truncated Poisson dist.
    U <- runif(n)
    t <- -log(1-U*(1-exp(-T)))
    T1 <- (T-t)
    X <- rpois(n,T1)+1
    return(X)
} ## thanks to Ted Harding, https://stat.ethz.ch/pipermail/r-help/2005-May/070678.html

