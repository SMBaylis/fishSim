#' birthdays(): takes an individual-data matrix, and adds one to each individual's age
#'
#' Nothing fancy: this is separated from the other functions to allow more flexible
#' assignments of movement, mating and mortality within each year.
#'
#' @param indiv A matrix of individuals, as from makeFounders(), move(), mate(), or mort().
#' @export

birthdays <- function(indiv = makeFounders() ) {
    indiv[,8] <- as.numeric(indiv[,8]) + 1
    return(indiv)
}

