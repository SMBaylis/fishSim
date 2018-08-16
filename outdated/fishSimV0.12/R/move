#' move(): markovian movement between breeding stocks.
#'
#' returns a pop-by-8 character matrix, defined in the makeFounders documentation.
#'
#' @param indiv Individual matrix, e.g. from makeFounders(). Can also be a non-founder matrix.
#' @param moveMat An s-by-s matrix describing the probability of moving from each stock to
#'                each other stock, with 'from' by row and 'to' by column.
#' @export

move <- function(indiv = makeFounders(), moveMat = admix.m) {
    if(nrow(moveMat) != ncol(moveMat)) warning("movement matrix is not square")
    if(nrow(moveMat) > length(unique(indiv[,7]))) warning("One or more stocks start empty")
    
    newStock <- matrix(data = NA, nrow = nrow(indiv))
    for(i in 1:length(newStock)) {
        newStock[i] <- sample(1:nrow(moveMat), 1, TRUE,
                              prob = moveMat[as.numeric(indiv[i,7]),])
    }
    indiv[,7] <- newStock
    return(indiv)
}
