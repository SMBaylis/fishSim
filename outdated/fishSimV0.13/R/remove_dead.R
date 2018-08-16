#' remove_dead(): removes the dead from an individual-data matrix.
#'
#' For larger simulations, the matrix 'indiv' may grow very large and slow down the simulation.
#' In these cases, run-times may be improved by periodically moving dead individuals into an
#' archive that is read and written less frequently than 'indiv'. See also archive_dead().
#' @param indiv A matrix of individuals, as from mort().
#' @export

remove_dead <- function(indiv = mort() ) {
    indiv <- indiv[is.na(indiv[,6]),]
    return(indiv)
}
