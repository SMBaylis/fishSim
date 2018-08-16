#' make_archive(): sets up an archive matrix, for storing simulation outputs.
#'
#' For larger simulations, the matrix 'indiv' may grow very large and slow down the simulation.
#' In these cases, run-times may be improved by periodically moving dead individuals into an
#' archive that is read and written less frequently than 'indiv'.
#' @return The function returns an empty 0-by-8 matrix, to which subsets of 'indiv' can be
#'         attached. See archive_dead() and remove_dead().
#' @export

make_archive <- function() {
    archive <- matrix(data = NA, nrow = 0, ncol = 8)
    return(archive)
}

