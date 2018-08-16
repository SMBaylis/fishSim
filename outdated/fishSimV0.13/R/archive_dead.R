#' archive_dead(): takes all the dead individuals and copies them to an archive matrix.
#'
#' For larger simulations, the matrix 'indiv' may grow very large and slow down the simulation.
#' In these cases, run-times may be improved by periodically moving dead individuals into an
#' archive that is read and written less frequently than 'indiv'.
#' @param indiv A matrix of individuals, as from makeFounders(), move(), mate(), or mort().
#' @param archive A matrix of individuals, probably from make_archive() or a previous call of
#'                archive_dead().
#' @export
#' @examples
#' archive <- make_archive()
#' ages <- min(as.numeric(indiv[,8])):max(as.numeric(indiv[,8]))
#' ageMort <- 0.1 + (0.2*1/(ages+1))  ## placeholder ageMort with (extreme) negative senescence
#' stocks <- c(0.3,0.3,0.4)  ## matches defaults in makeFounders
#' admix.m <- matrix(NA, nrow = length(stocks), ncol = length(stocks))
#' for (i in 1:nrow(admix.m)) {
#'      admix.m[i,] <- stocks*stocks[i]
#' } 
#' indiv <- makeFounders()
#' for(k in 1:30) {
#'    indiv <- move(indiv = indiv, moveMat = admix.m)
#'    indiv <- mate(indiv = indiv, osr = c(0.55,0.45), year = k)
#'    indiv <- mort(indiv = indiv, type = "age", ageMort = ageMort, year = k)
#'    indiv <- birthdays(indiv = indiv)
#'    archive <- archive_dead(indiv = indiv, archive = archive) # archives a copy of dead animals'
#'                                                              # data
#'    indiv <- remove_dead(indiv = indiv) # actually removes the dead from 'indiv'. 
#' }

archive_dead <- function(indiv = mort(), archive = make_archive() ) {
    archive <- rbind(archive, indiv[!is.na(indiv[,6]),] )
    return(archive)
}

