#' mort(): sets death year to a constant for some members of the population.
#' 
#' Members are chosen according to one of several defined mortality structures.
#' Mortality rates can bring the population to a specified size, or can be a flat probability,
#' or a probability that depends on age, stock, or age:stock.
#'
#' @param indiv A matrix of inidviduals, as from makeFounders(), move(), or mate().
#' @param type One of "simple", "flat", "age", "stock", or "ageStock".
#'             If type = "simple" and the living pop > maxPop, individuals are killed at random until
#'             the living pop == maxPop. Can easily be set to never cause extinctions.
#'             If type = "flat", individuals are killed with probability set in mortRate. Generates
#'             an exponential death curve.
#'             If type = "age", individuals are killed with probability for their age set in ageMort.
#'             If type = "stock", individuals are killed with probability for their stock set in
#'             stockMort.
#'             If type = "ageStock", individuals are killed with probability for their age:stock
#'             combination set in ageStockMort.
#' @param maxAge Sets an age at which animals *will* be killed before anything else happens. Allows
#'               a short age-specific mortality curve to be set, without checking if there are any
#'               individuals outside the range for each iteration.
#' @param maxPop If type = "simple", the population will be reduced to this number, if not already
#'               smaller. See 'type'.
#' @param mortRate Numeric mortality rate. See 'type'.
#' @param ageMort Numeric vector of mortality rates, one for each age, ordered 0:max(age). See 'type'.
#' @param stockMort Numeric vector of mortality rates, one for each stock, ordered 1:max(stock). Note
#'                  that stocks are numbered (as in makeFounders() ), not named. Because stocks are
#'                  stored as a character vector, stocks are converted via as.numeric() to associate
#'                  rates with stocks. This distinction is important in cases with >9 stocks. See
#'                  'type'.
#' @param ageStockMort A matrix of mortality rates, with age by row and stock by column. See 'ageMort'
#'                     and 'stockMort' for structure of rows and columns.
#' @export

mort <- function(indiv = makeFounders(), year = "-1", type = "simple", maxAge = Inf,
                 maxPop = 1000, mortRate, ageMort, stockMort, ageStockMort) {

    if (!(type %in% c("simple", "flat", "age", "stock", "ageStock"))) {
        stop("selected 'type' must be 'simple', 'flat', 'age', 'stock', or 'ageStock'")
    }

    geriatrics <- is.na(indiv[,6]) & as.numeric(indiv[,8]) > maxAge
    indiv[,6][geriatrics] <- year  ## marks geriatrics as dead this year, unless they
                                   ## already have a non-NA death year.
    if(type == "simple") {
        toKill <- nrow(indiv) - (maxPop + sum(!is.na(indiv[,6])))
        stillAlive <- nrow(indiv[is.na(indiv[,6]),])
        if(toKill > 0) {
            indiv[is.na(indiv[,6]),][sample.int(stillAlive, size = toKill),6] <- year
            ## magic one-liner takes the rows of indiv where death hasn't happened
            ## (i.e., indiv[is.na(indiv[,6]),] ), and within those rows, randomly selects
            ## a number of rows (the number is 'toKill') and assigns the current year as
            ## the death year for those rows.
        }
    } else if (type == "flat") {
        stillAlive <- nrow(indiv[is.na(indiv[,6]),])
        indiv[is.na(indiv[,6]),][runif(n=stillAlive, min=0, max=1)<mortRate, 6] <- year
        ## magic one-liner assigns a random number to every row where death hasn't happened
        ## and assigns the current year as the death year iff the random number is less than
        ## mortRate.
    } else if (type == "age") {
        stillAlive <- nrow(indiv[is.na(indiv[,6]),])
        for(i in min(as.numeric(indiv[,8])):max(as.numeric(indiv[,8]))) {
            ageProbs <- runif(n = nrow(indiv[is.na(indiv[,6])&as.numeric(indiv[,8])==i,,drop = FALSE]))
            ## for the individuals at age i, make a vector of runif numbers between 0 and 1
            if(length(ageProbs) > 1) {
                indiv[is.na(indiv[,6])&as.numeric(indiv[,8])==i,][ageProbs<ageMort[i+1],6] <- year
                ## if there are any individuals, compare the random number to the mortality rate for
                ## that age-class, and assign a death year if random number < mortality rate.
            } else if(length(ageProbs) == 1) {
                if(ageProbs<ageMort[i]) {
                    indiv[is.na(indiv[,6])&as.numeric(indiv[,8])==i,6] <- year
                }
            }  ## special error-handling for age-classes with only one individual.
        }
    } else if (type == "stock") {
        stillAlive <- nrow(indiv[is.na(indiv[,6]),])
        for(j in min(as.numeric(indiv[,7])):max(as.numeric(indiv[,7]))) {
            stockProbs <- runif(n=nrow(indiv[is.na(indiv[,6])&as.numeric(indiv[,7])==j,,drop=FALSE]))

            if(length(stockProbs) > 1) {
                indiv[is.na(indiv[,6])&as.numeric(indiv[,7])==j,][stockProbs<stockMort[j],6] <- year
            } else if(length(stockProbs == 1)) {
                if(stockProbs<stockMort[j]) {
                    indiv[is.na(indiv[,6])&as.numeric(indiv[,7])==j,6] <- year
                }
            }  ## almost exactly the same as age-specific mort, but no 'stock 0' synonymous with
               ## 'age 0'.
        }
    } else if (type == "ageStock") {
        stillAlive <- nrow(indiv[is.na(indiv[,6]),])
        for (i in min(as.numeric(indiv[,8])):max(as.numeric(indiv[,8]))) {
            for (j in min(as.numeric(indiv[,7])):max(as.numeric(indiv[,7]))) {
                ageStockProbs <- runif(n = nrow(indiv[is.na(indiv[,6])&
                                                      as.numeric(indiv[,8])==i&
                                                      as.numeric(indiv[,7])==j,,drop = FALSE]))
                if(length(ageStockProbs) > 1) {
                    indiv[is.na(indiv[,6])&
                          as.numeric(indiv[,8])==i&
                          as.numeric(indiv[,7])==j,][ageStockProbs<ageStockMort[i+1,j],6] <- year
                } else if(length(ageStockProbs) == 1) {
                    if(ageStockProbs<ageStockMort[i,j]) {
                        indiv[is.na(indiv[,6])&
                              as.numeric(indiv[,8])==i&
                              as.numeric(indiv[,7])==j,6] <- year
                    }
                }
            }
        } ## the same structure as 'age' and 'stock', but iterating through both.     
    }
    return(indiv)
}

