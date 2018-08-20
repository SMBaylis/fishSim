#' makeFounders(): makes a matrix of individuals to act as a founding population.
#'
#' Returns a pop-by-8 character matrix, with each row being an individual in the
#' founder population.
#' [,1] is a unique (uuid) identifier for each animal.
#' [,2] is sex, values "M" or "F".
#' [,3] is "founder" for all animals.
#' [,4] is "founder" for all animals.
#' [,5] is "-1" for all animals.
#' [,6] is NA for all animals.
#' [,7] is the stock membership for each animal.
#' [,8] is the age of each animal (in 'breeding seasons') at year 0.
#'
#' @param pop The size of the founder population.
#' @param osr A numeric vector describing the sex ratio, c([male], [female]). Must sum to 1.
#' @param stocks A numeric vector describing the probability that an individual
#'               is in each stock. Must sum to 1.
#' @param maxAge Numeric. The max age to which an animal may survive.
#' @param survCurv Numeric vector. Describes the probability within the founder cohort of belonging
#'                 to each age-class for age=-classes 1:maxAge. Must sum to 1. Cannot be blank.
#' @export

makeFounders <- function(pop = 1000, osr = c(0.5,0.5), stocks = c(0.3,0.3,0.4),
                         maxAge = 20, survCurv = 0.7^(1:maxAge)/sum(0.7^(1:maxAge))) {
    require(ids)
    
    if(sum(osr) != 1) warning("osr does not sum to 1")
    if(sum(stocks) != 1) warning("stocks do not sum to 1")
    if(sum(survCurv) != 1) warning("survCurv does not sum to 1")
    if(length(survCurv) != maxAge) warning("survCurv and maxAge imply different maximum ages")
    
    indiv <- matrix(data = NA ,nrow = pop, ncol = 8)
    indiv[,1] <- uuid(n = nrow(indiv), drop_hyphens=TRUE)  ## uuid IDs for each animal.
    indiv[,2] <- sample(c("M", "F"), nrow(indiv), TRUE, prob = osr)
                                        # assign sexes by probability
    indiv[,3] <- c(rep("founder", nrow(indiv))) ## has no father
    indiv[,4] <- c(rep("founder", nrow(indiv))) ## has no mother
    indiv[,5] <- c(rep(-1, nrow(indiv))) ## founders are born at year -1
    indiv[,6] <- c(rep(NA, nrow(indiv))) ## founders are not yet dead
    indiv[,7] <- sample(1:length(stocks), nrow(indiv), TRUE, prob = stocks)
                                        # assign founder stock membership
    indiv[,8] <- sample.int(maxAge, nrow(indiv), TRUE, prob = survCurv)

    return(indiv)
}

###############################################################################################

## Section 2: move, mate, and mort.
## Section 2a: movement, mating and births

#stocks <- c(0.3,0.3,0.4)  ## matches defaults in makeFounders
#admix.m <- matrix(NA, nrow = length(stocks), ncol = length(stocks))
#for (i in 1:nrow(admix.m)) {
#    admix.m[i,] <- stocks*stocks[i]
#}  ## probability of moving into a stock is proportional to the default size of that stock
   ## from makeFounders.

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

#############################################################################################

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

###############################################################################################

#' mate(): find male-female pairs within the same stock, and generate batches of offspring
#' until a recruitment quota is filled.
#'
#' Returns an individual matrix with added newborns at age 0.
#' 
#' @param indiv A matrix of individuals, e.g., generated in makeFounders() or output from move()
#' @param fecundity Numeric variable. mean number of recruits to generate, as a proportion of
#'                  the number of animals in 'indiv'. nrow(indiv)*fecundity is the annual
#'                  turnover - i.e., the number of animals killed in 'mortality' will equal
#'                  nrow(indiv)*fecundity in order to keep the population size constant.
#' @param batchSize Numeric. The mean number of offspring produced per pairing. Follows ~Poisson.
#'                  Used iff type = "flat".
#' @param osr Numeric vector with length two, c(male, female), giving the sex ratio at birth
#'            (recruitment). Used to assign sexes to new offspring.
#' @param year Intended to be used in a simulation loop - this will be the iteration number, and
#'             holds the 'birthyear' value to give to new recruits.
#' @param firstBreed Integer variable. The age at first breeding, default zero. The minimum age
#'                   at which individuals can breed. Applies to potential mothers and potential
#'                   fathers. Both firstBreed and 'fecundityCurve', 'maleCurve, and 'femaleCurve' are
#'                   capable of specifying an age at first breeding, and firstBreed takes precedence.
#' @param type The type of fecundity-age relationship to simulate. Must be one of "flat",
#'             "age", or "ageSex". If "flat", offspring batch sizes will be the same for
#'             all age:sex combinations. If "age", the number of offspring for each pairing
#'             is a function of 'fecundityCurve', with the less-fecund parent determining the
#'             batch size. If "ageSex", the number of offspring for each pairing is a function
#'             of 'maleCurve' and 'femaleCurve', with the less-fecund parent determining the
#'             batch size.
#' @param fecundityCurve Numeric vector describing the age-specific fecundity curve. One
#'                       value per age, over all ages from 0:max(indiv[,8]). Used if "type"
#'                       = "age". Note that 'firstBreed' can interfere with 'fecundityCurve'
#'                       by setting fecundities to zero for some age classes. Recommended
#'                       usage is to set 'firstBreed' to zero whenever 'fecundityCurve' is
#'                       specified.
#' @param maleCurve Numeric vector describing age-specific fecundity for males. One value per
#'                  age, over all ages from 0:max(indiv[,8]). Used if "type" = "ageSex".
#'                  Note that 'firstBreed' can interfere with 'fecundityCurve' by setting
#'                  fecundities to zero for some age classes. Recommended usage is to set
#'                  'firstBreed' to zero whenever 'fecundityCurve' is specified.
#' @param femaleCurve Numeric vector describing age-specific fecundity for females. One value
#'                    per age, over all ages from 0:max(indiv[,8]). Used if "type" = "ageSex".
#'                    Note that 'firstBreed' can interfere with 'fecundityCurve' by setting
#'                    fecundities to zero for some age classes. Recommended usage is to set
#'                    'firstBreed' to zero whenever 'fecundityCurve' is specified.
#' @export
 
mate <- function(indiv = makeFounders(), fecundity = 0.2, batchSize = 0.5, osr = c(0.5,0.5),
                   year = "-1", firstBreed = 0, type = "flat",
                   fecundityCurve, maleCurve, femaleCurve) {
    require(ids)
    if (!(type %in% c("flat", "age", "ageSex"))) {
        stop("'type' must be one of 'flat', 'age', or 'ageSex'.")
    }
    
    sprog.m <- matrix(data = NA, nrow = ceiling(nrow(indiv)*fecundity),
                      ncol = 8)
    mothers <- subset(indiv, indiv[,2] == "F")
    if(nrow(mothers) == 0) stop("There are no females in the population")
    fathers <- subset(indiv, indiv[,2] == "M")
    if(nrow(fathers) == 0) stop("There are no males in the population")
    ticker <- 1
    while(ticker <= nrow(sprog.m)) {
        drawMother <- mothers[sample(nrow(mothers), size = 1, replace = FALSE),]
                                        # Note: character vector, not matrix
        fathersInStock <- subset(fathers, fathers[,7] == drawMother[7])
        if(nrow(fathersInStock) > 1) {
            drawFather <- fathersInStock[sample(nrow(fathersInStock), size = 1, replace = FALSE),]
            if(drawMother[8] >= firstBreed & drawFather[8] >= firstBreed) {
                if(type == "flat") {
                    n.sprogs <- rTruncPoisson(n = 1, T = batchSize) ## age-dependent fecundity here?
                } else if (type == "age") {
                    n.sprogs <- rpois(n = 1,
                                           lambda = min(c(fecundityCurve[as.numeric(drawFather[8])+1],
                                                          fecundityCurve[as.numeric(drawMother[8])+1]))
                                      )
                } else if (type == "ageSex") {
                    n.sprogs <- rpois(n = 1,
                                      lambda = min(c(maleCurve[as.numeric(drawFather[8])+1],
                                                     femaleCurve[as.numeric(drawMother[8])+1]))
                                      )
                }
                batch <- matrix(data = NA, nrow = n.sprogs, ncol = ncol(indiv))
                batch[,1] <- uuid(n = nrow(batch), drop_hyphens = TRUE)
                batch[,2] <- sample(c("M", "F"), size = nrow(batch), replace = TRUE, prob = osr)
                batch[,3] <- drawFather[1]
                batch[,4] <- drawMother[1]
                batch[,5] <- year
                ## [,6] stays 'NA', because no death year yet.
                batch[,7] <- drawMother[7] ## always recruits into the parents' stock.
                batch[,8] <- 0 ## it's a zero-year-old until the first 'birthday' step.

                if(n.sprogs > 0) {
                    if ((ticker + n.sprogs) <= nrow(sprog.m)) {
                        sprog.m[ticker:(ticker+n.sprogs-1),] <- batch
                    } else if ((ticker + n.sprogs) > nrow(sprog.m)) {
                        sprog.m[ticker:nrow(sprog.m),] <- batch[1:(nrow(sprog.m)-(ticker-1)),1:8]
                    }
                }
            } else {
                n.sprogs <- 0
            }
        }
        ticker <- ticker+n.sprogs
    }
    indiv <- rbind(indiv, sprog.m)
    return(indiv)
}







###################################################################################################

## ages <- min(as.numeric(indiv[,8])):max(as.numeric(indiv[,8]))
## ageMort <- 0.1 + (0.2*1/(ages+1))  ## placeholder ageMort with (extreme) negative senescence

## stocks <- min(as.numeric(indiv[,7])):max(as.numeric(indiv[,7]))
## stockMort <- sample(c(0.1,0.5,0.2), length(stocks), replace = TRUE)
                                        # placeholder stockMort with variable mortality

## ageStockMort <- matrix(data = ageMort, nrow = length(ages), ncol = length(stocks)) 
## for(i in 1:ncol(ageStockMort)) ageStockMort[,i] <- ageStockMort[,i]*sample(c(1,2,1.3),1)
## placeholder ageStockMort with a simple multiplicative effect of stock on age-related mort.
## only works with these values because ageMort is never >= 0.5.

## Section 2b: mortality and aging.

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

#############################################################################################

#' sexSwitch(): induces sex-switching in a subset of the population.
#'
#' Given an individual-data matrix, this function stochastically switches the sex of individuals
#' with a set probability. Can be one-directional (i.e., only male-to-female or only female-to-male
#' switches), or bidirectional.
#'
#' @param indiv A matrix of individuals, as from makeFounders(), mate(), or mort().
#' @param direction One of "MF", "FM", or "both", indicating male-to-female, female-to-male,
#'                  or both directions.
#' @param prob A numeric switching probability.
#' @export

sexSwitch <- function(indiv = makeFounders(), direction = "both", prob = 0.0001) {
    if(!(direction %in% c("MF", "FM", "both"))) {
        stop("'direction' must be one of 'MF', 'FM', or 'both'.")
    }
    probs <- runif(nrow(indiv))
    if(direction == "both") {
        indiv[probs < prob & indiv[,2] == "M",2] <- "F" ## males to female
        indiv[probs < prob & indiv[,2] == "F",2] <- "M" ## females to male
    } else if (direction == "MF") {
        indiv[probs < prob & indiv[,2] == "M",2] <- "F" ## males to female
    } else if(direction == "FM") {
        indiv[probs < prob & indiv[,2] == "F",2] <- "M" ## females to male
    }
    return(indiv)
}

################################################################################################

#' birthdays(): takes an individual-data matrix, and adds one to each individual's age
#'
#' Nothing fancy: this is separated from the other functions to allow more flexible
#' assignments of movement, mating and mortality within each year. Only updates ages
#' for animals that are alive, so indiv[,8] will remain 'age at death' for all dead animals.
#'
#' @param indiv A matrix of individuals, as from makeFounders(), move(), mate(), or mort().
#' @export

birthdays <- function(indiv = makeFounders() ) {
    indiv[is.na(indiv[,6]) ,8] <- as.numeric(indiv[is.na(indiv[,6]),8]) + 1
    return(indiv)
}

###############################################################################################

## Section 2c: archiving dead animals to keep the active data matrix small.

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

################################################################################################

## an example implementation:
#start_time <- Sys.time()
## 1) setup:
#archive <- make_archive()

#ages <- min(as.numeric(indiv[,8])):max(as.numeric(indiv[,8]))
#ageMort <- 0.1 + (0.2*1/(ages+1))  ## placeholder ageMort with (extreme) negative senescence

#stocks <- c(0.3,0.3,0.4)  ## matches defaults in makeFounders
#admix.m <- matrix(NA, nrow = length(stocks), ncol = length(stocks))
#for (i in 1:nrow(admix.m)) {
#    admix.m[i,] <- stocks*stocks[i]
#} 

#indiv <- makeFounders(pop = 1000, osr = c(0.55, 0.45))

## 2) population processes and archiving
#for(k in 1:30) {
#    indiv <- move(indiv = indiv, moveMat = admix.m)
#    indiv <- mate(indiv = indiv, osr = c(0.55,0.45), year = k)
#    indiv <- mort(indiv = indiv, type = "age", ageMort = ageMort, year = k)
#    indiv <- birthdays(indiv = indiv)

#    archive <- archive_dead(indiv = indiv, archive = archive)
#    indiv <- remove_dead(indiv = indiv)
#}

## 3) archiving the living; sim close
#indiv[,6] <- c(rep("alive", nrow(indiv)))
#archive <- archive_dead(indiv = indiv, archive = archive)
#
#end_time <- Sys.time()
#end_time - start_time
#nrow(archive)
