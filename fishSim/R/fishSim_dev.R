#' makeFounders(): makes a matrix of individuals to act as a founding population.
#'
#' Returns a pop-by-8 character matrix, with each row being an individual in the
#' founder population.
#' [,1] is a unique (uuid) identifier for each animal.
#' [,2] is sex, values "M" or "F".
#' [,3] is "founder" for all animals.
#' [,4] is "founder" for all animals.
#' [,5] is the animal's birth year. Implicitly assumes that 'makeFounders' occurs at the
#'      very start of year 1, just after the 'birthdays' step of year 0.
#' [,6] is NA for all animals.
#' [,7] is the stock membership for each animal.
#' [,8] is the age of each animal (in 'breeding seasons') at the beginning of year 1,
#'      given that birthdays occur at the very end.
#' makeFounders() will throw a warning if osr, stocks, or survCurv do not sum to 1. It is
#' not strictly necessary that they sum to 1 (proportionality within each class is sufficient),
#' but error-checking and readability is easiest if they do sum to 1.
#'
#' @param pop The size of the founder population.
#' @param osr A numeric vector describing the sex ratio, c([male], [female]).
#' @param stocks A numeric vector describing the probability that an individual
#'               is in each stock.
#' @param maxAge Numeric. The max age to which an animal may survive.
#' @param survCurv Numeric vector. Describes the probability within the founder cohort of belonging
#'                 to each age-class for age=-classes 1:maxAge. Cannot be blank.
#' @seealso [fishSim::make_archive()]
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
##  indiv[,5] <- c(rep(-1, nrow(indiv))) ## founders are born at year -1
    indiv[,6] <- c(rep(NA, nrow(indiv))) ## founders are not yet dead
    indiv[,7] <- sample(1:length(stocks), nrow(indiv), TRUE, prob = stocks)
                                        # assign founder stock membership
    indiv[,8] <- sample.int(maxAge, nrow(indiv), TRUE, prob = survCurv)
    indiv[,5] <- 1 - as.numeric(indiv[,8]) ## back-infer birth year from age

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
#' returns a pop-by-8 character matrix, defined in the makeFounders documentation. Dev note:
#' this version will still apply movement to dead individuals, but most users would probably
#' expect the dead to stay in the stock where they died. Update that.
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
#' Returns an individual matrix with added newborns at age 0. This is the first
#' mate method, where the number of newborns is a set as a proportion of the adult population,
#' and matings happen between individuals that are older than the age at first-breeding
#' (optionally with breeding success-rates set by the age of the less-mature parent, etc.,) until
#' that 'quota' or newborns has been met, or all potential parents have reached breeding exhaustion.
#' A second mate method exists, where the number of newborns is derived from the fecundities of all
#' breeding females in the population - see altMate().
#' 
#' @param indiv A matrix of individuals, e.g., generated in makeFounders() or output from move()
#' @param fecundity Numeric variable. mean number of recruits to generate, as a proportion of
#'                  the number of animals in 'indiv'. nrow(indiv)*fecundity is the annual
#'                  turnover - i.e., the number of animals killed in 'mortality' will equal
#'                  nrow(indiv)*fecundity in order to keep the population size constant.
#' @param batchSize Numeric. The mean number of offspring produced per pairing. Follows ~Poisson.
#'                  Used iff type = "flat". Note that, if run within a loop that
#'                  goes [move -> mate -> mort -> birthdays], 'produced' is not the same as 'enters
#'                  age-class 1', as some mortality may occur at age 0.
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
#' @param maxClutch Numeric value giving the maximum clutch / litter / batch / whatever size.
#'                  Reduces larger clutches to this size, within each pairing.
#' @param exhaustMothers TRUE/FALSE value indicating whether mothers should become 'exhausted'
#'                         by one breeding attempt. If exhausted, an individual will not mate
#'                         again in this mate() call.
#' @param exhaustFathers TRUE/FALSE value indicating whether fathers should become 'exhausted'
#'                         by one breeding attempt. If exhausted, an individual will not mate
#'                         again in this mate() call.
#' @param fecundityCurve Numeric vector describing the age-specific fecundity curve. One
#'                       value per age, over all ages from 0:max(indiv[,8]). Used if "type"
#'                       = "age". Note that 'firstBreed' can interfere with 'fecundityCurve'
#'                       by setting fecundities to zero for some age classes. Recommended
#'                       usage is to set 'firstBreed' to zero whenever 'fecundityCurve' is
#'                       specified.
#' @param maleCurve Numeric vector describing age-specific fecundity for males. One value per
#'                  age, over all ages from 0:max(indiv[,8]). Used if "type" = "ageSex".
#'                  Note that 'firstBreed' can interfere with 'maleCurve' by setting
#'                  fecundities to zero for some age classes. Recommended usage is to set
#'                  'firstBreed' to zero whenever 'maleCurve' is specified.
#' @param femaleCurve Numeric vector describing age-specific fecundity for females. One value
#'                    per age, over all ages from 0:max(indiv[,8]). Used if "type" = "ageSex".
#'                    Note that 'firstBreed' can interfere with 'femaleCurve' by setting
#'                    fecundities to zero for some age classes. Recommended usage is to set
#'                    'firstBreed' to zero whenever 'femaleCurve' is specified.
#' @seealso [fishSim::altMate()]
#' @export

mate <- function(indiv = makeFounders(), fecundity = 0.2, batchSize = 0.5,
                 osr = c(0.5,0.5), year = "-1", firstBreed = 0, type = "flat", maxClutch = Inf,
                 exhaustMothers = FALSE, exhaustFathers = FALSE,
                 fecundityCurve, maleCurve, femaleCurve) {
    require(ids)
    if (!(type %in% c("flat", "age", "ageSex"))) {
        stop("'type' must be one of 'flat', 'age', or 'ageSex'.")
    }

    mothers <- subset(indiv, indiv[,2] == "F" & as.numeric(indiv[,8]) > firstBreed & is.na(indiv[,6]))
    if(nrow(mothers) == 0) stop("There are no females in the population")
    fathers <- subset(indiv, indiv[,2] == "M" & as.numeric(indiv[,8]) > firstBreed & is.na(indiv[,6]))
    if(nrow(fathers) == 0) stop("There are no males in the population")

    sprog.m <- matrix(data = NA, nrow = floor(nrow(indiv[is.na(indiv[,6]),])*fecundity),
                      ncol = 8)  ## Number of sprogs is a proportion of the number of
                                 ## *live* animals in the matrix.
    ticker <- 1
    while(ticker <= nrow(sprog.m)) {
        if(nrow(fathers) == 0) stop("All potential fathers are exhausted")
        if(nrow(mothers) == 0) stop("All potential mothers are exhausted")        
        drawMother <- mothers[sample(nrow(mothers), size = 1, replace = FALSE),]
                                        # Note: character vector, not matrix
        fathersInStock <- subset(fathers, fathers[,7] == drawMother[7])
        n.sprogs <- 0
        if(nrow(fathersInStock) > 1) {
            drawFather <- fathersInStock[sample(nrow(fathersInStock), size = 1, replace = FALSE),]
            if(as.numeric(drawMother[8]) >= firstBreed & as.numeric(drawFather[8]) >= firstBreed) {
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
                if(n.sprogs > maxClutch) n.sprogs <- maxClutch
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
            }
        }
        if(exhaustMothers == TRUE & n.sprogs > 0) {
            mothers <- mothers[mothers[,1] != drawMother[1] , , drop = FALSE]
        } ## remove exhausted mother from potential mothers.
        if(exhaustFathers == TRUE & n.sprogs > 0) {
            fathers <- fathers[fathers[,1] != drawFather[1] , , drop = FALSE]
        } ## remove exhausted father from potential fathers.
        ticker <- ticker+n.sprogs
    }
    indiv <- rbind(indiv, sprog.m)  
    return(indiv)
}


#' altMate(): breeding based on mature females, not quota-filling.
#'
#' Returns an individual matrix with added newborns at age 0. This is the second
#' mate method, where the number of offspring is derived from the number of mature females,
#' such that each mature female produces a number of offspring specified by a sampling distribution,
#' and fathers are randomly drawn from all mature males within the mother's stock. Note that in
#' this mating system, *maturity by age* is specified, rather than *fecundity by age*. A single
#' probability distribution sets the number of offspring for each female, but the probability
#' that an individidual female is mature may vary by age. The same maturity by age structure applies
#' for males. It is possible, in cases where the maturity-by-age slope is shallow, that an individual
#' may be 'mature' in one breeding season, but then 'not mature' the next season.
#' Another mate method exists, where the number of newborns is set as a proportion of the population
#' size, and mating occurs until the required number of offspring are generated (if possible, given
#' breeding constraints) - see mate().
#' 
#' @param indiv A matrix of individuals, e.g., generated in makeFounders() or output from move()
#' @param batchSize Numeric. The mean number of offspring produced per mature female. Follows
#'                  ~Poisson. Used for all maturity structures. Note that, if run within a loop that
#'                  goes [move -> altMate -> mort -> birthdays], 'produced' is not the same as
#'                  'enters age-class 1', as some individuals will die at age 0.
#' @param osr Numeric vector with length two, c(male, female), giving the sex ratio at birth
#'            (recruitment). Used to assign sexes to new offspring.
#' @param year Intended to be used in a simulation loop - this will be the iteration number, and
#'             holds the 'birthyear' value to give to new recruits.
#' @param firstBreed Integer variable. The age at first breeding, default zero. The minimum age
#'                   at which individuals can breed. Applies to potential mothers and potential
#'                   fathers. 'firstBreed', 'maturityCurve', 'maleCurve, and 'femaleCurve' are
#'                   all capable of specifying an age at first breeding, and 'firstBreed' takes
#'                   precedence.
#' @param type The type of maturity-age relationship to simulate. Must be one of "flat",
#'             "age", or "ageSex". If "flat", the probability of parenthood is the same for
#'             all age:sex combinations above firstBreed. If "age", the probability that an
#'             individual is sexually mature is age-specific, set in 'maturityCurve'. If "ageSex",
#'             the probability that an individual is sexually mature is age- and sex-specific,
#'             set for males in 'maleCurve' and for females in 'femaleCurve'.
#' @param maxClutch Numeric value giving the maximum clutch / litter / batch / whatever size.
#'                  Reduces larger clutches to this size, for each breeding female.
#' @param singlePaternity TRUE/FALSE value indicating whether all the offspring produced by
#'                        female in a year should have the same father. Default TRUE. If
#'                        FALSE, each offspring will have a randomly-drawn father from within
#'                        the mother's stock. Note that this can lead to rapid exhaustion of
#'                        fathers if exhaustFathers = TRUE.
#' @param exhaustFathers TRUE/FALSE value indicating whether fathers should become 'exhausted'
#'                       by one breeding attempt. If exhausted, an individual will only mate with
#'                       one female, though may father more than one offspring - see
#'                       'singlePaternity' and 'batchSize'.
#' @param maturityCurve Numeric vector describing the age-specific maturity curve. One
#'                      value per age, over all ages from 0:max(indiv[,8]). Used if "type"
#'                      = "age". Note that 'firstBreed' can interfere with 'maturityCurve'
#'                      by setting maturities to zero for some age classes. Recommended
#'                      usage is to set 'firstBreed' to zero whenever 'maturityCurve' is
#'                      specified.
#' @param maleCurve Numeric vector describing age-specific maturity for males. One value per
#'                  age, over all ages from 0:max(indiv[,8]). Used if "type" = "ageSex".
#'                  Note that 'firstBreed' can interfere with 'maleCurve' by setting
#'                  maturities to zero for some age classes. Recommended usage is to set
#'                  'firstBreed' to zero whenever 'maleCurve' is specified.
#' @param femaleCurve Numeric vector describing age-specific maturity for females. One value
#'                    per age, over all ages from 0:max(indiv[,8]). Used if "type" = "ageSex".
#'                    Note that 'firstBreed' can interfere with 'femaleCurve' by setting
#'                    maturities to zero for some age classes. Recommended usage is to set
#'                    'firstBreed' to zero whenever 'femaleCurve' is specified.
#' @seealso [fishSim::mate()]
#' @export

altMate <- function(indiv = makeFounders(), batchSize = 0.5,
                 osr = c(0.5,0.5), year = "-1", firstBreed = 0, type = "flat", maxClutch = Inf,
                 singlePaternity = TRUE, exhaustFathers = FALSE,
                 maturityCurve, maleCurve, femaleCurve) {
    require(ids)
    if (!(type %in% c("flat", "age", "ageSex"))) {
        stop("'type' must be one of 'flat', 'age', or 'ageSex'.")
    }
    mothers <- subset(indiv, indiv[,2] == "F" & as.numeric(indiv[,8]) > firstBreed & is.na(indiv[,6]))
    if(nrow(mothers) == 0) warning("There are no mature females in the population")
    fathers <- subset(indiv, indiv[,2] == "M" & as.numeric(indiv[,8]) > firstBreed & is.na(indiv[,6]))
    if(nrow(fathers) == 0) warning("There are no mature males in the population")

    if (type == "flat") {

        clutch <- rpois(n = nrow(mothers), lambda = batchSize)
        mothers <- subset(mothers, clutch > 0)  ## identify the mothers that truly breed,
        clutch <- clutch[clutch>0]              ## how many offspring each produces,
        
    } else if (type == "age") {

        mothers <- mothers[runif(nrow(mothers)) < maturityCurve[as.numeric(mothers[,8])],
                         , drop = FALSE] ## trims 'mothers' to those that pass a random
                                         ## maturity test.
        fathers <- fathers[runif(nrow(fathers)) < maturityCurve[as.numeric(fathers[,8])],
                         , drop = FALSE] ## trims 'fathers' to those that pass a random
                                         ## maturity test.
        
        clutch <- rpois(n = nrow(mothers), lambda = batchSize)
        mothers <- subset(mothers, clutch > 0) ## trims 'mothers' to those that truly breed
        clutch <- clutch[clutch>0]
        
    } else if (type == "ageSex") {

        mothers <- mothers[runif(nrow(mothers)) < femaleCurve[as.numeric(mothers[,8])] ,
                          ,drop = FALSE] ## trims 'mothers' to just those that pass a random
                                         ## maturity test.
        
        fathers <- fathers[runif(nrow(fathers)) < maleCurve[as.numeric(fathers[,8])] ,
                          ,drop = FALSE]  ## trims 'fathers' to just those that pass a random
                                          ## maturity test.
        clutch <- rpois(n = nrow(mothers), lambda = batchSize)
        mothers <- subset(mothers, clutch > 0)
        clutch <- clutch[clutch>0]
    }

    clutch[clutch > maxClutch] <- maxClutch ## delimits clutch sizes to not exceed maxClutch
    sprog.m <- matrix(data = NA, nrow = 0, ncol = 8) ## left empty if no-one breeds.
    
    for (s in unique(mothers[,7])) { ## s for 'stock'.
        mothersInStock <- mothers[mothers[,7] == s , , drop = FALSE]
        clutchInStock <- clutch[mothers[,7] == s]
        fathersInStock <- fathers[fathers[,7] == s , , drop = FALSE]
        if(nrow(fathersInStock) == 0) {
            warning (paste("There were no mature males in stock ",
                           s, ", so ", nrow(mothersInStock),
                           " mature females did not produce offspring",
                           sep = ""))
            sprog.stock <- matrix(data = NA, nrow = 0, ncol = 8)
        } else if(nrow(fathersInStock > 0)) {
            sprog.stock <- matrix(data = NA, nrow = sum(clutchInStock), ncol = 8)
            ticker <- 1
            for (m in 1:nrow(mothersInStock)) { ## m for 'mothers'
                if(nrow(fathersInStock) == 0) {
                    warning(paste("All fathers in stock ", s, " are exhausted.", sep = ""))
                } else {
                    sprog.stock[ticker:(ticker+clutchInStock[m]-1), 4] <- mothersInStock[m,1]
                    ## assign mother
                    if (singlePaternity == TRUE) {
                        sprog.stock[ticker:(ticker+clutchInStock[m]-1), 3] <-
                            fathersInStock[sample(1:nrow(fathersInStock), 1), 1]
                    } else if (singlePaternity == FALSE) {
                        if(nrow(fathersInStock) >= clutchInStock[m]) {
                            sprog.stock[ticker:(ticker+clutchInStock[m]-1), 3] <-
                                fathersInStock[sample(1:nrow(fathersInStock), clutchInStock[m]), 1]
                        } else {
                            sprog.stock[ticker:(ticker+nrow(fathersInStock)-1),3] <-
                                fathersInStock[ ,1]
                        }
                    } ## Assign father(s).
                    ## Note the potential conflict here with 'exhaustFathers' - but maybe not
                    ## of concern, because there can't be many species with multiple paternity
                    ## that also exhaust fathers after one mating attempt.
                    if(exhaustFathers == TRUE) {
                        fathersInStock <- fathersInStock[!fathersInStock[,1] %in% sprog.stock[,3],
                                                        ,drop = FALSE]
                        ## remove the used fathers from 'fathersInStock'                        
                    }
                    ticker <- ticker+clutchInStock[m]  ## increment ticker
                }
            }
        }
        sprog.stock <- sprog.stock[!is.na(sprog.stock[,3]), , drop = FALSE]
     
        sprog.stock[,1] <- uuid(n = nrow(sprog.stock), drop_hyphens = TRUE)
        sprog.stock[,2] <- sample(c("M", "F"), nrow(sprog.stock), TRUE, prob = osr)
        sprog.stock[,5] <- year
        sprog.stock[,7] <- s
        sprog.stock[,8] <- 0
        sprog.m <- rbind(sprog.m, sprog.stock)
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
#' @seealso [fishSim::archive_dead(), fishSim::remove_dead()]
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
#' @seealso [fishSim::remove_dead(), fishSim::make_archive()]
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
#' @seealso [fishSim::archive_dead(), fishSim::make_archive()]
#' @export

remove_dead <- function(indiv = mort() ) {
    indiv <- indiv[is.na(indiv[,6]), , drop = FALSE]
    return(indiv)
}

################################################################################################

#' check_growthrate(): estimate population growth under some mate() and mort() conditions.
#'
#' In complex models, population trends may not be immediately clear from the settings. Yet
#' it can be important to know population trends in advance: growing populations take progressively
#' more computing time, and the relationship dynamics between individuals across generations are
#' different for growing vs. shrinking populations. check_growthrate() provides an estimate of
#' the long-term population growth rate under some altMate() and mort() settings, assuming one
#' altMate() and one mort() per cycle. Estimation is via Leslie matrices.
#' check_growthrates only functions for some altMate() and mort() structures. Specifically, in
#' altMate, 'type' must be one of 'flat', 'age', or 'ageSex', and in mort, if 'type' is 'flat'
#' or 'age', a single growth rate will be returned, but if 'stock' or 'ageStock', one growth
#' rate will be returned per stock. If 'type' is 'simple' for mort(), the growthrate is forced to
#' zero, so check_growthrate() does not explicitly handle this case.
#' Some conditions can cause check_growthrate to fail or provide inaccurate estimates. If mature
#' females go unmated through lack of available fathers (for instance, if exhaustFathers = TRUE
#' in mate and N(mature females) > N(mature males) ), the Leslie matrix approach will provide an
#' over-estimate of the growth rate. In mate(), batchSize is the mean number of offspring per female,
#' but if maxClutch is also set to a value other than Inf, the *effective* mean number of offspring
#' per female is estimated by simulation. The mean number of female offspring per female per year
#' is assumed to be half of the effective mean number of offspring per female unless osr is specified,
#' in which case the proportion of female offspring is taken from osr.
#' If your model involves a variation not handled by check_growthrates(), you may find it simplest
#' to run your simulation for a few generations with a smallish founder population and empirically
#' estimate the growth rate.
#'
#' @param forceY1 optionally force first-year mortality to a specific value. Defaults to NA. If non-NA,
#'                should be a numeric value between 0 and 1. If NA, ignored. Must be the first
#'                argument so that uniroot()-based solutions for null growth will work.
#' @param mateType the value of 'type' used in the altMate() call. Must be one of 'flat', 'age',
#'                 or 'ageSex'. If 'flat', 'batchSize' must be provided. If 'age', 'maturityCurve'
#'                 and 'batchSize' must be provided. If 'ageSex', 'femaleCurve' and 'batchSize' must
#'                 be provided. Defaults to 'flat'.
#' @param mortType the value of 'type' used in the mort() call. Must be one of 'flat', 'age',
#'                 'stock', or 'ageStock'. If 'flat', 'mortRate' must be provided. If 'age',
#'                 'ageMort' must be provided. If 'stock', 'stockMort' must be provided. If
#'                 'ageStock', 'ageStockMort' must be provided. Defaults to 'flat'.
#' @param batchSize the value of 'batchSize' used in the altMate() call. Cannot be blank.
#' @param firstBreed the value of 'firstBreed' used in the altMate() call. Defaults to 0.
#' @param maxClutch the value of 'maxClutch' used in the altMate() call. Defaults to Inf. If non-Inf,
#'                  _effective_ batchSize is estimated as the mean of 1000000 draws from the
#'                  distribution of batchSize, subsetted to those <= maxAge.
#' @param osr the value of 'osr' used in the altMate() call. Female proportion is used as a
#'            multiplier on the fecundities. Defaults to c(0.5, 0.5).
#' @param maturityCurve the value of 'maturityCurve' used in the altMate() call. check_growthrates()
#'                      only uses female fecundities in its estimates, so femaleCurve is
#'                      equivalent to maturityCurve in check_growthrates(), but maturityCurve is
#'                      used when mateType is 'age'. If both mortality and maturity are specified
#'                      as vectors, they can be of different lengths. If the maturity vector is
#'                      shorter, it is 'padded' to the same length as the mortality vector by
#'                      repeating the last value in the vector.
#' @param femaleCurve the value of 'femaleCurve' used in the altMate() call. check_growthrates()
#'                    only uses female fecundities in its estimates, so femaleCurve is
#'                    equivalent to maturityCurve in check_growthrates(), but femaleCurve is used
#'                    when 'mateType' is 'ageSex'. If both mortality and maturity are specified
#'                    as vectors, they can be of different lengths. If the maturity vector is
#'                    shorter, it is 'padded' to the same length as the mortality vector by
#'                    repeating the last value in the vector.
#' @param maxAge the value of 'maxAge' used in the mort() call. Defaults to Inf.
#' @param mortRate the value of 'mortRate' used in the mort() call
#' @param ageMort the value of 'ageMort' used in the mort() call. If both mortality and maturity are
#'                specified as vectors, they can be of different lengths. If the mortality vector is
#'                shorter, it is 'padded' to the same length as the maturity vector by repeating the
#'                last value in the vector.
#' @param stockMort the value of 'stockMort' used in the mort() call
#' @param ageStockMort the value of 'ageStockMort' used in the mort() call. If both mortality and
#'                maturity are specified as vectors, they can be of different lengths. If the
#'                mortality vector is shorter, it is 'padded' to the same length as the maturity
#'                vector by repeating the last value in the vector.
#' @seealso [fishSim::PoNG()]
#' @export

check_growthrate <- function(forceY1 = NA, mateType = "flat", mortType = "flat", batchSize,
                             firstBreed = 0, maxClutch = Inf, osr = c(0.5, 0.5), maturityCurve,
                             femaleCurve, maxAge = Inf, mortRate, ageMort, stockMort, ageStockMort) {

    if(!(mateType %in% c("flat", "age", "ageSex"))) {
        stop("'mateType' must be one of 'flat', 'age', or 'ageSex'.")
    }
    
    if(!(mortType %in% c("flat", "age", "stock", "ageStock"))) {
        stop("'mortType' must be one of 'flat', 'age', 'stock', or 'ageStock'.")
    }

    if(missing(batchSize)) stop("'batchSize' must be specified.")

    ## re-calculate batchSize if maxClutch is non-Inf    
    if(batchSize != Inf) {
        batches <- rpois(1000000, lambda = batchSize)
        batchSize <- mean(batches[batches <= maxClutch])
    }

    if(mateType == "flat") {
        if(mortType == "flat") {
            mat <- matrix(data = 0, nrow = length(0:firstBreed)+1, ncol = length(0:firstBreed)+1)
            mat[1,((2+firstBreed):ncol(mat))] <- batchSize*osr[2]  ##for all mateType == "flat"
            for(i in 1:ncol(mat) ) {
                if((i+1) <= nrow(mat)) mat[i+1,i] <- 1-mortRate
            }
            mat[nrow(mat), ncol(mat)] <- 1-mortRate            
            ## build a Leslie matrix
        } else if(mortType == "age") {
            mat <- matrix(data = 0, nrow = length(ageMort)+1, ncol = length(ageMort)+1)
            mat[1,((2+firstBreed):ncol(mat))] <- batchSize*osr[2]  ##for all mateType == "flat"
            for(i in 1:ncol(mat) ) {
                if((i+1) <= nrow(mat)) mat[i+1, i] <- 1-ageMort[i]
            }
            mat[nrow(mat), ncol(mat)] <- 1-ageMort[length(ageMort)]
            ## build a Leslie matrix
        } else if(mortType == "stock") {
            mat <- matrix(data = 0, nrow = length(0:firstBreed)+1, ncol = length(0:firstBreed)+1)
            mat.l <- lapply(seq_len(length(stockMort)), function(X) mat) ## list of empty 'mat'
            for(s in 1:length(stockMort)) {
                mat.l[[s]][1,((2+firstBreed):ncol(mat.l[[s]]))] <- batchSize*osr[2]
                for( i in 1:ncol(mat.l[[s]]) ) {
                    if((i+1) <= nrow(mat.l[[s]])) mat.l[[s]][i+1,i] <- 1-stockMort[s]
                }
                mat.l[[s]][nrow(mat.l[[s]]), ncol(mat.l[[s]])] <- 1-stockMort[s]
            }
            ## build a list of Leslie matrices
        } else if(mortType == "ageStock") {
            mat <- matrix(data = 0, nrow = length(ageStockMort[,1])+1,
                          ncol = length(ageStockMort[,1])+1)
            mat.l <- lapply(seq_len(ncol(ageStockMort)), function(X) mat) ## list of empty 'mat'
            for(s in 1:ncol(ageStockMort)) {
                mat.l[[s]][1,((2+firstBreed):ncol(mat.l[[s]]))] <- batchSize*osr[2]
                for( i in 1:ncol(mat.l[[s]]) ) {
                    if((i+1) <= nrow(mat.l[[s]])) mat.l[[s]][i+1,i] <- 1-ageStockMort[i,s]
                }
                mat.l[[s]][nrow(mat.l[[s]]),ncol(mat.l[[s]])] <- 1-ageStockMort[nrow(ageStockMort),s]
            }
            ## build a list of Leslie matrices
        }
    } else if( mateType == "age") {
        if(firstBreed > 0) maturityCurve[1:(firstBreed-1)] <- 0 ## truncates maturity by firstBreed
        if(mortType == "flat") {
            mat <- matrix(data = 0, nrow = length(maturityCurve), ncol = length(maturityCurve))
            mat[1,] <- maturityCurve*batchSize*osr[2]
            for(i in 1:ncol(mat) ) {
                if((i+1) <= nrow(mat)) mat[i+1,i] <- 1-mortRate
            }
            mat[nrow(mat), ncol(mat)] <- 1-mortRate
            ## build a Leslie matrix
        } else if(mortType == "age") {
            mat <- matrix(data = 0, nrow = max(c(length(maturityCurve), length(ageMort))),
                          ncol = max(c(length(maturityCurve), length(ageMort))))
            mat[1,(1:length(maturityCurve))] <- maturityCurve*batchSize*osr[2]
            mat[1,(length(maturityCurve):
                   ncol(mat))] <- maturityCurve[length(maturityCurve)]*batchSize*osr[2]
            ## 'padding' of maturity is done here.
            for(i in 1:length(ageMort)) {
                if((i+1) <= nrow(mat)) mat[i+1,i] <- 1-ageMort[i]
            }
            for(i in length(ageMort):ncol(mat) ) {
                if((i+1) <= nrow(mat)) mat[i+1,i] <- 1-ageMort[length(ageMort)]
            }
            mat[nrow(mat), ncol(mat)] <- 1 - ageMort[length(ageMort)]
            ## build a Leslie matrix
        } else if(mortType == "stock") {
            mat <- matrix(data = 0, nrow = length(maturityCurve), ncol = length(maturityCurve))
            mat.l <- lapply(seq_len(length(stockMort)), function(X) mat) ## list of empty 'mat'
            for(s in 1:length(stockMort)) {
                mat.l[[s]][1,(1:length(maturityCurve))] <- maturityCurve*batchSize*osr[2]
                mat.l[[s]][1,(length(maturityCurve):
                           ncol(mat.l[[s]]))] <- maturityCurve[length(maturityCurve)]*batchSize*osr[2]
                for(i in 1:ncol(mat.l[[s]])) {
                    if((i+1) <= nrow(mat.l[[s]])) mat.l[[s]][i+1,i] <- 1-stockMort[s]
                }
                mat.l[[s]][nrow(mat.l[[s]]), ncol(mat.l[[s]])] <- 1 - stockMort[s]
            }
            ## build a list of Leslie matrices
        } else if(mortType == "ageStock") {
            mat <- matrix(data = 0, nrow = max(c(length(maturityCurve), nrow(ageStockMort))),
                          ncol = max(c(length(maturityCurve), nrow(ageStockMort))))
            mat.l <- lapply(seq_len(ncol(ageStockMort)), function(X) mat) ## list of empty 'mat'
            for(s in 1:ncol(ageStockMort)) {
                mat.l[[s]][1,(1:length(maturityCurve))] <- maturityCurve*batchSize*osr[2]
                mat.l[[s]][1,(length(maturityCurve):
                           ncol(mat.l[[s]]))] <- maturityCurve[length(maturityCurve)]*batchSize*osr[2]
                for(i in 1:ncol(mat.l[[s]])) {
                    if((i+1) <= nrow(mat.l[[s]])) mat.l[[s]][i+1,i] <- 1-ageStockMort[i,s]
                }
                for(i in nrow(ageStockMort):ncol(mat.l[[s]])) {
                    if((i+1) <= nrow(mat.l[[s]])) mat.l[[s]][i+1,i] <- 1-ageStockMort[nrow(ageMort),s]
                }
                mat.l[[s]][nrow(mat.l[[s]]), ncol(mat.l[[s]])] <- 1-ageStockMort[nrow(ageStockMort),s]
            }
            ## build a list of Leslie matrices
        }
    } else if( mateType == "ageSex") {
        if(mortType == "flat") {
            mat <- matrix(data = 0, nrow = length(femaleCurve), ncol = length(femaleCurve))
            mat[1,] <- femaleCurve*batchSize*osr[2]
            for(i in 1:ncol(mat) ) {
                if((i+1) <= nrow(mat)) mat[i+1,i] <- 1-mortRate
            }
            mat[nrow(mat), ncol(mat)] <- 1-mortRate
            ## build a Leslie matrix
        } else if(mortType == "age") {
            mat <- matrix(data = 0, nrow = max(c(length(femaleCurve), length(ageMort))),
                          ncol = max(c(length(femaleCurve), length(ageMort))))
            mat[1,(1:length(femaleCurve))] <- femaleCurve*batchSize*osr[2]
            mat[1,(length(femaleCurve):
                   ncol(mat))] <- femaleCurve[length(femaleCurve)]*batchSize*osr[2]
            ## 'padding' of maturity is done here.
            for(i in 1:length(ageMort)) {
                if((i+1) <= nrow(mat)) mat[i+1,i] <- 1-ageMort[i]
            }
            for(i in length(ageMort):ncol(mat) ) {
                if((i+1) <= nrow(mat)) mat[i+1,i] <- 1-ageMort[length(ageMort)]
            }
            mat[nrow(mat), ncol(mat)] <- 1 - ageMort[length(ageMort)]
            ## build a Leslie matrix
        } else if(mortType == "stock") {
            mat <- matrix(data = 0, nrow = length(femaleCurve), ncol = length(femaleCurve))
            mat.l <- lapply(seq_len(length(stockMort)), function(X) mat) ## list of empty 'mat'
            for(s in 1:length(stockMort)) {
                mat.l[[s]][1,(1:length(femaleCurve))] <- femaleCurve*batchSize*osr[2]
                mat.l[[s]][1,(length(femaleCurve):
                           ncol(mat.l[[s]]))] <- femaleCurve[length(femaleCurve)]*batchSize*osr[2]
                for(i in 1:ncol(mat.l[[s]])) {
                    if((i+1) <= nrow(mat.l[[s]])) mat.l[[s]][i+1,i] <- 1-stockMort[s]
                }
                mat.l[[s]][nrow(mat.l[[s]]), ncol(mat.l[[s]])] <- 1 - stockMort[s]
            }
            ## build a list of Leslie matrices
        } else if(mortType == "ageStock") {
            mat <- matrix(data = 0, nrow = max(c(length(femaleCurve), nrow(ageStockMort))),
                          ncol = max(c(length(femaleCurve), nrow(ageStockMort))))
            mat.l <- lapply(seq_len(ncol(ageStockMort)), function(X) mat) ## list of empty 'mat'
            for(s in 1:ncol(ageStockMort)) {
                mat.l[[s]][1,(1:length(femaleCurve))] <- femaleCurve*batchSize*osr[2]
                mat.l[[s]][1,(length(femaleCurve):
                           ncol(mat.l[[s]]))] <- femaleCurve[length(femaleCurve)]*batchSize*osr[2]
                for(i in 1:ncol(mat.l[[s]])) {
                    if((i+1) <= nrow(mat.l[[s]])) mat.l[[s]][i+1,i] <- 1-ageStockMort[i,s]
                }
                for(i in nrow(ageStockMort):ncol(mat.l[[s]])) {
                    if((i+1) <= nrow(mat.l[[s]])) mat.l[[s]][i+1,i] <- 1-ageStockMort[nrow(ageMort),s]
                }
                mat.l[[s]][nrow(mat.l[[s]]), ncol(mat.l[[s]])] <- 1-ageStockMort[nrow(ageStockMort),s]
            }
            ## build a list of Leslie matrices
        }
    }
    if(!is.na(forceY1)) {
        if(mortType %in% c("flat", "age")) mat[2,1] <- 1-forceY1
        if(mortType %in% c("stock", "ageStock")) for(i in 1:length(mat.l)) mat.l[[i]][2,1] <- 1-forceY1
    }  ## optionally forces year 1 mortality rate to a specific value.
        
    if(mortType %in% c("flat", "age")) {
        return(eigen(mat)$values[1]) ## this is a numeric growth rate.
    }
    if(mortType %in% c("stock", "ageStock")) {
        outs <- c(rep(NA, length(mat.l)))
        for(i in 1:length(outs)) outs[i] <- eigen(mat.l[[i]])$values[1]
        return(outs)
    }
}

##################################################################################################

#' PoNG(): find a Point of No Growth (by tweaking first-year mortality).
#'
#' This is essentially a wrapper for check_growthrate(), taking all the same inputs. It returns
#' a plot of projected growth-rates across all possible first-year survival rates and the numeric
#' first-year survival rate at which zero population growth is expected (via uniroot, or something).
#' The reasoning behind this utility is that often in biological systems, adult survival rates and
#' fecundities may be quite well-characterised, and population growth rates may also be well
#' characterised, but first-year survival may be nearly impossible to assess. This utility allows
#' a value of first-year survival to be chosen such that the population size does not change,
#' while leaving all well-characterised adult survival parameters unchanged. It is also useful for
#' answering questions of the form: 'how high would our first-year survival have to be, in order
#' for this population to _not_ be in decline?', which will surely be familiar in applied management
#' situations. Note that PoNG() uses some brute-force methods on the back end, so
#' it's not terribly efficient. It takes around 5 - 10 seconds per stock on a fairly-modern
#' (vintage 2018) laptop.
#' 
#' @param mateType the value of 'type' used in the altMate() call. Must be one of 'flat', 'age',
#'                 or 'ageSex'. If 'flat', 'batchSize' must be provided. If 'age', 'maturityCurve'
#'                 and 'batchSize' must be provided. If 'ageSex', 'femaleCurve' and 'batchSize' must
#'                 be provided. Defaults to 'flat'.
#' @param mortType the value of 'type' used in the mort() call. Must be one of 'flat', 'age',
#'                 'stock', or 'ageStock'. If 'flat', 'mortRate' must be provided. If 'age',
#'                 'ageMort' must be provided. If 'stock', 'stockMort' must be provided. If
#'                 'ageStock', 'ageStockMort' must be provided. Defaults to 'flat'.
#' @param batchSize the value of 'batchSize' used in the altMate() call. Cannot be blank.
#' @param firstBreed the value of 'firstBreed' used in the altMate() call. Defaults to 0.
#' @param maxClutch the value of 'maxClutch' used in the altMate() call. Defaults to Inf. If non-Inf,
#'                  _effective_ batchSize is estimated as the mean of 1000000 draws from the
#'                  distribution of batchSize, subsetted to those <= maxAge.
#' @param osr the value of 'osr' used in the altMate() call. Female proportion is used as a
#'            multiplier on the fecundities. Defaults to c(0.5, 0.5).
#' @param maturityCurve the value of 'maturityCurve' used in the altMate() call. check_growthrates()
#'                      only uses female fecundities in its estimates, so femaleCurve is
#'                      equivalent to maturityCurve in check_growthrates(), but maturityCurve is
#'                      used when mateType is 'age'. If both mortality and maturity are specified
#'                      as vectors, they can be of different lengths. If the maturity vector is
#'                      shorter, it is 'padded' to the same length as the mortality vector by
#'                      repeating the last value in the vector.
#' @param femaleCurve the value of 'femaleCurve' used in the altMate() call. check_growthrates()
#'                    only uses female fecundities in its estimates, so femaleCurve is
#'                    equivalent to maturityCurve in check_growthrates(), but femaleCurve is used
#'                    when 'mateType' is 'ageSex'. If both mortality and maturity are specified
#'                    as vectors, they can be of different lengths. If the maturity vector is
#'                    shorter, it is 'padded' to the same length as the mortality vector by
#'                    repeating the last value in the vector.
#' @param maxAge the value of 'maxAge' used in the mort() call. Defaults to Inf.
#' @param mortRate the value of 'mortRate' used in the mort() call
#' @param ageMort the value of 'ageMort' used in the mort() call. If both mortality and maturity are
#'                specified as vectors, they can be of different lengths. If the mortality vector is
#'                shorter, it is 'padded' to the same length as the maturity vector by repeating the
#'                last value in the vector.
#' @param stockMort the value of 'stockMort' used in the mort() call
#' @param ageStockMort the value of 'ageStockMort' used in the mort() call. If both mortality and
#'                maturity are specified as vectors, they can be of different lengths. If the
#'                mortality vector is shorter, it is 'padded' to the same length as the maturity
#'                vector by repeating the last value in the vector.
#' @seealso [fishSim::check_growthrate()]
#' @export
#' @examples
#' batchSize = 0.8
#' firstBreed = 1
#' mortRate = 0.2
#' PoNG(batchSize = batchSize, firstBreed = firstBreed, mortRate = mortRate)
#'
#' mortType = "stock"
#' stockMort = c(0.2, 0.3, 0.5)
#' firstBreed = 1
#' batchSize = 0.9
#' PoNG(mortType = "stock", batchSize = batchSize, firstBreed = firstBreed, stockMort = stockMort)
#'  ## note that only two of the stocks return a valid PoNG - with 0.5 mortality, stock 3 cannot
#'  ## reach null growth with any first-year survival rate between 0 and 1.

PoNG <- function(mateType = "flat", mortType = "flat", batchSize, firstBreed = 0,
                 maxClutch = Inf, osr = c(0.5, 0.5), maturityCurve, femaleCurve,
                 maxAge = Inf, mortRate, ageMort, stockMort, ageStockMort) {

    testPoints <- c(seq(0,1,0.01))
    outs <- c(rep(NA, length(testPoints)))
    if(mortType %in% c("flat", "age")) {
        for( i in 1:length(testPoints)) {
            outs[i] <- check_growthrate(mateType = mateType, mortType = mortType, batchSize=batchSize,
                                        firstBreed = firstBreed, maxClutch = maxClutch, osr = osr,
                                        maturityCurve = maturityCurve, femaleCurve = femaleCurve,
                                        maxAge = maxAge, forceY1 = testPoints[i], mortRate = mortRate,
                                        ageMort = ageMort, stockMort = stockMort,
                                        ageStockMort = ageStockMort)
        }  ## will need to make a special case for stock-structured populations.
        plot(x = testPoints, y = outs, type = "l", col = "red", xlab = "First-year mortality rate",
             ylab = "Population growth rate")
        abline(h = 1, lty = 2)
        if(!all(Re(outs) < 1) & !all(Re(outs) > 1) ) {
##            PNG <- testPoints[abs(1 - outs) == min(abs(1 - outs))]
            PNG <- uniroot(function(x) {Re(check_growthrate(forceY1 = x,
                                                         mateType = mateType, mortType = mortType,
                                            batchSize = batchSize, firstBreed = firstBreed,
                                            maxClutch = maxClutch, osr = osr,
                                            maturityCurve = maturityCurve,
                                            femaleCurve = femaleCurve, maxAge = maxAge,
                                            mortRate = mortRate, ageMort = ageMort,
                                            stockMort = stockMort, ageStockMort = ageStockMort))-1},
                           interval = c(0,1), tol = 0.001)
            abline(v = PNG$root, lty = 2)
            return(PNG) ## Point of No Growth
        }
    } else if (mortType %in% c("stock", "ageStock")) {
        if(mortType == "stock") {
            testPoints.m <- matrix(data = testPoints, nrow = length(testPoints),
                                   ncol = length(stockMort))
            outs.m <- matrix(data = outs, nrow = length(testPoints), ncol = length(stockMort))
        }
        if(mortType == "ageStock") {
            testPoints.m <- matrix(data = testPoints, nrow = length(testPoints),
                                   ncol = ncol(ageStockMort))
            outs.m <- matrix(data = data, nrow = length(testPoints), ncol = ncol(ageStockMort))
        }
        for( i in 1:length(testPoints)) {
            outs.m[i,] <- Re(check_growthrate(mateType = mateType, mortType = mortType,
                                        batchSize=batchSize,
                                        firstBreed = firstBreed, maxClutch = maxClutch, osr = osr,
                                        maturityCurve = maturityCurve, femaleCurve = femaleCurve,
                                        maxAge = maxAge, forceY1 = testPoints[i], mortRate = mortRate,
                                        ageMort = ageMort, stockMort = stockMort,
                                        ageStockMort = ageStockMort))
        } ## This is the special case for stock-structured populations.
        plot(x = testPoints.m[,1], y = outs.m[,1], type = "l", col = "black",
             xlim = c(min(testPoints.m), max(testPoints.m)),
             ylim = c(min(outs.m), max(outs.m)),
             xlab = "First-year mortality rate",
             ylab = "Population growth rate")
        for(i in 2:ncol(testPoints.m)) {
            lines(x = testPoints.m[,i], y = outs.m[,i], col = i)
        }
        PNGs <- c(rep(NA, ncol(outs.m)))
        for(i in 1:length(PNGs)) {
            if(!all(Re(outs.m[,i]) < 1) & !all(Re(outs.m[,i]) > 1) ) {
##                PNGs[i] <- testPoints.m[(abs(1 - outs.m[,i]) == min(abs(1 - outs.m[,i]))),i]
                PNGs[i] <- uniroot(function(x) {(Re(check_growthrate(forceY1 = x,
                                            mateType = mateType, mortType = mortType,
                                            batchSize = batchSize, firstBreed = firstBreed,
                                            maxClutch = maxClutch, osr = osr,
                                            maturityCurve = maturityCurve,
                                            femaleCurve = femaleCurve, maxAge = maxAge,
                                            mortRate = mortRate, ageMort = ageMort,
                                            stockMort = stockMort,
                                            ageStockMort = ageStockMort)[i])-1)},
                                   interval = c(0,1), tol = 0.001)$root
                abline(v = PNGs, lty = 2)
                abline(h = 1, lty = 2)
            }
        }
        return(PNGs) ## Points of No Growth
    }
}







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
