#' mate(): find male-female pairs within the same stock, and generate batches of offspring
#' until a recruitment quota is filled.
#'
#' Returns an individual matrix with added newborns at age 0.
#' 
#' @param indiv A matrix of individuals, e.g., generated in makeFounders() or output from move()
#' @param fecundity mean number of recruits to generate, as a proportion of the number
#'                  of animals in 'indiv'. nrow(indiv)*fecundity is the annual turnover -
#'                  i.e., the number of animals killed in 'mortality' will equal
#'                  nrow(indiv)*fecundity in order to keep the population size constant.
#' @param batchSize The mean number of offspring produced per pairing. Follows ~Poisson.
#' @param osr The sex ratio at birth (recruitment). Used to assign sexes to new offspring.
#' @param year Intended to be used in a simulation loop - this will be the iteration number, and
#'             holds the 'birthyear' value to give to new recruits.
#' @export

mate <- function(indiv = makeFounders(), fecundity = 0.2, batchSize = 0.5, osr = c(0.5,0.5),
                 year = "-1") {
    require(ids)
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
        if(nrow(fathersInStock) > 1) { ##add exemption case if nrow(fathersInStock) = 1 later
            drawFather <- fathersInStock[sample(nrow(fathersInStock), size = 1, replace = FALSE),]
            n.sprogs <- rTruncPoisson(n = 1, T = batchSize) ## age-dependent fecundity could go here.
            batch <- matrix(data = NA, nrow = n.sprogs, ncol = ncol(indiv))
            batch[,1] <- uuid(n = nrow(batch), drop_hyphens = TRUE)
            batch[,2] <- sample(c("M", "F"), size = nrow(batch), replace = TRUE, prob = osr)
            batch[,3] <- drawFather[1]
            batch[,4] <- drawMother[1]
            batch[,5] <- year
            ## [,6] stays 'NA', because no death year yet.
            batch[,7] <- drawMother[7] ## always recruits into the parents' stock.
            batch[,8] <- 0 ## it's a zero-year-old until the first 'birthday' step.

            if ((ticker + n.sprogs) <= nrow(sprog.m)) {
                sprog.m[ticker:(ticker+n.sprogs-1),] <- batch
            } else if ((ticker + n.sprogs) > nrow(sprog.m)) {
                sprog.m[ticker:nrow(sprog.m),] <- batch[1:(nrow(sprog.m)-(ticker-1)),1:8]
            }
        }
        ticker <- ticker+n.sprogs
    }
    indiv <- rbind(indiv, sprog.m)
    return(indiv)
}

