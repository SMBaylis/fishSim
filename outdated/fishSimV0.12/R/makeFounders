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
