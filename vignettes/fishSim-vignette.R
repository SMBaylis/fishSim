## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(mc.cores = 2)

## -----------------------------------------------------------------------------
library(fishSim)
indiv <- makeFounders(stocks = c(1))
head(indiv)

## -----------------------------------------------------------------------------
nrow(indiv)
indiv <- mate(indiv, year = 1)
nrow(indiv) ## 200 newborns
tail(indiv) ## newborns are added to the end

## -----------------------------------------------------------------------------
indiv <- makeFounders(stocks = c(1))
nrow(indiv) ## 1000 founders
indiv <- altMate(indiv, firstBreed = 2, batchSize = 0.9, year = 1)
tail(indiv)

## -----------------------------------------------------------------------------
nrow(indiv)
indiv <- mort(indiv, year = 1, type = "flat", mortRate = 0.2)
nrow(indiv)
head(indiv, n = 15)

## -----------------------------------------------------------------------------
tail(indiv)
indiv <- birthdays(indiv)
tail(indiv)

## -----------------------------------------------------------------------------
## non-lethal sampling of 5 females
indiv <- capture(indiv, n = 5, year = 1, fatal = FALSE, sex = "F")
## lethal sampling of 5 animals, from either sex
indiv <- capture(indiv, n = 5, year = 1, fatal = TRUE)
## non-lethal sampling of 5 one-year-old females
indiv <- capture(indiv, n = 5, year = 1, fatal = FALSE, sex = "F", age = 1)

## -----------------------------------------------------------------------------
check_growthrate(mateType = "flat", mortType = "flat",
                 batchSize = 0.9, firstBreed = 2, mortRate = 0.2)

## with the current settings, we expect our population to grow by about 1.8% annually.

PoNG(mateType = "flat", mortType = "flat", batchSize = 0.9, firstBreed = 2,
     mortRate = 0.2)

## if first-year mortality was about 0.305, rather than 0.2, our population
## would have null growth. You can also read off a range of possible growth rates
## from the plot.

## -----------------------------------------------------------------------------
indiv <- makeFounders(stocks = c(1))
ageMort <- c(0.305, rep(0.2, 100)) ## age-specific mortality is 0.305 for first-years,
## 0.2 for all older age-classes.

for (y in 1:60) {
    indiv <- altMate(indiv, firstBreed = 2, batchSize = 0.9, year = y) ## y for year
    indiv <- mort(indiv, year = y, type = "age", ageMort = ageMort)    ## age-specific mort
    indiv <- birthdays(indiv)
}
tail(indiv) ## a population with 60 years of births, deaths, and birthdays
nrow(indiv[is.na(indiv$DeathY),]) ## the currently-alive population size. Note that population
                               ## growth only *averages* zero, and variability occurs!

## -----------------------------------------------------------------------------
## mark 200 individuals alive at the end of the simulation as 'captured'.
indiv <- capture(indiv, n = 200, year = 60, fatal = FALSE)

## look up each animal's ancestors and look for shared ancestors between each
## pair of sampled animals:
pairs <- findRelativesPar(indiv = indiv, sampled = TRUE, nCores = 2)

POPs <- pairs[pairs$OneTwo == 1,] ## Parent-Offspring pairs
GGPs <- pairs[pairs$OneThree == 1,] ## Grandparent-Grandoffspring pairs
HSPs <- pairs[pairs$TwoTwo == 1,] ## Half-sibling pairs
FSPs <- pairs[pairs$TwoTwo == 2,] ## Full Sibling pairs (self-comparisons
## are automatically excluded)
FCPs <- pairs[pairs$ThreeThree == 2 & pairs$TwoTwo != 1,] ## Full Cousin pairs

## look at the number of shared ancestors at each ancestral
## generation, for one of the half-sibling pairs.
lookAtPair(HSPs[1,])

relatives <- namedRelatives(pairs) ## shows the number of pairs of each relationship type
relatives

## -----------------------------------------------------------------------------
## set up founders with three stocks: two that each contain 30% of the population,
## and one that contains the remaining 40%.

indiv <- makeFounders(pop = 1000, stocks = c(0.3, 0.3, 0.4))

## set up archive - just a matrix with zero rows and eight columns

archive <- make_archive()


## -----------------------------------------------------------------------------
## Markovian movement matrix
stocks <- c(0.3, 0.3, 0.4)
admix.m <- matrix(NA, nrow = length(stocks), ncol = length(stocks))
for(i in 1:nrow(admix.m)) {
    admix.m[i,] <- stocks*stocks[i]
}
## admix.m shows movement proportional to starting population sizes.
admix.m
## let's tweak those numbers so that animals tend to stay where
## they are, and not move around so much.
admix.m <- matrix(c(0.23, 0.03, 0.04, 0.03, 0.23, 0.04, 0.04, 0.04, 0.32),
                  nrow = length(stocks), ncol = length(stocks), byrow = FALSE)
admix.m

## Age- and stock-dependent survival
ageStockMort <- matrix(c(0.47, 0.37, 0.27, rep(0.23, 97),
                         0.45, 0.35, 0.25, 0.20, rep(0.22, 96),
			 0.45, 0.3, 0.3, 0.19, rep(0.2, 96)),
                       ncol = length(stocks), nrow = 100)
head(ageStockMort)

## Sex-specific maturity curves

maleCurve <- c(0,0,0,0.1,0.5,0.8,0.85,0.9,0.95,rep(1, 91))
femaleCurve <- c(0,0,0.5,0.9,0.95,rep(1,95))
## maleCurve and femaleCurve should both be long enough that no individuals
## will outlive the curve.
head(maleCurve)
head(femaleCurve)


## -----------------------------------------------------------------------------
check_growthrate(mateType = "ageSex", mortType = "ageStock", batchSize = 1.6,
                 femaleCurve = femaleCurve,
                 ageStockMort = ageStockMort)
## Not bad. Two of the three populations are increasing. The third will probably be kept viable
## by immigration from the other two. Note that I fiddled 'batchSize' (which is the mean number
## of offspring per mature female per breeding attempt) a bit.

## -----------------------------------------------------------------------------
for (k in 1:100) {
    ## very rarely, switch some males to females
    indiv <- sexSwitch(indiv = indiv, direction = "MF", prob = 1e-04)
    ## move animals according to the markovian matrix we set up before
    indiv <- move(indiv = indiv, moveMat = admix.m)
    ## mate animals using the age-specific, sex-specific curves we set up before
    indiv <- altMate(indiv = indiv, batchSize = 1.6, type = "ageSex", maleCurve = maleCurve,
                     femaleCurve = femaleCurve, year = k)
    ## kill animals on the basis of their age and stock, as set up before
    indiv <- mort(indiv = indiv, year = k, type = "ageStock", ageStockMort = ageStockMort)
    if(k %% 10 == 0) {
        archive <- archive_dead(indiv = indiv, archive = archive)
	    indiv <- remove_dead(indiv = indiv)
	    cat( sprintf( '\rIteration %i complete', k))
	## occasionally sample from the population and do some clean-up. And report on progress.
    }
    if(k %in% c(94:100)) {
	indiv <- capture(indiv, n = 40, fatal = TRUE, year = k)
    }
    indiv <- birthdays(indiv = indiv)
}
archive <- rbind(archive, indiv)
indiv <- archive  ## merge 'indiv' and 'archive', since they were only separated for speed.


## -----------------------------------------------------------------------------
## look up each animal's ancestors and look for shared ancestors between each
## pair of sampled animals:
pairs <- findRelativesPar(indiv = indiv, sampled = TRUE, nCores = 2)

POPs <- pairs[pairs$OneTwo == 1,] ## Parent-Offspring pairs
GGPs <- pairs[pairs$OneThree == 1,] ## Grandparent-Grandoffspring pairs
HSPs <- pairs[pairs$TwoTwo == 1,] ## Half-sibling pairs
FSPs <- pairs[pairs$TwoTwo == 2,] ## Full Sibling pairs (self-comparisons
## are automatically excluded)
FCPs <- pairs[pairs$ThreeThree == 2 & pairs$TwoTwo != 1,] ## Full Cousin pairs

## look at the number of shared ancestors at each ancestral
## generation, for one of the half-sibling pairs.
lookAtPair(HSPs[1,])

relatives <- namedRelatives(pairs) ## shows the number of pairs of each relationship type
relatives

## ----eval = FALSE-------------------------------------------------------------
#  indiv <- makeFounders(stocks = c(1))
#  ageMort <- c(0.305, rep(0.2, 100)) ## age-specific mortality is 0.305 for first-years,
#  ## 0.2 for all older age-classes.
#  
#  for (y in 1:60) {
#      indiv <- altMate(indiv, firstBreed = 2, batchSize = 0.9, year = y)
#      ## 'mate' happens first, affecting the whole population
#  
#      indiv <- capture(indiv, n = 5, year = y, fatal = FALSE)
#      indiv_unmarked <- indiv[is.na(indiv$SampY),] ## all the unmarked animals
#      indiv_marked <- indiv[!is.na(indiv$SampY),]  ## all the marked animals
#  
#      indiv_unmarked <- mort(indiv_unmarked, year = y, type = "age", ageMort = ageMort)
#          ## normal mort for unmarked animals
#      indiv_marked <- mort(indiv_marked, year = y, type = "age", ageMort = ageMort*2)
#          ## double mort for marked animals
#      indiv <- rbind(indiv_unmarked, indiv_marked)
#          ## stick them back together and proceed as normal
#      rm(indiv_marked, indiv_unmarked) ## clean up
#  
#      indiv <- birthdays(indiv)
#  }

