library(fishSim)
library(tinytest)
library(ids)

## Unit tests for fishSim:

## 1: move

### deterministic movement moves everyone from the expected place to the expected place

## move 100% from 1 to 2
indiv <- makeFounders( stocks = c(1, 0, 0)) ## everyone starts in stock 1
## 'from' by row, 'to' by column, 'age' by aisle
moveMat <- matrix( data = c(0,1,0,0,0,0,0,0,0), ncol = 3, nrow = 3, byrow = TRUE)
indiv2 <- move( indiv, moveMat = moveMat)
expect_true( all( indiv2$Stock == 2))

## everyone moves up 1 stock:
indiv <- makeFounders( stocks = c(0.3, 0.3, 0.4)) ## mix of starting stocks
## move everyone up one stock (1 -> 2, 2 -> 3, 3 -> 1)
moveMat <- matrix( data = c(0,1,0,0,0,1,1,0,0), ncol = 3, nrow = 3, byrow = TRUE)
indiv2 <- move( indiv, moveMat = moveMat)
expect_true( all( indiv2$Stock[ indiv$Stock == 1] == 2)) ## 1 goes to 2
expect_true( all( indiv2$Stock[ indiv$Stock == 2] == 3)) ## 2 goes to 3
expect_true( all( indiv2$Stock[ indiv$Stock == 3] == 1)) ## 3 goes to 1

### sex-specific deterministic movement moves everyone to and from the expected places, by sex
## all males move to stock 1, all females move to stock 2
indiv <- makeFounders( stocks = c(0.3, 0.3, 0.4)) ## mix of starting stocks
moveMat_M <- matrix( data = c(1,0,0,1,0,0,1,0,0), ncol = 3, nrow = 3, byrow = TRUE)
moveMat_F <- matrix( data = c(0,1,0,0,1,0,0,1,0), ncol = 3, nrow = 3, byrow = TRUE)
indiv2 <- move( indiv, moveMat_M = moveMat_M, moveMat_F = moveMat_F)
expect_true( all( indiv2$Stock[ indiv2$Sex == "M"] == 1))
expect_true( all( indiv2$Stock[ indiv2$Sex == "F"] == 2))

### dead animals don't move
## repeat the sex-specific case because it's the most complicated, but make some animals dead first
indiv <- makeFounders( stocks = c(0.3, 0.3, 0.4)) ## mix of starting stocks
moveMat_M <- matrix( data = c(1,0,0,1,0,0,1,0,0), ncol = 3, nrow = 3, byrow = TRUE)
moveMat_F <- matrix( data = c(0,1,0,0,1,0,0,1,0), ncol = 3, nrow = 3, byrow = TRUE)
indiv2 <- mort( type = "flat", mortRate = 0.5)
indiv3 <- move( indiv2, moveMat_M = moveMat_M, moveMat_F = moveMat_F)
## not all animals have moved according to the rule...
expect_false( all( indiv2$Stock[ indiv2$Sex == "M"] == 1))
expect_false( all( indiv2$Stock[ indiv2$Sex == "F"] == 2))
## but all the _live_ ones have
expect_true( all( indiv3$Stock[ indiv3$Sex == "M" & is.na(indiv3$DeathY)] == 1))
expect_true( all( indiv3$Stock[ indiv3$Sex == "F" & is.na(indiv3$DeathY)] == 2))

### deterministic _non_-movement doesn't move anyone
indiv <- makeFounders( stocks = c(0.3, 0.3, 0.4)) ## mix of starting stocks
## move everyone up one stock (1 -> 2, 2 -> 3, 3 -> 1)
moveMat <- matrix( data = c(1,0,0,0,1,0,0,0,1), ncol = 3, nrow = 3, byrow = TRUE)
indiv2 <- move( indiv, moveMat = moveMat)
expect_true( all( indiv$Stock == indiv2$Stock)) ## everyone stays in the same place

### age-specific movement moves the right ages to the right places
#### move all three year-olds to stock 1 and 5 year-olds to stock3; leave everyone else in place
indiv <- makeFounders(pop = 10000, minAge = 1, maxAge = 5)
moveMat <- array(data = c(1,0,0,0,1,0,0,0,1, ## zero year-olds (if any) stay put
                          1,0,0,0,1,0,0,0,1, ## one year-olds stay put
                          1,0,0,0,1,0,0,0,1, ## two year-olds stay put
                          1,1,1,0,0,0,0,0,0, ## three year-olds all go to stock 1
                          1,0,0,0,1,0,0,0,1, ## four year-olds stay put
                          0,0,0,0,0,0,1,1,1), ## five year-olds all go to stock 3
                 dim = c(3,3,6))
indiv2 <- move( indiv, moveMat = moveMat)
## expect all 0, 1, 2, 4 year-olds to have the same stock as before
expect_true(
    all( indiv2$Stock[ indiv$AgeLast %in% c(0,1,2,4)] == indiv$Stock[ indiv$AgeLast %in% c(0,1,2,4)] )
) ## all the ones that shouldn't have moved haven't moved
expect_true( all( indiv2$Stock[ indiv$AgeLast == 3] == 1 )) ## all 3 year-olds have moved to stock 1
expect_true( all( indiv2$Stock[ indiv$AgeLast == 5] == 3 )) ## all 5 year-olds have moved to stock 3

### ... including zero year-olds
indiv <- makeFounders(pop = 10000, minAge = 0, maxAge = 5)
moveMat <- array(data = c(0,0,0,1,1,1,0,0,0, ## zero year-olds move to stock 2
                          1,0,0,0,1,0,0,0,1, ## one year-olds stay put
                          1,0,0,0,1,0,0,0,1, ## two year-olds stay put
                          1,1,1,0,0,0,0,0,0, ## three year-olds all go to stock 1
                          1,0,0,0,1,0,0,0,1, ## four year-olds stay put
                          0,0,0,0,0,0,1,1,1), ## five year-olds all go to stock 3
                 dim = c(3,3,6))
indiv2 <- move( indiv, moveMat = moveMat)
## expect all 0, 1, 2, 4 year-olds to have the same stock as before
expect_true(
    all( indiv2$Stock[ indiv$AgeLast %in% c(1,2,4)] == indiv$Stock[ indiv$AgeLast %in% c(1,2,4)] )
) ## all the ones that shouldn't have moved haven't moved
expect_true( all( indiv2$Stock[ indiv$AgeLast == 3] == 1 )) ## all 3 year-olds have moved to stock 1
expect_true( all( indiv2$Stock[ indiv$AgeLast == 5] == 3 )) ## all 5 year-olds have moved to stock 3
expect_true( all( indiv2$Stock[ indiv$AgeLast == 0] == 2 )) ## all 5 year-olds have moved to stock 3

## 2: mate

### only animals at or above firstBreed breed
indiv <- makeFounders()
firstBreed <- 2
indiv2 <- mate( indiv, year = 1, firstBreed = firstBreed)
parents <- indiv2[ indiv2$Me %in% c(indiv2$Mum, indiv2$Dad),]
expect_true( all( parents$AgeLast >= firstBreed))

### dead animals don't breed
indiv <- makeFounders()
indiv2 <- mort( indiv, year = 999, mortRate = 0.3, maxAge = 7)
indiv3 <- mate( indiv2, year = 999)
parents <- indiv3$Me[ indiv3$Me %in% c(indiv3$Mum, indiv3$Dad)]
deaders <- indiv3$Me[ !is.na( indiv3$DeathY)]

expect_false( any( deaders %in% parents))
expect_false( any( parents %in% deaders))


### breeders include all eligible age-classes
### This one is tricky, because breeding or non-breeding is non-deterministic
### for individuals. Have just a few age-classes and large number of offspring
### with a small batch-size, and we should be OK.

#### type = 'flat'
#### type = 'age'
#### type = 'ageSex'

### breeders include all eligible age-classes, even if age-classes are different by sex
#### type = 'flat'
#### type = 'age'
#### type = 'ageSex'

### n recruits is a reasonable approximation of [available breeders]x[fecundity]

### 'firstBreed' takes precedence over fecundityCurve, maleCurve, and femaleCurve

### fathers only have offspring with one mother if exhaustFathers = TRUE

### mothers only have offspring with one father if exhaustMothers = TRUE

## 2a: altMate

### only animals above firstBreed breed

### dead animals don't breed

### breeders include all eligible age-classes
#### type = 'flat'
#### type = 'age'
#### type = 'ageSex'

### breeders include all eligible age-classes, even if age-classes are different by sex
#### type = 'flat'
#### type = 'age'
#### type = 'ageSex'

### 'firstBreed' takes precedence over fecundityCurve, maleCurve, and femaleCurve

### fathers only have offspring with one mother if exhaustFathers = TRUE

### mothers only have offspring with one father if singlePaternity = TRUE

### mothers can have offspring with multiple fathers if singlePaternity = FALSE

## 3: mort

### if type = 'simple', living pop is brought down to exactly maxPop
#### if living pop is already less than maxPop, no animals are killed

### no animals above maxAge survive, but some _at_ maxAge survive

### ppn of animals killed is within 15% of mortRate if type = "flat"

### if type = "age", mortality exactly follows a (0, 1, 0, 1, ...) ageMort curve

### if type = "stock", mortality exactly follows a (0, 1, 0) stockMort

### if type = "ageStock", mortality exactly follows a checkerboard mort matrix

## 4: birthdays

### dead animals don't age


## check_growthrates and PoNG

### forcing Y1 to 1 gives a greater estimate than forcing Y1 to 0.5

### increasing batchSize gives a higher estimate

### shifting maturityCurve earlier gives a greater estimate

### shifting mortRate higher gives a lower estimate

### setting mortRate to 1 gives a growth rate of -Inf (or something equally sensible)

### setting mortRate to 0 gives growthRate == fecundity, or similar

### setting osr to more-female gives higher growth rates

### reducing maxClutch reduces growthrate

### check_growthrates is within 15% of real growthrate if:



































check_growthrate(mateType = "flat", mortType = "flat",
                 batchSize = 0.9, firstBreed = 2, mortRate = 0.2)

## with the current settings, we expect our population to grow by about 1.8% annually.

PoNG(mateType = "flat", mortType = "flat", batchSize = 0.9, firstBreed = 2,
     mortRate = 0.2)

## this is the 'basic scenario' from the vignette. indiv is experiencing population
## growth (10x over 60 years), but PoNG suggests it shouldn't.
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

## first, check whether the loop as-written obeys the rules of check_growthrate and PoNG
## uses altMate, mort, and birthdays.
## mateType and mortType are 'flat'
mtrace(check_growthrate)











## this is a script file for exploring the basic properties of Leslie matrices.
## Because I don't currently know how they work.

## Coding up the American Bison example

mat <- matrix(data = c(0,0,0.42,0.60,0,0,0,0.75,0.8120), nrow = 3, ncol = 3, byrow = TRUE)
              ## gives near-stable pop size.
mat <- matrix(data = c(0,0,1.1,0.90,0,0,0,0.90,0.12), nrow = 3, ncol = 3, byrow = TRUE)
              ## delta-pop is 1.066 per time-step. See eigen(mat)
pop <- c(5000,0,0)  ## initialise population here
## ordered: calves, yearlings, adults

output <- matrix(data = NA, nrow = 201, ncol = 3)
output[1, ] <- pop

for (i in 2:201) {
    output[i, 1] <- floor(output[(i-1), 3] * mat[1,3])
    output[i, 2] <- floor(output[(i-1), 1] * mat[2,1])
    output[i, 3] <- floor((output[(i-1), 2] * mat[3,2]) + (output[(i-1), 3] * mat[3,3]))
}

plot(output[,1], type = "l", col = "blue", ylim = c(min(output), max(output)),
     ylab = "N individuals", xlab = "Time") ## calves
lines(output[,2], col = "purple")
lines(output[,3], col = "red")
legend(x = 1, y = max(output), legend = c("adults", "yearlings", "calves"),
       lty = 1, col = c("red", "purple", "blue"))

## ...and, the pigeon example.

mat <- matrix(data = c(0,3.5,1.5,0.5,
                       0.4,0,0,0,
                       0,0.5,0,0,
                       0,0,0.3,0,
                       0,0,0,0), nrow = 5, ncol = 4, byrow = TRUE)
pop <- c(2500,1000,300,150)

niter = 30
output <- matrix(data = NA, nrow = niter+1, ncol = 4)
output[1,] <- pop

for(i in 2:nrow(output)) {
    output[i,1] <- floor(output[(i-1), 2]*mat[1,2]) + floor(output[(i-1), 3]*mat[1,3]) +
                   floor(output[(i-1), 4]*mat[1,4])
    output[i,2] <- floor(output[(i-1), 1] * mat[2,1])
    output[i,3] <- floor(output[(i-1), 2] * mat[3,2])
    output[i,4] <- floor(output[(i-1), 3] * mat[4,3]) ## 100% terminal mortality, so no loop.
}

plot(output[,1], type = "l", col = "blue", ylim = c(min(output), max(output)),
     ylab = "N individuals", xlab = "Time") ## calves
lines(output[,2], col = "purple")
lines(output[,3], col = "red")
