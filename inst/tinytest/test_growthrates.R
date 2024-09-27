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
indiv <- makeFounders( stocks = c(1), maxAge = 5)
indiv2 <- mate( indiv, year = 999, type = "flat", firstBreed = 2)
## should have parents of all ages from 2 -> 5
parents <- indiv2[ indiv2$Me %in% c(indiv2$Mum, indiv2$Dad),]
expect_true( all ( c(2,3,4,5) %in% parents$AgeLast))
expect_false( any ( c(0,1) %in% parents$AgeLast))

## check that again with zero year-olds in indiv, in case of OBOEs
indiv <- makeFounders( stocks = c(1), maxAge = 5, minAge = 0)
indiv2 <- mate( indiv, year = 999, type = "flat", firstBreed = 2)
## should have parents of all ages from 2 -> 5
parents <- indiv2[ indiv2$Me %in% c(indiv2$Mum, indiv2$Dad),]
expect_true( all ( c(2,3,4,5) %in% parents$AgeLast))
expect_false( any ( c(0,1) %in% parents$AgeLast))

#### type = 'age'
indiv <- makeFounders( stocks = c(1), maxAge = 5)
## fecundityCurve for ages from 0:5
fecundityCurve <- c(0,0,1,1,1,1)
indiv2 <- mate( indiv, year = 999, type = "age", firstBreed = 2, fecundityCurve = fecundityCurve)
## should have parents of all ages from 2 -> 5
parents <- indiv2[ indiv2$Me %in% c(indiv2$Mum, indiv2$Dad),]
expect_true( all ( c(2,3,4,5) %in% parents$AgeLast))
expect_false( any ( c(0,1) %in% parents$AgeLast))

## check that again with zero year-olds in indiv, in case of OBOEs
indiv <- makeFounders( stocks = c(1), maxAge = 5, minAge = 0)
## fecundityCurve for ages from 0:5
fecundityCurve <- c(0,0,1,1,1,1)
indiv2 <- mate( indiv, year = 999, type = "age", firstBreed = 2, fecundityCurve = fecundityCurve)
## should have parents of all ages from 2 -> 5
parents <- indiv2[ indiv2$Me %in% c(indiv2$Mum, indiv2$Dad),]
expect_true( all ( c(2,3,4,5) %in% parents$AgeLast))
expect_false( any ( c(0,1) %in% parents$AgeLast))

#### type = 'ageSex'
indiv <- makeFounders( stocks = c(1), maxAge = 5)
## fecundityCurve for ages from 0:5
maleCurve <- c(0,0,1,1,1,1) ## males start breeding age 2
femaleCurve <- c(0,0,0,1,1,1) ## females start breeding age 3
indiv2 <- mate( indiv, year = 999, type = "ageSex", maleCurve = maleCurve, femaleCurve = femaleCurve)
## should have male parents of all ages from 2 -> 5
dads <- indiv2[ indiv2$Me %in% c(indiv2$Dad),]
## should have female parents of all ages from 3 -> 5
mums <- indiv2[ indiv2$Me %in% c(indiv2$Mum),]

## check ages for dads
expect_true( all ( c(2,3,4,5) %in% dads$AgeLast))
expect_false( any ( c(0,1) %in% dads$AgeLast))
## check ages for mums
expect_true( all ( c(3,4,5) %in% mums$AgeLast))
expect_false( any ( c(0,1,2) %in% mums$AgeLast))

### n recruits is [population size]x[fecundity]

indiv <- makeFounders( stocks = c(1))
fec <- 0.4
indiv2 <- mate( indiv, fecundity = fec, year = 999, type = "flat", firstBreed = 2)
recruits <- indiv2[ indiv2$BirthY == 999,]
expectedN <- floor(nrow(indiv) * fec)
minExp <- floor( expectedN * 0.95) ## minus 5%
maxExp <- ceiling( expectedN * 1.05) ## plus 5%
expect_true( nrow(recruits) > minExp)
expect_true( nrow(recruits) < maxExp)

### 'firstBreed' takes precedence over fecundityCurve, maleCurve, and femaleCurve
#### type = 'flat'
indiv <- makeFounders( stocks = c(1), maxAge = 5)
indiv2 <- mate( indiv, year = 999, type = "flat", firstBreed = 4)
## should have parents only from ages 4 -> 5
parents <- indiv2[ indiv2$Me %in% c(indiv2$Mum, indiv2$Dad),]
expect_true( all ( c(4,5) %in% parents$AgeLast))
expect_false( any ( c(0,1,2,3) %in% parents$AgeLast))

#### type = 'age'
indiv <- makeFounders( stocks = c(1), maxAge = 5)
## fecundityCurve for ages from 0:5
fecundityCurve <- c(0,0,1,1,1,1)
indiv2 <- mate( indiv, year = 999, type = "age", firstBreed = 4, fecundityCurve = fecundityCurve)
## should have parents only from ages 4 -> 5
parents <- indiv2[ indiv2$Me %in% c(indiv2$Mum, indiv2$Dad),]
expect_true( all ( c(4,5) %in% parents$AgeLast))
expect_false( any ( c(0,1,2,3) %in% parents$AgeLast))

#### type = 'ageSex'
indiv <- makeFounders( stocks = c(1), maxAge = 5)
## fecundityCurve for ages from 0:5
maleCurve <- c(0,0,1,1,1,1) ## males start breeding age 2
femaleCurve <- c(0,0,0,1,1,1) ## females start breeding age 3
indiv2 <- mate( indiv, year = 999, type = "ageSex", firstBreed = 4,
               maleCurve = maleCurve, femaleCurve = femaleCurve)
## should have dads only from ages 4 -> 5
dads <- indiv2[ indiv2$Me %in% c(indiv2$Dad),]
## should have mums only from ages 4 -> 5
mums <- indiv2[ indiv2$Me %in% c(indiv2$Mum),]

## check ages for dads
expect_true( all ( c(4,5) %in% dads$AgeLast))
expect_false( any ( c(0,1,2,3) %in% dads$AgeLast))
## check ages for mums
expect_true( all ( c(4,5) %in% mums$AgeLast))
expect_false( any ( c(0,1,2,3) %in% mums$AgeLast))

### fathers only have offspring with one mother if exhaustFathers = TRUE
indiv <- makeFounders(stocks = c(1))
indiv2 <- mate( indiv, year = 999, fecundity = 0.3, batchSize = 2, exhaustFathers = TRUE)
recruits <- indiv2[ indiv2$BirthY == 999,]
dads <- unique(recruits$Dad)
mums <- unique(recruits$Mum)
couples <- data.frame( Dad = recruits$Dad, Mum = recruits$Mum,
                      pairID = paste( recruits$Dad, recruits$Mum),
                      DadsNMates = NA, MumsNMates = NA)
for( r in 1:nrow( couples)) {
    thisDad <- couples$Dad[r]
    thisMum <- couples$Mum[r]
    allThisDad <- couples[ couples$Dad == thisDad,]
    allThisMum <- couples[ couples$Mum == thisMum,]
    couples$DadsNMates[r] <- length( unique( allThisDad$Mum))
    couples$MumsNMates[r] <- length( unique( allThisMum$Dad))
}
expect_true( all( couples$DadsNMates == 1)) ## passes now (Sept 24, 2024)

### mothers only have offspring with one father if exhaustMothers = TRUE
indiv <- makeFounders(stocks = c(1))
indiv2 <- mate( indiv, year = 999, fecundity = 0.3, batchSize = 2, exhaustMothers = TRUE)
recruits <- indiv2[ indiv2$BirthY == 999,]
dads <- unique(recruits$Dad)
mums <- unique(recruits$Mum)
couples <- data.frame( Dad = recruits$Dad, Mum = recruits$Mum,
                      pairID = paste( recruits$Dad, recruits$Mum),
                      DadsNMates = NA, MumsNMates = NA)
for( r in 1:nrow( couples)) {
    thisDad <- couples$Dad[r]
    thisMum <- couples$Mum[r]
    allThisDad <- couples[ couples$Dad == thisDad,]
    allThisMum <- couples[ couples$Mum == thisMum,]
    couples$DadsNMates[r] <- length( unique( allThisDad$Mum))
    couples$MumsNMates[r] <- length( unique( allThisMum$Dad))
}
expect_true( all( couples$MumsNMates == 1)) ## passes now (Sept 24, 2024)




## 2a: altMate
### only animals above firstBreed breed
indiv <- makeFounders()
firstBreed <- 2
indiv2 <- altMate( indiv, year = 1, firstBreed = firstBreed)
parents <- indiv2[ indiv2$Me %in% c(indiv2$Mum, indiv2$Dad),]
expect_true( all( parents$AgeLast >= firstBreed))

### dead animals don't breed
indiv <- makeFounders()
indiv2 <- mort( indiv, year = 999, mortRate = 0.3, maxAge = 7)
indiv3 <- altMate( indiv2, year = 999)
parents <- indiv3$Me[ indiv3$Me %in% c(indiv3$Mum, indiv3$Dad)]
deaders <- indiv3$Me[ !is.na( indiv3$DeathY)]

expect_false( any( deaders %in% parents))
expect_false( any( parents %in% deaders))

### breeders include all eligible age-classes
#### type = 'flat'
indiv <- makeFounders( stocks = c(1), maxAge = 5)
indiv2 <- altMate( indiv, year = 999, type = "flat", firstBreed = 2)
## should have parents of all ages from 2 -> 5
parents <- indiv2[ indiv2$Me %in% c(indiv2$Mum, indiv2$Dad),]
expect_true( all ( c(2,3,4,5) %in% parents$AgeLast))
expect_false( any ( c(0,1) %in% parents$AgeLast))

## check that again with zero year-olds in indiv, in case of OBOEs
indiv <- makeFounders( stocks = c(1), maxAge = 5, minAge = 0)
indiv2 <- altMate( indiv, year = 999, type = "flat", firstBreed = 2)
## should have parents of all ages from 2 -> 5
parents <- indiv2[ indiv2$Me %in% c(indiv2$Mum, indiv2$Dad),]
expect_true( all ( c(2,3,4,5) %in% parents$AgeLast))
expect_false( any ( c(0,1) %in% parents$AgeLast))

#### type = 'age'
indiv <- makeFounders( stocks = c(1), maxAge = 5)
## maturityCurve for ages from 0:5
maturityCurve <- c(0,0,1,1,1,1)
indiv2 <- altMate( indiv, year = 999, type = "age", maturityCurve = maturityCurve)
## should have parents of all ages from 2 -> 5
parents <- indiv2[ indiv2$Me %in% c(indiv2$Mum, indiv2$Dad),]
expect_true( all ( c(2,3,4,5) %in% parents$AgeLast))
expect_false( any ( c(0,1) %in% parents$AgeLast))

## check that again with zero year-olds in indiv, in case of OBOEs
indiv <- makeFounders( stocks = c(1), maxAge = 5, minAge = 0)
## maturityCurve for ages from 0:5
maturityCurve <- c(0,0,1,1,1,1)
indiv2 <- altMate( indiv, year = 999, type = "age", maturityCurve = maturityCurve)
## should have parents of all ages from 2 -> 5
parents <- indiv2[ indiv2$Me %in% c(indiv2$Mum, indiv2$Dad),]
expect_true( all ( c(2,3,4,5) %in% parents$AgeLast))
expect_false( any ( c(0,1) %in% parents$AgeLast))

#### type = 'ageSex'
indiv <- makeFounders( stocks = c(1), maxAge = 5)
## maleCurve and femaleCurve for ages from 0:5
maleCurve <- c(0,0,1,1,1,1) ## males start breeding age 2
femaleCurve <- c(0,0,0,1,1,1) ## females start breeding age 3
indiv2 <- mate( indiv, year = 999, type = "ageSex", maleCurve = maleCurve, femaleCurve = femaleCurve)
## should have male parents of all ages from 2 -> 5
dads <- indiv2[ indiv2$Me %in% c(indiv2$Dad),]
## should have female parents of all ages from 3 -> 5
mums <- indiv2[ indiv2$Me %in% c(indiv2$Mum),]

## check ages for dads
expect_true( all ( c(2,3,4,5) %in% dads$AgeLast))
expect_false( any ( c(0,1) %in% dads$AgeLast))
## check ages for mums
expect_true( all ( c(3,4,5) %in% mums$AgeLast))
expect_false( any ( c(0,1,2) %in% mums$AgeLast))

### 'firstBreed' takes precedence over fecundityCurve, maleCurve, and femaleCurve
#### type = 'flat'
indiv <- makeFounders( stocks = c(1), maxAge = 5)
indiv2 <- altMate( indiv, year = 999, type = "flat", firstBreed = 4)
## should have parents only from ages 4 -> 5
parents <- indiv2[ indiv2$Me %in% c(indiv2$Mum, indiv2$Dad),]
expect_true( all ( c(4,5) %in% parents$AgeLast))
expect_false( any ( c(0,1,2,3) %in% parents$AgeLast))

#### type = 'age'
indiv <- makeFounders( stocks = c(1), maxAge = 5)
## fecundityCurve for ages from 0:5
maturityCurve <- c(0,0,1,1,1,1)
indiv2 <- altMate( indiv, year = 999, type = "age", firstBreed = 4, maturityCurve = maturityCurve)
## should have parents only from ages 4 -> 5
parents <- indiv2[ indiv2$Me %in% c(indiv2$Mum, indiv2$Dad),]
expect_true( all ( c(4,5) %in% parents$AgeLast))
expect_false( any ( c(0,1,2,3) %in% parents$AgeLast))

#### type = 'ageSex'
indiv <- makeFounders( stocks = c(1), maxAge = 5)
## fecundityCurve for ages from 0:5
maleCurve <- c(0,0,1,1,1,1) ## males start breeding age 2
femaleCurve <- c(0,0,0,1,1,1) ## females start breeding age 3
indiv2 <- altMate( indiv, year = 999, type = "ageSex", firstBreed = 4,
               maleCurve = maleCurve, femaleCurve = femaleCurve)
## should have dads only from ages 4 -> 5
dads <- indiv2[ indiv2$Me %in% c(indiv2$Dad),]
## should have mums only from ages 4 -> 5
mums <- indiv2[ indiv2$Me %in% c(indiv2$Mum),]

## check ages for dads
expect_true( all ( c(4,5) %in% dads$AgeLast))
expect_false( any ( c(0,1,2,3) %in% dads$AgeLast))
## check ages for mums
expect_true( all ( c(4,5) %in% mums$AgeLast))
expect_false( any ( c(0,1,2,3) %in% mums$AgeLast))

### fathers only have offspring with one mother if exhaustFathers = TRUE
indiv <- makeFounders(stocks = c(1))
indiv2 <- altMate( indiv, year = 999, batchSize = 2, exhaustFathers = TRUE)
recruits <- indiv2[ indiv2$BirthY == 999,]
dads <- unique(recruits$Dad)
mums <- unique(recruits$Mum)
couples <- data.frame( Dad = recruits$Dad, Mum = recruits$Mum,
                      pairID = paste( recruits$Dad, recruits$Mum),
                      DadsNMates = NA, MumsNMates = NA)
for( r in 1:nrow( couples)) {
    thisDad <- couples$Dad[r]
    thisMum <- couples$Mum[r]
    allThisDad <- couples[ couples$Dad == thisDad,]
    allThisMum <- couples[ couples$Mum == thisMum,]
    couples$DadsNMates[r] <- length( unique( allThisDad$Mum))
    couples$MumsNMates[r] <- length( unique( allThisMum$Dad))
}
expect_true( all( couples$DadsNMates == 1))

### mothers only have offspring with one father if singlePaternity = TRUE
indiv <- makeFounders(stocks = c(1))
indiv2 <- altMate( indiv, year = 999, batchSize = 2, singlePaternity = TRUE)
recruits <- indiv2[ indiv2$BirthY == 999,]
dads <- unique(recruits$Dad)
mums <- unique(recruits$Mum)
couples <- data.frame( Dad = recruits$Dad, Mum = recruits$Mum,
                      pairID = paste( recruits$Dad, recruits$Mum),
                      DadsNMates = NA, MumsNMates = NA)
for( r in 1:nrow( couples)) {
    thisDad <- couples$Dad[r]
    thisMum <- couples$Mum[r]
    allThisDad <- couples[ couples$Dad == thisDad,]
    allThisMum <- couples[ couples$Mum == thisMum,]
    couples$DadsNMates[r] <- length( unique( allThisDad$Mum))
    couples$MumsNMates[r] <- length( unique( allThisMum$Dad))
}
expect_true( all( couples$MumsNMates == 1))

### mothers can have offspring with multiple fathers if singlePaternity = FALSE
indiv <- makeFounders(stocks = c(1))
indiv2 <- altMate( indiv, year = 999, batchSize = 2, singlePaternity = FALSE)
recruits <- indiv2[ indiv2$BirthY == 999,]
dads <- unique(recruits$Dad)
mums <- unique(recruits$Mum)
couples <- data.frame( Dad = recruits$Dad, Mum = recruits$Mum,
                      pairID = paste( recruits$Dad, recruits$Mum),
                      DadsNMates = NA, MumsNMates = NA)
for( r in 1:nrow( couples)) {
    thisDad <- couples$Dad[r]
    thisMum <- couples$Mum[r]
    allThisDad <- couples[ couples$Dad == thisDad,]
    allThisMum <- couples[ couples$Mum == thisMum,]
    couples$DadsNMates[r] <- length( unique( allThisDad$Mum))
    couples$MumsNMates[r] <- length( unique( allThisMum$Dad))
}
expect_false( all( couples$MumsNMates == 1))




## 3: mort
### if type = 'simple', living pop is brought down to exactly maxPop
indiv <- makeFounders()
indiv2 <- mort(indiv, type = 'simple', maxPop = 789)
expect_equal( sum( is.na(indiv2$DeathY)), 789)

#### if living pop is equal to maxPop, no animals are killed
indiv <- makeFounders(pop = 788)
indiv2 <- mort(indiv, type = 'simple', maxPop = 788)
expect_equal( sum( !is.na(indiv2$DeathY)), 0)

#### if living pop is less than maxPop, no animals are killed
indiv <- makeFounders(pop = 788)
indiv2 <- mort(indiv = indiv, type = 'simple', maxPop = 1000)
expect_equal( sum( !is.na(indiv2$DeathY)), 0)

### no animals above maxAge survive, but some _at_ maxAge survive
indiv <- makeFounders(maxAge = 8)
maxAge <- 6
indiv2 <- mort(indiv, maxAge = maxAge)
survivors <- indiv2[ is.na(indiv2$DeathY),]
expect_equal( max(survivors$AgeLast), maxAge)

#### ... and that's even true if there are some zero year-olds in the mix
indiv <- makeFounders(minAge = 0, maxAge = 8)
maxAge <- 6
indiv2 <- mort(indiv, maxAge = maxAge)
survivors <- indiv2[ is.na(indiv2$DeathY),]
expect_equal( max(survivors$AgeLast), maxAge)

### ppn of animals killed is within a few % of mortRate if type = "flat"
pop <- 10000
mortRate <- 0.2
indiv <- makeFounders(pop = pop)
indiv2 <- mort(indiv, type = "flat", mortRate = mortRate)
E_toll <- pop*mortRate
O_toll <- sum( !is.na(indiv2$DeathY))
expect_true( (E_toll / O_toll) < 1.03)
expect_true( (E_toll / O_toll) > 0.97)

### if type = "age", mortality exactly follows a (1, 0, 1, ...) ageMort curve
indiv <- makeFounders(maxAge = 6)
ageMort <- c(1,0,1,0,1,0,1) ## one value for each age from 0 to 6
indiv2 <- mort( indiv, type = "age", ageMort = ageMort)
deaders <- indiv2[ !is.na(indiv2$DeathY),]
survivors <- indiv2[ is.na(indiv2$DeathY),]
expect_true( all( deaders$AgeLast %in% c(0, 2, 4, 6)))
expect_true( all( survivors$AgeLast %in% c(1, 3, 5)))

### ...and that's still the case if there are 0 year-olds in the mix
indiv <- makeFounders(maxAge = 6, minAge = 0)
ageMort <- c(1,0,1,0,1,0,1) ## one value for each age from 0 to 6
indiv2 <- mort( indiv, type = "age", ageMort = ageMort)
deaders <- indiv2[ !is.na(indiv2$DeathY),]
survivors <- indiv2[ is.na(indiv2$DeathY),]
expect_true( all( deaders$AgeLast %in% c(0, 2, 4, 6)))
expect_true( all( survivors$AgeLast %in% c(1, 3, 5)))

### if type = "stock", mortality exactly follows a (1, 0, 1) stockMort
indiv <- makeFounders()
stockMort <- c(1,0,1) ## one value for each stock from 1 to 3
indiv2 <- mort( indiv, type = "stock", stockMort = stockMort)
deaders <- indiv2[ !is.na(indiv2$DeathY),]
survivors <- indiv2[ is.na(indiv2$DeathY),]
expect_true( all( deaders$Stock %in% c(1,3)))
expect_true( all( survivors$Stock %in% c(2)))

### if type = "ageStock", mortality exactly follows a checkerboard mort matrix
indiv <- makeFounders(maxAge = 6, minAge = 0)
ageMort1 <- c(1,0,1,0,1,0,1) ## one value for each age from 0 to 6
ageMort2 <- c(0,1,0,1,0,1,0) ## the opposite pattern
ageStockMort <- matrix( c(ageMort1, ageMort2, ageMort1), nrow = length(ageMort1), ncol = 3)
indiv2 <- mort( indiv, type = "ageStock", ageStockMort = ageStockMort)

aliveS1 <- indiv2[ is.na(indiv2$DeathY) & indiv2$Stock == 1,]
aliveS2 <- indiv2[ is.na(indiv2$DeathY) & indiv2$Stock == 2,]
aliveS3 <- indiv2[ is.na(indiv2$DeathY) & indiv2$Stock == 3,]

deadS1 <- indiv2[ (!is.na(indiv2$DeathY)) & indiv2$Stock == 1,]
deadS2 <- indiv2[ (!is.na(indiv2$DeathY)) & indiv2$Stock == 2,]
deadS3 <- indiv2[ (!is.na(indiv2$DeathY)) & indiv2$Stock == 3,]

expect_true( all( deadS1$AgeLast %in% c(0, 2, 4, 6)))
expect_true( all( deadS2$AgeLast %in% c(1, 3, 5)))
expect_true( all( deadS3$AgeLast %in% c(0, 2, 4, 6)))

expect_true( all( aliveS1$AgeLast %in% c(1, 3, 5)))
expect_true( all( aliveS2$AgeLast %in% c(0, 2, 4, 6)))
expect_true( all( aliveS3$AgeLast %in% c(1, 3, 5)))

## 4: birthdays

### dead animals don't age
indiv <- makeFounders( maxAge = 6, minAge = 0)
indiv2 <- mort(indiv, year = 999, type = "flat", mortRate = 0.4)
indiv3 <- birthdays( indiv2)
isDead <- ! is.na(indiv2$DeathY)

expect_true( all( indiv3$AgeLast[ isDead] == indiv2$AgeLast[ isDead]))

### all alive animals _do_ age
expect_true( all( indiv3$AgeLast[ !isDead] == (indiv2$AgeLast[ !isDead] +1) ))



## 5: check_growthrates

### forcing Y1 to 1 gives a lesser estimate than forcing Y1 to 0.5
unforced <- check_growthrate( batchSize = 2, mortRate = 0.5, forceY1 = 0.5)
forced <- check_growthrate( batchSize = 2, mortRate = 0.5, forceY1 = 1)
expect_true( forced < unforced)

### increasing batchSize gives a higher estimate
smallbatch <- check_growthrate( batchSize = 1, mortRate = 0.5)
largebatch <- check_growthrate( batchSize = 2, mortRate = 0.5)
expect_true( largebatch > smallbatch)

### shifting maturityCurve earlier gives a greater estimate
bloomLate <- c(0,0,0,0,.5,1,1) ## 50% maturity at age 4; full maturity at 5+
bloomEarly <- c(0,0,.5,1,1,1,1) ## 50% maturity at age 2; full maturity at 3+

lateMat <- check_growthrate( batchSize = 1, mortRate = 0.2, mateType = "age",
                            maturityCurve = bloomLate)
earlyMat <- check_growthrate( batchSize = 1, mortRate = 0.2, mateType = "age",
                             maturityCurve = bloomEarly)
expect_true( earlyMat > lateMat)

### shifting mortRate higher gives a lower estimate
lomort <- check_growthrate( batchSize = 2, mortRate = 0.1)
himort <- check_growthrate( batchSize = 2, mortRate = 0.4)
expect_true( lomort > himort)

### setting mortRate to 1 gives a growth rate of 0
allmort <- check_growthrate( batchSize = 2, mortRate = 1)
expect_true( allmort == 0)

### setting mortRate to 0 gives growthRate == fecundity
nomort <- check_growthrate( batchSize = 2, mortRate = 0)
nomort ## is 1.618 (AKA the bloody Golden Ratio).
## turns out that my instinct was true for ?continuous-time
## (i.e., zero year-olds are included in the breeders)
allBreedNoMort <- check_growthrate( batchSize = 2, mortRate = 0, firstBreed = 0)
expect_true( round( allBreedNoMort, digits = 2) == 2)

### setting osr to more-female gives higher growth rates
manyfem <- check_growthrate( batchSize = 2, mortRate = 0.3, osr = c(0.2, 0.8))
fewfem <- check_growthrate( batchSize = 2, mortRate = 0.3, osr = c(0.8, 0.2))
expect_true( manyfem > fewfem)

### reducing maxClutch reduces growthrate
smallClutch <- check_growthrate( batchSize = 2, mortRate = 0.3, maxClutch = 3)
largeClutch <- check_growthrate( batchSize = 2, mortRate = 0.3, maxClutch = 6)
expect_true( largeClutch > smallClutch)

## 6: PoNG

### setting y1 survival to the value suggested by PoNG really does
### result in near-zero population growth

### ... with flat breeding and...
### ... flat mortality
batchSize = 2
mortRate = 0.3
ng <- PoNG( mateType = "flat", mortType = "flat", batchSize = batchSize,
           mortRate = mortRate)
## quick loop to test
indiv <- makeFounders(stocks = c(1))
for( y in 1:10) {
    indiv <- altMate( indiv, batchSize = batchSize, type = "flat")
    indiv <- mort( indiv, year = y, type = "age", ageMort = c(ng, rep(mortRate, 99)))
    indiv <- birthdays(indiv)
}
outsize <- nrow(indiv[ is.na(indiv$DeathY),] )
outsize
expect_true( (800 < outsize) & 1200 > outsize)

### ... age-specific mortality
batchSize = 2
ageMort = c(0.5, rep(0.3, 5), rep(0.4, 5), rep(0.5, 5)) ## senescence and immaturity
ng <- PoNG( mateType = "flat", mortType = "age", batchSize = batchSize,
           ageMort = ageMort)
## quick loop to test
indiv <- makeFounders(stocks = c(1))
for( y in 1:10) {
    indiv <- altMate( indiv, batchSize = batchSize, type = "flat")
    indiv <- mort( indiv, year = y, type = "age", ageMort = c(ng, ageMort[-1]))
    indiv <- birthdays(indiv)
}
outsize <- nrow(indiv[ is.na(indiv$DeathY),] )
expect_true( (800 < outsize) & 1200 > outsize)

### ... stock-specific mortality
batchSize = 2
stockMort = c(0.3, 0.4, 0.2)
ng <- PoNG( mateType = "flat", mortType = "stock", batchSize = batchSize,
           stockMort = stockMort)
ng <- ng ## no $root here, because PoNG only returns the values
## quick loop to test
indiv <- makeFounders(stocks = c(1))
ageStockMort <- matrix(data = c(ng, rep(stockMort,30)), byrow = TRUE, ncol = length(ng))
for( y in 1:10) {
    indiv <- altMate( indiv, batchSize = batchSize, type = "flat")
    indiv <- mort( indiv, year = y, type = "ageStock", ageStockMort = ageStockMort)
    indiv <- birthdays(indiv)
}
outsize <- nrow(indiv[ is.na(indiv$DeathY),] )
expect_true( (800 < outsize) & 1200 > outsize)

### ... age/stock specific mortality
batchSize = 2
ageStockMort <- matrix(data = c(0.3,0.4,0.2,0.35,0.4,0.3,0.38,0.4,0.35,0.4,0.4,0.4, rep(0.4,78)),
                       byrow = TRUE, ncol = 3)
ng <- PoNG( mateType = "flat", mortType = "ageStock", batchSize = batchSize,
           ageStockMort = ageStockMort)
ng <- ng ## no $root here, because PoNG only returns the values
## quick loop to test
indiv <- makeFounders(stocks = c(0.3, 0.3, 0.4))
ageStockMort <- as.matrix(rbind(ng, ageStockMort[-1,]))
for( y in 1:10) {
    indiv <- altMate( indiv, batchSize = batchSize, type = "flat")
    indiv <- mort( indiv, year = y, type = "ageStock", ageStockMort = ageStockMort)
    indiv <- birthdays(indiv)
}
outsize <- nrow(indiv[ is.na(indiv$DeathY),] )
expect_true( (800 < outsize) & 1200 > outsize)

### ... with age-specific breeding and...
### ... flat mortality
batchSize = 2
mortRate = 0.3
maturityCurve = c(0,0,0.25,0.75, rep(1, 30))

ng <- PoNG( mateType = "age", mortType = "flat", batchSize = batchSize,
           mortRate = mortRate, maturityCurve = maturityCurve)
## quick loop to test
indiv <- makeFounders(stocks = c(1))
for( y in 1:10) {
    indiv <- altMate( indiv, batchSize = batchSize, type = "age",
                     maturityCurve = c(maturityCurve))
    indiv <- mort( indiv, year = y, type = "age", ageMort = c(ng, rep(mortRate, 30)))
    indiv <- birthdays(indiv)
}
outsize <- nrow(indiv[ is.na(indiv$DeathY),] )
expect_true( (800 < outsize) & 1200 > outsize)

### ... with age- and sex-specific breeding and...
### ... flat mortality
batchSize = 2
mortRate = 0.3
maleCurve = c(0,0,0.25,0.75, rep(1, 30)) ## for males
femaleCurve = c(0,0,0.5,0.8,rep(1, 30)) ## for females

ng <- PoNG( mateType = "ageSex", mortType = "flat", batchSize = batchSize,
           mortRate = mortRate, maturityCurve = maturityCurve, femaleCurve = femaleCurve)
## quick loop to test
indiv <- makeFounders(stocks = c(1))
for( y in 1:10) {
    indiv <- altMate( indiv, batchSize = batchSize, type = "ageSex",
                     maleCurve = maleCurve, femaleCurve = femaleCurve)
    indiv <- mort( indiv, year = y, type = "age", ageMort = c(ng, rep(mortRate, 30)))
    indiv <- birthdays(indiv)
}
outsize <- nrow(indiv[ is.na(indiv$DeathY),] )
expect_true( (800 < outsize) & 1200 > outsize)
