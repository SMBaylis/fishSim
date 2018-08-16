## this is a script to make and document the current dev version of fishSim

name <- "fishSim0.0.0.8000" ## update this, obvs
# SB_description <- list("maintainer" = "'Shane Baylis' <shane.baylis@csiro.au>",
                       
library(devtools)
library(roxygen2)

setwd("~/Dropbox/R/localPackages/fishSim")
create(name)

## now, copy over your latest dev functions into "./[version]/R/"
currentfiles <- "~/Dropbox/CSIRO/Simulations/fishSim_dev.R"
newlocation <- paste("./",name,"/R", sep = "")

file.copy(from = currentfiles, to = newlocation,
          overwrite = TRUE)

## and then document the new version.
setwd(paste("./",name,sep = ""))
document()
