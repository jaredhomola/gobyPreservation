############################################
### Creating the .rda file from raw data ###
############################################

library(geomorph)
library(tidyverse)

dat.dorsal <- readland.tps("./Ethanol Experiment/gobyPreservation/extData/dorsal.TPS",
                           specID = "ID",
                           readcurves = FALSE)

dat.lateral <- readland.tps("./Ethanol Experiment/gobyPreservation/extData/lateral-unbent.TPS",
                           specID = "ID",
                           readcurves = FALSE)

gobyData <- read.csv("./Ethanol Experiment/gobyPreservation/extData/ethanolExperimentData.csv")

dorsalLinks <- read.csv("./Ethanol Experiment/gobyPreservation/extData/links.dorsal.csv")
lateralLinks <- read.csv("./Ethanol Experiment/gobyPreservation/extData/links.lateral.csv")



### Save as a .rda file ###
save(dat.dorsal, dat.lateral, gobyData, dorsalLinks, lateralLinks, file = "./Ethanol Experiment/gobyPreservation/data/gobyPreservation.rda")
