############################################
### Creating the .rda file from raw data ###
############################################

library(geomorph)
library(tidyverse)

dat.dorsal <- readland.tps("./extData/dorsal.TPS",
                           specID = "ID",
                           readcurves = FALSE)

dat.lateral <- readland.tps("./extData/lateral-unbent.TPS",
                           specID = "ID",
                           readcurves = FALSE)

gobyData <- read.csv("./extData/ethanolExperimentData.csv")

dorsalLinks <- read.csv("./extData/links.dorsal.csv")
lateralLinks <- read.csv("./extData/links.lateral.csv")



### Save as a .rda file ###
save(dat.dorsal, dat.lateral, gobyData, dorsalLinks, lateralLinks, file = "./data/gobyPreservation.rda")
