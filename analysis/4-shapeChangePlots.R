#################################################################
######       Ethanol experiment morphometric analyses      ######
######################     June  2020     #######################
#################################################################

library(geomorph)
library(ggrepel)
library(plyr)
library(stringi)
library(tidyverse)

## Establish factors
pop = as.factor(c(rep(c("CEC"), 1620), rep(c("MSK"), 1620)))
ind = as.factor(c(rep(1:30, each = 18, times = 3), rep(1:30, each = 18, times = 3)))
rep = as.factor(c(rep(1:3, each = 540), rep(1:3, each = 540)))
day = as.factor(c(rep(rep(c(0, 1, 2, 7, 14, 21, 28, 35, 42, 49, 56, 63, 70, 77, 84, 105, 126, 154), times = 60), times = 3)))

######## Dorsal ########
### Read data from fully appended tps file (created in TPSutil)
dat.dorsal <-
  readland.tps("./Ethanol Experiment/Experiment Photos/Experiment-photos-dorsal/individual-dorsal-appended/individual-appended-dorsal-all.TPS",
               specID = "ID",
               readcurves = FALSE)

### Perform Procrutes transformation and construct data frame ###
proc.dorsal <- gpagen(dat.dorsal)
gdf.dorsal <- geomorph.data.frame(proc.dorsal, shape = proc.dorsal$coords) 

## Perform Generalized Procrustes transformation and correct for allometric relationships ##
fit.simple.allo <- procD.lm(coords ~ log(Csize), 
                            data = gdf.dorsal, print.progress = FALSE)
shape.resid <- arrayspecs(fit.simple.allo$residuals,
                          p=dim(gdf.dorsal$coords)[1], k=dim(gdf.dorsal$coords)[2]) # allometry-adjusted residuals
adj.shape <- shape.resid + array(mshape(gdf.dorsal$coords), dim(shape.resid)) # Size adjusted shapes
proc.dorsal <- gpagen(adj.shape)

dorsal.df <- data.frame(two.d.array(proc.dorsal$coords)) %>% 
  mutate(pop = pop, 
         ind = ind, 
         rep = rep, 
         day = day) %>% 
  group_by(pop, ind, day) %>% 
  summarize_if(is.numeric, funs(mean))

array.dorsal <- arrayspecs(dorsal.df[,4:15], 
                           6, 
                           2)

gdf.dorsal <- geomorph.data.frame(shape = array.dorsal,
                                  pop = dorsal.df$pop,
                                  ind = dorsal.df$ind,
                                  day = dorsal.df$day) 

### Subset by pop and day
group <- factor(paste0(gdf.dorsal$pop, gdf.dorsal$day))
gdf.popDay <- coords.subset(gdf.dorsal$shape, group = group)
pop <- factor(gdf.dorsal$pop)
gdf.pop <- coords.subset(gdf.dorsal$shape, group = pop)


### Calculate mean shapes
mean.cec0 <- mshape(gdf.popDay$CEC0)
mean.msk0 <- mshape(gdf.popDay$MSK0)
mean.cec14 <- mshape(gdf.popDay$CEC14)
mean.msk14 <- mshape(gdf.popDay$MSK14)
mean.cec154 <- mshape(gdf.popDay$CEC154)
mean.msk154 <- mshape(gdf.popDay$MSK154)

mean.cec <- mshape(gdf.pop$CEC)
mean.msk <- mshape(gdf.pop$MSK)

## Define links
#links.dorsal <- define.links(mean.cec0)
#write.csv(links.dorsal, "links.dorsal2.csv")
links.dorsal <- read.csv("links.dorsal2.csv")[,2:3]

## Population through time (First listed is gray, second listed is black)
plotRefToTarget(mean.cec0, mean.msk0, method="points", links = links.dorsal, mag = 2) 
plotRefToTarget(mean.cec14, mean.msk14, method="points", links = links.dorsal, mag = 2) 
plotRefToTarget(mean.cec154, mean.msk154, method="points", links = links.dorsal, mag = 2) 

## Within-population through time (First listed is gray, second listed is black)
plotRefToTarget(mean.cec0, mean.cec14, method="points", links = links.dorsal, mag = 1) 
plotRefToTarget(mean.cec14, mean.cec154, method="points", links = links.dorsal, mag = 1) 
plotRefToTarget(mean.cec0, mean.cec154, method="points", links = links.dorsal, mag = 1) 
plotRefToTarget(mean.msk0, mean.msk14, method="points", links = links.dorsal, mag = 1)
plotRefToTarget(mean.msk14, mean.msk154, method="points", links = links.dorsal, mag = 1)
plotRefToTarget(mean.msk0, mean.msk154, method="points", links = links.dorsal, mag = 1)

## Overall between population difference ##
plotRefToTarget(mean.cec, mean.msk, method="points", links = links.dorsal, mag = 2)




#### LATERAL #####
dat.side <-
  readland.tps("./Ethanol Experiment/Experiment Photos/Experiment-photos-side-UNBEND/unbentAll.TPS",
               specID = "ID",
               readcurves = FALSE)

dat.side <- dat.side[c(1:17),,]

### Perform Procrutes transformation and construct data frame ###
proc.side <- gpagen(dat.side)
gdf.side <- geomorph.data.frame(proc.side, shape = proc.side$coords) 

## Perform Generalized Procrustes transformation and correct for allometric relationships ##
fit.simple.allo <- procD.lm(coords ~ log(Csize), 
                            data = gdf.side, print.progress = FALSE)
shape.resid <- arrayspecs(fit.simple.allo$residuals,
                          p=dim(gdf.side$coords)[1], k=dim(gdf.side$coords)[2]) # allometry-adjusted residuals
adj.shape <- shape.resid + array(mshape(gdf.side$coords), dim(shape.resid)) # Size adjusted shapes
proc.side <- gpagen(adj.shape)

side.df <- data.frame(two.d.array(proc.side$coords)) %>% 
  mutate(pop = pop, 
         ind = ind, 
         rep = rep, 
         day = day) %>% 
  group_by(pop, ind, day) %>% 
  summarize_if(is.numeric, funs(mean))

array.side <- arrayspecs(side.df[,4:37], 
                           17, 
                           2)

gdf.side <- geomorph.data.frame(shape = array.side,
                                  pop = side.df$pop,
                                  ind = side.df$ind,
                                  day = side.df$day) 

### Subset by pop and day
group <- factor(paste0(gdf.side$pop, gdf.side$day))
gdf.popDay <- coords.subset(gdf.side$shape, group = group)
pop <- factor(gdf.side$pop)
gdf.pop <- coords.subset(gdf.side$shape, group = pop)

### Calculate mean shapes
mean.cec0 <- mshape(gdf.popDay$CEC0)
mean.msk0 <- mshape(gdf.popDay$MSK0)
mean.cec14 <- mshape(gdf.popDay$CEC14)
mean.msk14 <- mshape(gdf.popDay$MSK14)
mean.cec154 <- mshape(gdf.popDay$CEC154)
mean.msk154 <- mshape(gdf.popDay$MSK154)

mean.cec <- mshape(gdf.pop$CEC)
mean.msk <- mshape(gdf.pop$MSK)

## Define links
#links.side <- define.links(mean.cec0)
#write.csv(links.side, "links.side.csv")
links.side <- read.csv("links.side.csv")[,2:3]


## Population through time (First listed is gray, second listed is black)
plotRefToTarget(mean.cec0, mean.msk0, method="points", links = links.side, mag = 2) 
plotRefToTarget(mean.cec14, mean.msk14, method="points", links = links.side, mag = 2) 
plotRefToTarget(mean.cec154, mean.msk154, method="points", links = links.side, mag = 2) 

## Within-population through time (First listed is gray, second listed is black)
plotRefToTarget(mean.cec0, mean.cec14, method="points", links = links.side, mag = 1) 
plotRefToTarget(mean.cec14, mean.cec154, method="points", links = links.side, mag = 1) 
plotRefToTarget(mean.cec0, mean.cec154, method="points", links = links.side, mag = 1) 
plotRefToTarget(mean.msk0, mean.msk14, method="points", links = links.side, mag = 1)
plotRefToTarget(mean.msk14, mean.msk154, method="points", links = links.side, mag = 1)
plotRefToTarget(mean.msk0, mean.msk154, method="points", links = links.side, mag = 1)

## Overall between population difference ##
plotRefToTarget(mean.cec, mean.msk, method="points", links = links.side, mag = 2)
