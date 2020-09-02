#################################################################
######       Ethanol experiment morphometric analyses      ######
######################     June  2020     #######################
#################################################################

library(geomorph)
library(ggrepel)
library(plyr)
library(stringi)
library(gobyPreservation)
library(tidyverse)
data("gobyPreservation")

## Establish factors
pop = as.factor(c(rep(c("CEC"), 1620), rep(c("MSK"), 1620)))
ind = as.factor(c(rep(1:30, each = 18, times = 3), rep(1:30, each = 18, times = 3)))
rep = as.factor(c(rep(1:3, each = 540), rep(1:3, each = 540)))
day = as.factor(c(rep(rep(c(0, 1, 2, 7, 14, 21, 28, 35, 42, 49, 56, 63, 70, 77, 84, 105, 126, 154), times = 60), times = 3)))

######## Dorsal ########
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

## Population through time (First listed is gray, second listed is black)
plotRefToTarget(mean.cec0, mean.msk0, method="points", links = dorsalLinks[,2:3], mag = 2) 
plotRefToTarget(mean.cec14, mean.msk14, method="points", links = dorsalLinks[,2:3], mag = 2) 
plotRefToTarget(mean.cec154, mean.msk154, method="points", links = dorsalLinks[,2:3], mag = 2) 

## Within-population through time (First listed is gray, second listed is black)
plotRefToTarget(mean.cec0, mean.cec14, method="points", links = dorsalLinks[,2:3], mag = 1) 
plotRefToTarget(mean.cec14, mean.cec154, method="points", links = dorsalLinks[,2:3], mag = 1) 
plotRefToTarget(mean.cec0, mean.cec154, method="points", links = dorsalLinks[,2:3], mag = 1) 
plotRefToTarget(mean.msk0, mean.msk14, method="points", links = dorsalLinks[,2:3], mag = 1)
plotRefToTarget(mean.msk14, mean.msk154, method="points", links = dorsalLinks[,2:3], mag = 1)
plotRefToTarget(mean.msk0, mean.msk154, method="points", links = dorsalLinks[,2:3], mag = 1)

## Overall between population difference ##
plotRefToTarget(mean.cec, mean.msk, method="points", links = dorsalLinks[,2:3], mag = 2)




#### LATERAL #####
dat.lateral <- dat.lateral[c(1:17),,]

### Perform Procrutes transformation and construct data frame ###
proc.lateral <- gpagen(dat.lateral)
gdf.lateral <- geomorph.data.frame(proc.lateral, shape = proc.lateral$coords) 

## Perform Generalized Procrustes transformation and correct for allometric relationships ##
fit.simple.allo <- procD.lm(coords ~ log(Csize), 
                            data = gdf.lateral, print.progress = FALSE)
shape.resid <- arrayspecs(fit.simple.allo$residuals,
                          p=dim(gdf.lateral$coords)[1], k=dim(gdf.lateral$coords)[2]) # allometry-adjusted residuals
adj.shape <- shape.resid + array(mshape(gdf.lateral$coords), dim(shape.resid)) # Size adjusted shapes
proc.lateral <- gpagen(adj.shape)

lateral.df <- data.frame(two.d.array(proc.lateral$coords)) %>% 
  mutate(pop = pop, 
         ind = ind, 
         rep = rep, 
         day = day) %>% 
  group_by(pop, ind, day) %>% 
  summarize_if(is.numeric, funs(mean))

array.lateral <- arrayspecs(lateral.df[,4:37], 
                           17, 
                           2)

gdf.lateral <- geomorph.data.frame(shape = array.lateral,
                                  pop = lateral.df$pop,
                                  ind = lateral.df$ind,
                                  day = lateral.df$day) 

### Subset by pop and day
group <- factor(paste0(gdf.lateral$pop, gdf.lateral$day))
gdf.popDay <- coords.subset(gdf.lateral$shape, group = group)
pop <- factor(gdf.lateral$pop)
gdf.pop <- coords.subset(gdf.lateral$shape, group = pop)

### Calculate mean shapes
mean.cec0 <- mshape(gdf.popDay$CEC0)
mean.msk0 <- mshape(gdf.popDay$MSK0)
mean.cec14 <- mshape(gdf.popDay$CEC14)
mean.msk14 <- mshape(gdf.popDay$MSK14)
mean.cec154 <- mshape(gdf.popDay$CEC154)
mean.msk154 <- mshape(gdf.popDay$MSK154)

mean.cec <- mshape(gdf.pop$CEC)
mean.msk <- mshape(gdf.pop$MSK)


## Population through time (First listed is gray, second listed is black)
plotRefToTarget(mean.cec0, mean.msk0, method="points", links = lateralLinks[,2:3], mag = 2) 
plotRefToTarget(mean.cec14, mean.msk14, method="points", links = lateralLinks[,2:3], mag = 2) 
plotRefToTarget(mean.cec154, mean.msk154, method="points", links = lateralLinks[,2:3], mag = 2) 

## Within-population through time (First listed is gray, second listed is black)
plotRefToTarget(mean.cec0, mean.cec14, method="points", links = lateralLinks[,2:3], mag = 1) 
plotRefToTarget(mean.cec14, mean.cec154, method="points", links = lateralLinks[,2:3], mag = 1) 
plotRefToTarget(mean.cec0, mean.cec154, method="points", links = lateralLinks[,2:3], mag = 1) 
plotRefToTarget(mean.msk0, mean.msk14, method="points", links = lateralLinks[,2:3], mag = 1)
plotRefToTarget(mean.msk14, mean.msk154, method="points", links = lateralLinks[,2:3], mag = 1)
plotRefToTarget(mean.msk0, mean.msk154, method="points", links = lateralLinks[,2:3], mag = 1)

## Overall between population difference ##
plotRefToTarget(mean.cec, mean.msk, method="points", links = lateralLinks[,2:3], mag = 2)
