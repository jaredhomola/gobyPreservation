#################################################################
######       Ethanol experiment morphometric analyses      ######
######################     June  2020     #######################
#################################################################

library(geomorph)
library(raster)
library(gobyPreservation)
library(tidyverse)
data("gobyPreservation")


######## Nested ANOVA approach #########
#establish factors
pop = as.factor(c(rep(c("CEC"), 1620), rep(c("MSK"), 1620)))
ind = as.factor(c(rep(1:30, each = 18, times = 3), rep(1:30, each = 18, times = 3)))
rep = as.factor(c(rep(1:3, each = 540), rep(1:3, each = 540)))
day = as.factor(c(rep(rep(c(0, 1, 2, 7, 14, 21, 28, 35, 42, 49, 56, 63, 70, 77, 84, 105, 126, 154), times = 60), times = 3)))



###### Dorsal analysis #####
### Set data up for analysis ###
## Set up geomorph dataframe ##
proc.dorsal <- gpagen(dat.dorsal)

gdf.dorsal <- geomorph.data.frame(proc.dorsal,
                                  pop = pop,
                                  ind = ind,
                                  rep = rep,
                                  day = day,
                                  shape = proc.dorsal$coords)

## Remove allometry-associated variation
fit.simple.allo <- procD.lm(coords ~ log(Csize), 
                            data = gdf.dorsal, print.progress = FALSE)
shape.resid <- arrayspecs(fit.simple.allo$residuals,
                          p=dim(gdf.dorsal$coords)[1], k=dim(gdf.dorsal$coords)[2]) # allometry-adjusted residuals
adj.shape <- shape.resid + array(mshape(gdf.dorsal$coords), dim(shape.resid)) # Size adjusted shapes
proc.dorsal <- gpagen(adj.shape)

gdf.dorsal <- geomorph.data.frame(proc.dorsal,
                                  pop = pop,
                                  ind = ind,
                                  rep = rep,
                                  day = day,
                                  shape = proc.dorsal$coords)

## Use ANOVA to quantify measurement error
Nested_ANOVA.dorsal <- procD.lm(coords ~ pop/ind/day/rep, 
                                data = gdf.dorsal, 
                                seed = NULL, 
                                RRPP = TRUE, 
                                iter = 999)
summary(Nested_ANOVA.dorsal)


######## lateral ########
### Read data from fully appended tps file (created in TPSutil)
dat.lateral <- dat.lateral[c(1:17),,]

### Set data up for analysis ###
## Set up geomorph dataframe ##
proc.lateral <- gpagen(dat.lateral)

gdf.lateral <- geomorph.data.frame(proc.lateral,
                                  pop = pop,
                                  ind = ind,
                                  rep = rep,
                                  day = day,
                                  shape = proc.lateral$coords)

## Remove allometry-associated variation
fit.simple.allo <- procD.lm(coords ~ log(Csize), 
                            data = gdf.lateral, print.progress = FALSE)
shape.resid <- arrayspecs(fit.simple.allo$residuals,
                          p=dim(gdf.lateral$coords)[1], k=dim(gdf.lateral$coords)[2]) # allometry-adjusted residuals
adj.shape <- shape.resid + array(mshape(gdf.lateral$coords), dim(shape.resid)) # Size adjusted shapes
proc.lateral <- gpagen(adj.shape)


gdf.lateral <- geomorph.data.frame(proc.lateral,
                                pop = pop,
                                ind = ind,
                                rep = rep,
                                day = day,
                                shape = proc.lateral$coords)

## Use ANOVA to quantify measurement error
Nested_ANOVA.lateral <- procD.lm(coords ~ pop/ind/day/rep, 
                                 data = gdf.lateral, 
                                 seed = NULL, 
                                 RRPP = TRUE, 
                                 iter = 999)
summary(Nested_ANOVA.lateral)



######## PCA centroid approach #########
#### Dorsal ####
### Principal components analysis ###
pc.dorsal.df <- data.frame(two.d.array(proc.dorsal$coords))
pcaPoints.dorsal <- prcomp(pc.dorsal.df)$x

PCA.dorsal.df <- as.data.frame(pcaPoints.dorsal[,1:3]) # First three PCA axes
PCA.dorsal.df$ind <- ind
PCA.dorsal.df$rep <- rep
PCA.dorsal.df$day <- day
PCA.dorsal.df$pop <- pop

## Calculate population centroids by day
repCentroids <- aggregate(cbind(PC1, PC2) ~ day*ind*pop, PCA.dorsal.df, mean)

centDist.dorsal.tib <-
  left_join(repCentroids, PCA.dorsal.df, by = c("ind", "day", "pop")) %>%
  rename(PC1.cent = PC1.x,
    PC2.cent = PC2.x,
    PC1 = PC1.y,
    PC2 = PC2.y) %>%
  select(day, ind, rep, pop, PC1.cent, PC2.cent, PC1, PC2) %>%
  mutate(distCent = pointDistance(.[, 5:6], .[, 7:8], lonlat = FALSE)) %>%
  group_by(day, ind, pop) %>%
  summarize(meanDistCent.dorsal = mean(distCent))


#### lateral ####
### Principal components analysis ###
pc.lateral.df <- data.frame(two.d.array(proc.lateral$coords))
pcaPoints.lateral <- prcomp(pc.lateral.df)$x

PCA.lateral.df <- as.data.frame(pcaPoints.lateral[,1:3]) # First three PCA axes
PCA.lateral.df$ind <- ind
PCA.lateral.df$rep <- rep
PCA.lateral.df$day <- day
PCA.lateral.df$pop <- pop

## Calculate population centroids by day
repCentroids <- aggregate(cbind(PC1, PC2) ~ day*ind*pop, PCA.lateral.df, mean)

centDist.lateral.tib <-
  left_join(repCentroids, PCA.lateral.df, by = c("ind", "day", "pop")) %>%
  rename(PC1.cent = PC1.x,
         PC2.cent = PC2.x,
         PC1 = PC1.y,
         PC2 = PC2.y) %>%
  select(day, ind, rep, pop, PC1.cent, PC2.cent, PC1, PC2) %>%
  mutate(distCent = pointDistance(.[, 5:6], .[, 7:8], lonlat = FALSE)) %>%
  group_by(day, ind, pop) %>%
  summarize(meanDistCent.lateral = mean(distCent))

centDist.all.tib <- centDist.dorsal.tib %>% 
  left_join(centDist.lateral.tib) %>% 
  pivot_longer(cols = meanDistCent.dorsal:meanDistCent.lateral,
               values_to = "distCent",
               names_to = "perspective") %>% 
  mutate(perspective = recode(perspective, 
                              meanDistCent.dorsal = "Dorsal",
                              meanDistCent.lateral = "Lateral")) %>% 
  mutate(perspective = as.factor(perspective))
  
centDist.all.tib %>% 
  ggplot(aes(distCent)) +
  geom_histogram(bins = 60) +
  xlab("Mean distance to individual centroid in PCA space") +
  ylab("Frequency") +
  facet_grid(rows = vars(perspective)) +
  theme_bw()  


