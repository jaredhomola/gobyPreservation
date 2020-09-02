#################################################################
######       Ethanol experiment morphometric analyses      ######
######################     June  2020     #######################
#################################################################

library(geomorph)
library(raster)
library(tidyverse)


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
Nested_ANOVA.dorsal <- procD.lm(coords ~ pop/ind/day/rep, data = gdf.dorsal, seed = NULL, RRPP = TRUE, iter=999)
summary(Nested_ANOVA.dorsal)


######## Side ########
### Read data from fully appended tps file (created in TPSutil)
dat.side <- dat.side[c(1:17),,]

### Set data up for analysis ###
## Set up geomorph dataframe ##
proc.side <- gpagen(dat.side)

gdf.side <- geomorph.data.frame(proc.side,
                                  pop = pop,
                                  ind = ind,
                                  rep = rep,
                                  day = day,
                                  shape = proc.side$coords)

## Remove allometry-associated variation
fit.simple.allo <- procD.lm(coords ~ log(Csize), 
                            data = gdf.side, print.progress = FALSE)
shape.resid <- arrayspecs(fit.simple.allo$residuals,
                          p=dim(gdf.side$coords)[1], k=dim(gdf.side$coords)[2]) # allometry-adjusted residuals
adj.shape <- shape.resid + array(mshape(gdf.side$coords), dim(shape.resid)) # Size adjusted shapes
proc.side <- gpagen(adj.shape)


gdf.side <- geomorph.data.frame(proc.side,
                                pop = pop,
                                ind = ind,
                                rep = rep,
                                day = day,
                                shape = proc.side$coords)

## Use ANOVA to quantify measurement error
Nested_ANOVA.side <- procD.lm(coords ~ pop/ind/day/rep, data = gdf.side, seed = NULL, RRPP = TRUE, iter=999)
summary(Nested_ANOVA.side)



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


#### Side ####
### Principal components analysis ###
pc.side.df <- data.frame(two.d.array(proc.side$coords))
pcaPoints.side <- prcomp(pc.side.df)$x

PCA.side.df <- as.data.frame(pcaPoints.side[,1:3]) # First three PCA axes
PCA.side.df$ind <- ind
PCA.side.df$rep <- rep
PCA.side.df$day <- day
PCA.side.df$pop <- pop

## Calculate population centroids by day
repCentroids <- aggregate(cbind(PC1, PC2) ~ day*ind*pop, PCA.side.df, mean)

centDist.side.tib <-
  left_join(repCentroids, PCA.side.df, by = c("ind", "day", "pop")) %>%
  rename(PC1.cent = PC1.x,
         PC2.cent = PC2.x,
         PC1 = PC1.y,
         PC2 = PC2.y) %>%
  select(day, ind, rep, pop, PC1.cent, PC2.cent, PC1, PC2) %>%
  mutate(distCent = pointDistance(.[, 5:6], .[, 7:8], lonlat = FALSE)) %>%
  group_by(day, ind, pop) %>%
  summarize(meanDistCent.side = mean(distCent))

centDist.all.tib <- centDist.dorsal.tib %>% 
  left_join(centDist.side.tib) %>% 
  pivot_longer(cols = meanDistCent.dorsal:meanDistCent.side,
               values_to = "distCent",
               names_to = "perspective") %>% 
  mutate(perspective = recode(perspective, 
                              meanDistCent.dorsal = "Dorsal",
                              meanDistCent.side = "Lateral")) %>% 
  mutate(perspective = as.factor(perspective))
  
centDist.all.tib %>% 
  ggplot(aes(distCent)) +
  geom_histogram(bins = 60) +
  xlab("Mean distance to individual centroid in PCA space") +
  ylab("Frequency") +
  facet_grid(rows = vars(perspective)) +
  theme_bw()  

ggsave("./Ethanol Experiment/plots/errorDistribution.pdf",
       width = 7, height = 6, units = "in")
