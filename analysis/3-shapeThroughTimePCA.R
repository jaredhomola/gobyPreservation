#################################################################
######       Ethanol experiment morphometric analyses      ######
######################     June  2020     #######################
#################################################################

library(geomorph)
library(ggrepel)
library(raster)
library(abind)
library(GeometricMorphometricsMix)
library(gobyPreservation)
library(tidyverse)
data("gobyPreservation")


## Establish factors
pop = as.factor(c(rep(c("CEC"), 1620), rep(c("MSK"), 1620)))
ind = as.factor(c(rep(1:30, each = 18, times = 3), rep(31:60, each = 18, times = 3)))
rep = as.factor(c(rep(1:3, each = 540), rep(1:3, each = 540)))
day = as.factor(c(rep(rep(c(0, 1, 2, 7, 14, 21, 28, 35, 42, 49, 56, 63, 70, 77, 84, 105, 126, 154), times = 60), times = 3)))

day.vec <- c(0, 1, 2, 7, 14, 21, 28, 35, 42, 49, 56, 63, 70, 77, 84, 105, 126, 154)


###### Dorsal analysis #####
## Perform Generalized Procrustes transformation
proc.dorsal <- gpagen(dat.dorsal)

gdf.dorsal <- geomorph.data.frame(proc.dorsal,
                                  shape = proc.dorsal$coords) 
fit.simple.allo <- procD.lm(coords ~ log(Csize), 
                            data = gdf.dorsal, print.progress = FALSE)
shape.resid <- arrayspecs(fit.simple.allo$residuals,
                          p=dim(gdf.dorsal$coords)[1], k=dim(gdf.dorsal$coords)[2]) # allometry-adjusted residuals
adj.shape <- shape.resid + array(mshape(gdf.dorsal$coords), dim(shape.resid)) # Size adjusted shapes
proc.dorsal <- gpagen(adj.shape)

## Check for outliers
plotOutliers(proc.dorsal$coords)

### Principal components analysis ###
## Set up data frame and perform analysis
dorsal.df <- data.frame(two.d.array(proc.dorsal$coords)) %>% 
  mutate(pop = pop, 
         ind = ind, 
         rep = rep, 
         day = day) %>% 
  group_by(pop, ind, day) %>% 
  summarize_if(is.numeric, funs(mean))

pcaObj <- prcomp(dorsal.df[,4:15])
pcaPoints.dorsal <- pcaObj$x
summary(pcaObj)

PCA.dorsal.df <- as.data.frame(pcaPoints.dorsal[,1:3]) # First three PCA axes
PCA.dorsal.df$ind <- dorsal.df$ind
PCA.dorsal.df$day <- dorsal.df$day
PCA.dorsal.df$pop <- dorsal.df$pop

## Calculate population centroids by day
centroids.dorsal.12 <- aggregate(cbind(PC1, PC2) ~ day*pop, PCA.dorsal.df, mean)
centroids.dorsal.23 <- aggregate(cbind(PC2, PC3) ~ day*pop, PCA.dorsal.df, mean)

## PCA plot by population - overall
ggplot(PCA.dorsal.df, aes(x = PC1, y = PC2, color = as.factor(day))) +
  geom_point(alpha = .2, aes(shape = pop)) +
  geom_point(data = centroids.dorsal.12, size = 5, aes(shape = pop)) +
  stat_ellipse(alpha = 0.4) +
  geom_text_repel(data=centroids.dorsal.12, fontface = "bold", aes(label=day), size=6, box.padding = .4) +
  ylab("PC2 (20.1%)") +
  xlab("PC1 (45.9%)") +
  theme_classic(base_size = 18) +
  theme(legend.position = "none") +
  NULL

## PCA plot by population - centroids only
ggplot(PCA.dorsal.df, aes(x = PC1, y = PC2)) +
  geom_point(data = centroids.dorsal.12, size = 5, aes(shape = pop)) +
  geom_text_repel(data=centroids.dorsal.12, fontface = "bold", aes(label=day), size=6, box.padding = .4) +
  ylab("PC2 (20.1%)") +
  xlab("PC1 (45.9%)") +
  theme_classic(base_size = 18) +
  theme(legend.position = "none") +
  NULL

## Incremental centroid distance changes
CEC.cents <- centroids.dorsal.12[1:18,]
MSK.cents <- centroids.dorsal.12[19:36,]

CEC.centDist.tmp <- as.matrix(dist(CEC.cents[-1]))
MSK.centDist.tmp <- as.matrix(dist(MSK.cents[-1]))

CEC.centDist <- as.data.frame(CEC.centDist.tmp[row(CEC.centDist.tmp) == (col(CEC.centDist.tmp) - 1)])
MSK.centDist <- as.data.frame(MSK.centDist.tmp[row(MSK.centDist.tmp) == (col(MSK.centDist.tmp) - 1)])

dist.df <- cbind(CEC.centDist, MSK.centDist)
names(dist.df) <- c("cecDist", "mskDist")
dist.df$timePeriod <- c(1, 2, 7, 14, 21, 28, 35, 42, 49, 56, 63, 70, 77, 84, 105, 126, 154)
dist.df.dorsal <- dist.df %>%
  mutate(
    dayMove.dorsal.cec = cecDist / (timePeriod - lag(timePeriod, default = 0)),
    dayMove.dorsal.msk = mskDist / (timePeriod - lag(timePeriod, default = 0))
  ) %>% 
  select(-cecDist, -mskDist)

## Between-population centroid distance by day
distCent.dorsal <- as.data.frame(pointDistance(CEC.cents[, 3:4], 
                                        MSK.cents[, 3:4], 
                                        lonlat = FALSE))
names(distCent.dorsal) <- "distDorsal"
distCent.dorsal$timePeriod <- day.vec




#### Repeated measures modeling for geometric morphometrics
varPairs <- combn(unique(dorsal.df$day), 2, simplify = TRUE) %>%
  as.character() %>% as.data.frame() %>%
  group_by(grp = str_c('Column', rep(1:2, length.out = n()))) %>%
  mutate(rn = row_number()) %>%
  rename(value = 1) %>% 
  ungroup %>%
  pivot_wider(names_from = grp, values_from = value) %>%
  select(-rn)

gobyRM <- function(pop1, pop2, day1, day2, ...) {
  dat1 <- pcaPoints.dorsal %>% 
    as_tibble() %>% 
    mutate(pop = dorsal.df$pop,
           ind = dorsal.df$ind,
           day = dorsal.df$day) %>% 
    filter(pop == pop1, day == day1) %>% 
    ungroup() %>% 
    select(-c(pop, ind, day)) %>% 
    select(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8)
  
  dat2 <- pcaPoints.dorsal %>% 
    as_tibble() %>% 
    mutate(pop = dorsal.df$pop,
           ind = dorsal.df$ind,
           day = dorsal.df$day) %>% 
    filter(pop == pop2, day == day2) %>% 
    ungroup() %>% 
    select(-c(pop, ind, day)) %>% 
    select(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8)
  
  suppressWarnings(repeated_measures_test(dat1, dat2, rnames = FALSE))
}

results <- tibble(1:nrow(varPairs)) %>% 
  rename(EuclideanD = `1:nrow(varPairs)`) %>% 
  mutate(HotellingT2 = as.numeric(as.character(varPairs$Column1))) %>% 
  mutate(Fstat = as.numeric(as.character(varPairs$Column1))) %>% 
  mutate(p_val = as.numeric(as.character(varPairs$Column1))) %>% 
  mutate(time1 = as.numeric(as.character(varPairs$Column1))) %>% 
  mutate(time2 = as.numeric(as.character(varPairs$Column1)))
  
for(z in (1:nrow(varPairs))) {
  RM.out <- gobyRM("MSK", "CEC", varPairs$Column1[z], varPairs$Column2[z])
  results$EuclideanD[z] <- RM.out %>% as.data.frame %>% .[1,]
  results$HotellingT2[z] <- RM.out %>% as.data.frame %>% .[2,]
  results$Fstat[z] <- RM.out %>% as.data.frame %>% .[3,]
  results$p_val[z] <- RM.out %>% as.data.frame %>% .[4,]
  results$time1[z] <- varPairs$Column1[z]
  results$time2[z] <- varPairs$Column2[z]
}

results %>% 
  filter(time1 == 0) %>% 
  ggplot(aes(x = time2, y = EuclideanD)) +
  geom_point()

results <- results %>% 
  mutate(adj.p.value = p.adjust(p_val, method='fdr', n = nrow(.)))
  

## Comparisons within MSK and CEC
for(z in (1:nrow(varPairs))) {
  RM.out <- gobyRM("CEC", "CEC", varPairs$Column1[z], varPairs$Column2[z])
  results$EuclideanD.CEC[z] <- RM.out %>% as.data.frame %>% .[1,]
  results$HotellingT2.CEC[z] <- RM.out %>% as.data.frame %>% .[2,]
  results$Fstat.CEC[z] <- RM.out %>% as.data.frame %>% .[3,]
  results$p_val.CEC[z] <- RM.out %>% as.data.frame %>% .[4,]
  results$time1[z] <- varPairs$Column1[z]
  results$time2[z] <- varPairs$Column2[z]
}
results.CEC <- results

## Comparisons within MSK and CEC
for(z in (1:nrow(varPairs))) {
  RM.out <- gobyRM("MSK", "MSK", varPairs$Column1[z], varPairs$Column2[z])
  results$EuclideanD.MSK[z] <- RM.out %>% as.data.frame %>% .[1,]
  results$HotellingT2.MSK[z] <- RM.out %>% as.data.frame %>% .[2,]
  results$Fstat.MSK[z] <- RM.out %>% as.data.frame %>% .[3,]
  results$p_val.MSK[z] <- RM.out %>% as.data.frame %>% .[4,]
  results$time1[z] <- varPairs$Column1[z]
  results$time2[z] <- varPairs$Column2[z]
}
results.MSK <- results

results.CEC.MSK <- left_join(results.CEC, results.MSK) %>% 
  distinct(time1, .keep_all = TRUE) %>% 
  mutate(CEC.adj.p.value = p.adjust(p_val.CEC, method='fdr', n = nrow(.))) %>%
  mutate(MSK.adj.p.value = p.adjust(p_val.MSK, method='fdr', n = nrow(.)))

write.csv(results.CEC.MSK, "./Ethanol Experiment/RM.withinpop-dorsal.csv")


## Between-population modeling per day 
#### Question: How would a Procrustes ANOVA that compares that populations change across time?
array.dorsal <- arrayspecs(dorsal.df[,4:15], 6, 2)
gdf <- geomorph.data.frame(coords = array.dorsal, 
                           pop = dorsal.df$pop, 
                           day = dorsal.df$day)

varPairs.complete <- unique(dorsal.df$day) %>%
  as_tibble %>% 
  rename(Column1 = value) %>% 
  mutate(Column2 = unique(dorsal.df$day)) %>% 
  rbind(varPairs)

stats <- as_tibble(varPairs.complete) %>% 
  mutate(pval = NA,
         r2 = NA,
         Z = NA)
i = 0

## For all day-population combinations
for(z in (1:nrow(varPairs.complete))) {
  i = i+1
  coords1 <- gdf$coords[,,gdf$day == varPairs.complete$Column1[z]][,,1:30] ### CEC is column 1
  coords2 <- gdf$coords[,,gdf$day == varPairs.complete$Column2[z]][,,31:60] ### MSK is column 2
  coords.both <- abind(coords1, coords2, along = 3)

  pop.test <- as.factor(c(rep("CEC", 30), rep("MSK", 30)))
  gdf.test <- geomorph.data.frame(coords = coords.both, pop = pop.test)
  
  fit <- procD.lm(coords ~ pop,
                  data = gdf.test)
  
  stats$pval[i] <- fit$aov.table$`Pr(>F)`[1]
  stats$r2[i] <- fit$aov.table$Rsq[1]
  stats$Z[i] <- fit$aov.table$Z[1]
}

stats <- stats %>% 
  mutate(qval = p.adjust(pval, method='fdr', n = nrow(.)))
  

### Full matrix
stats1 <- stats %>% select(Column1, Column2, Z)
Zscores <- stats1 %>% pivot_wider(names_from = Column1, values_from = Z)




##### Lateral landmarks #####
## Read data from fully appended tps file (created in TPSutil)
dat.lateral <- dat.lateral[c(1:17),,]

## Perform Generalized Procrustes transformation
proc.lateral <- gpagen(dat.lateral)
gdf.lateral <- geomorph.data.frame(proc.lateral,
                                  shape = proc.lateral$coords) 
fit.simple.allo <- procD.lm(coords ~ log(Csize), 
                            data = gdf.lateral, print.progress = FALSE)
shape.resid <- arrayspecs(fit.simple.allo$residuals,
                          p=dim(gdf.lateral$coords)[1], k=dim(gdf.lateral$coords)[2]) # allometry-adjusted residuals
adj.shape <- shape.resid + array(mshape(gdf.lateral$coords), dim(shape.resid)) # Size adjusted shapes
proc.lateral <- gpagen(adj.shape)

## Check for outliers
plotOutliers(proc.lateral$coords)

### Principal components analysis ###
## Set up data frame and perform analysis
lateral.df <- data.frame(two.d.array(proc.lateral$coords)) %>% 
  mutate(pop = pop, 
         ind = ind, 
         rep = rep, 
         day = day) %>% 
  group_by(pop, ind, day) %>% 
  summarize_if(is.numeric, funs(mean))

pcaObj <- prcomp(lateral.df[,4:15])
pcaPoints.lateral <- pcaObj$x
summary(pcaObj)

PCA.lateral.df <- as.data.frame(pcaPoints.lateral[,1:3]) # First three PCA axes
PCA.lateral.df$ind <- lateral.df$ind
PCA.lateral.df$day <- lateral.df$day
PCA.lateral.df$pop <- lateral.df$pop

## Calculate population centroids by day
centroids.lateral.12 <- aggregate(cbind(PC1, PC2) ~ day*pop, PCA.lateral.df, mean)
centroids.lateral.23 <- aggregate(cbind(PC2, PC3) ~ day*pop, PCA.lateral.df, mean)

## PCA plot by population - overall
ggplot(PCA.lateral.df, aes(x = PC1, y = PC2, color = as.factor(day))) +
  geom_point(alpha = .2, aes(shape = pop)) +
  geom_point(data = centroids.lateral.12, size = 5, aes(shape = pop)) +
  stat_ellipse(alpha = 0.4) +
  geom_text_repel(data=centroids.lateral.12, fontface = "bold", aes(label=day), size=6, box.padding = .4) +
  ylab("PC2 (21.3%)") +
  xlab("PC1 (33.2%)") +
  theme_classic(base_size = 18) +
  theme(legend.position = "none") +
  NULL

## PCA plot by population - centroids only
ggplot(PCA.lateral.df, aes(x = PC1, y = PC2)) +
  geom_point(data = centroids.lateral.12, size = 5, aes(shape = pop)) +
  geom_text_repel(data=centroids.lateral.12, fontface = "bold", aes(label=day), size=6, box.padding = .4) +
  ylab("PC2 (21.3%)") +
  xlab("PC1 (33.2%)") +
  theme_classic(base_size = 18) +
  theme(legend.position = "none") +
  NULL

## Incremental centroid distance changes
CEC.cents <- centroids.lateral.12[1:18,]
MSK.cents <- centroids.lateral.12[19:36,]

CEC.centDist.tmp <- as.matrix(dist(CEC.cents[-1]))
MSK.centDist.tmp <- as.matrix(dist(MSK.cents[-1]))

CEC.centDist <- as.data.frame(CEC.centDist.tmp[row(CEC.centDist.tmp) == (col(CEC.centDist.tmp) - 1)])
MSK.centDist <- as.data.frame(MSK.centDist.tmp[row(MSK.centDist.tmp) == (col(MSK.centDist.tmp) - 1)])

dist.df <- cbind(CEC.centDist, MSK.centDist)
names(dist.df) <- c("cecDist", "mskDist")
dist.df$timePeriod <- c(1, 2, 7, 14, 21, 28, 35, 42, 49, 56, 63, 70, 77, 84, 105, 126, 154)
dist.df.tmp <- dist.df %>%
  mutate(
    dayMove.lateral.cec = cecDist / (timePeriod - lag(timePeriod, default = 0)),
    dayMove.lateral.msk = mskDist / (timePeriod - lag(timePeriod, default = 0))
  ) %>% 
  select(-cecDist, -mskDist) %>% 
  left_join(dist.df.dorsal) %>% 
  pivot_longer(cols = dayMove.lateral.cec:dayMove.lateral.msk,
               values_to = "dayMove",
               names_to = "popView") %>% 
  separate(popView, into = c(NA, "view", "pop")) %>% 
  pivot_longer(cols = dayMove.dorsal.cec:dayMove.dorsal.msk,
               values_to = "dayMove2",
               names_to = "popView") %>% 
  separate(popView, into = c(NA, "view2", "pop2"))

dist.df.all <- as_tibble(c(dist.df.tmp$timePeriod, dist.df.tmp$timePeriod)) %>% 
  mutate(dayMove = c(dist.df.tmp$dayMove, dist.df.tmp$dayMove2)) %>% 
  mutate(pop = c(dist.df.tmp$pop, dist.df.tmp$pop2)) %>% 
  mutate(view = c(dist.df.tmp$view, dist.df.tmp$view2)) %>% 
  distinct(value, dayMove, pop, view) %>% 
  mutate(view = recode(view, 
                       dorsal = "Dorsal",
                       lateral = "Lateral")) 

dist.df.all %>% 
  ggplot(aes(y = dayMove, x = value, shape = pop, linetype = pop)) +
  geom_line(size = 0.4) +
  geom_point(size = 2) +
  scale_linetype_discrete(breaks=c("A", "C")) +
  ylab("PCA centroid distance moved per day") +
  xlab("Time period endpoint (days)") +
  facet_grid(rows = vars(view)) +
  geom_hline(yintercept = 0, size = 0.1) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none") +
  NULL

## Between-population centroid distance by day
distCent.lateral <- as.data.frame(pointDistance(CEC.cents[, 3:4], 
                                        MSK.cents[, 3:4], 
                                        lonlat = FALSE))
names(distCent.lateral) <- "distlateral"
distCent.lateral$timePeriod <- day.vec

distCent.all <- left_join(distCent.dorsal, distCent.lateral) %>% 
  pivot_longer(cols = c("distDorsal", "distlateral"),
               names_to = "view",
               values_to = "dist") %>% 
  mutate(view = recode(view, 
                       distDorsal = "Dorsal",
                       distlateral = "Lateral")) 

distCent.all %>% 
  ggplot(aes(y = dist, x = timePeriod)) +
  geom_smooth(span = 1, col = "black") +
  geom_point(size = 2) +
  scale_linetype_discrete(breaks=c("A", "C")) +
  ylab("LkMi vs. CC centroid distance") +
  xlab("Days post-preservation") +
  facet_grid(rows = vars(view),
             scales = "free_y") +
  theme_bw(base_size = 14) +
  NULL

## Between-population modeling per day 
varPairs <- combn(unique(lateral.df$day), 2, simplify = TRUE) %>%
  as.character() %>% as.data.frame() %>%
  group_by(grp = str_c('Column', rep(1:2, length.out = n()))) %>%
  mutate(rn = row_number()) %>%
  rename(value = 1) %>% 
  ungroup %>%
  pivot_wider(names_from = grp, values_from = value) %>%
  select(-rn)

gobyRM <- function(pop1, pop2, day1, day2, ...) {
  dat1 <- pcaPoints.lateral %>% 
    as_tibble() %>% 
    mutate(pop = lateral.df$pop,
           ind = lateral.df$ind,
           day = lateral.df$day) %>% 
    filter(pop == pop1, day == day1) %>% 
    ungroup() %>% 
    select(-c(pop, ind, day)) %>% 
    select(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8)
  
  dat2 <- pcaPoints.lateral %>% 
    as_tibble() %>% 
    mutate(pop = lateral.df$pop,
           ind = lateral.df$ind,
           day = lateral.df$day) %>% 
    filter(pop == pop2, day == day2) %>% 
    ungroup() %>% 
    select(-c(pop, ind, day)) %>% 
    select(PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8)
  
  suppressWarnings(repeated_measures_test(dat1, dat2, rnames = FALSE))
}

results <- tibble(1:nrow(varPairs)) %>% 
  rename(EuclideanD = `1:nrow(varPairs)`) %>% 
  mutate(HotellingT2 = as.numeric(as.character(varPairs$Column1))) %>% 
  mutate(Fstat = as.numeric(as.character(varPairs$Column1))) %>% 
  mutate(p_val = as.numeric(as.character(varPairs$Column1))) %>% 
  mutate(time1 = as.numeric(as.character(varPairs$Column1))) #%>% 
# mutate(time2 = as.numeric(as.character(varPairs$Column1)))

## Comparisons within CEC and MSK
for(z in (1:nrow(varPairs))) {
  RM.out <- gobyRM("CEC", "CEC", varPairs$Column1[z], varPairs$Column2[z])
  results$EuclideanD.CEC[z] <- RM.out %>% as.data.frame %>% .[1,]
  results$HotellingT2.CEC[z] <- RM.out %>% as.data.frame %>% .[2,]
  results$Fstat.CEC[z] <- RM.out %>% as.data.frame %>% .[3,]
  results$p_val.CEC[z] <- RM.out %>% as.data.frame %>% .[4,]
  results$time1[z] <- varPairs$Column1[z]
  results$time2[z] <- varPairs$Column2[z]
}
results.CEC <- results

for(z in (1:nrow(varPairs))) {
  RM.out <- gobyRM("MSK", "MSK", varPairs$Column1[z], varPairs$Column2[z])
  results$EuclideanD.MSK[z] <- RM.out %>% as.data.frame %>% .[1,]
  results$HotellingT2.MSK[z] <- RM.out %>% as.data.frame %>% .[2,]
  results$Fstat.MSK[z] <- RM.out %>% as.data.frame %>% .[3,]
  results$p_val.MSK[z] <- RM.out %>% as.data.frame %>% .[4,]
  results$time1[z] <- varPairs$Column1[z]
  results$time2[z] <- varPairs$Column2[z]
}
results.MSK <- results

results.CEC.MSK <- left_join(results.CEC, results.MSK) %>% 
  distinct(time1, .keep_all = TRUE) %>% 
  mutate(CEC.adj.p.value = p.adjust(p_val.CEC, method='fdr', n = nrow(.))) %>%
  mutate(MSK.adj.p.value = p.adjust(p_val.MSK, method='fdr', n = nrow(.)))

write.csv(results.CEC.MSK, "./Ethanol Experiment/RM.withinpop-lateral.csv")


## Between-population modeling per day 
array.lateral <- arrayspecs(lateral.df[,4:37], 17, 2)
gdf <- geomorph.data.frame(coords = array.lateral, 
                           pop = dorsal.df$pop, 
                           day = dorsal.df$day)

stats <- as_tibble(day.vec) %>% 
  mutate(pval = NA,
         r2 = NA,
         F = NA)
i = 0

for(z in day.vec) {
  i = i+1
  coords.test <- gdf$coords[,,gdf$day == z]
  pop.test <- as.factor(c(rep("CEC", 30), rep("MSK", 30)))
  gdf.test <- geomorph.data.frame(coords = coords.test, pop = pop.test)
  
  fit <- procD.lm(coords ~ pop,
                  data = gdf.test)
  
  stats$pval[i] <- fit$aov.table$`Pr(>F)`[1]
  stats$r2[i] <- fit$aov.table$Rsq[1]
  stats$F[i] <- fit$aov.table$F[1]
  
}

stats <- stats %>% 
  mutate(qval = p.adjust(pval, method='fdr', n = nrow(.)))

## Between-population modeling per day 
#### Question: How would a Procrustes ANOVA that compares that populations change across time?
array.lateral <- arrayspecs(lateral.df[,4:37], 17, 2)
gdf <- geomorph.data.frame(coords = array.lateral, 
                           pop = lateral.df$pop, 
                           day = lateral.df$day)

varPairs <- combn(unique(lateral.df$day), 2, simplify = TRUE) %>%
  as.character() %>% as.data.frame() %>%
  group_by(grp = str_c('Column', rep(1:2, length.out = n()))) %>%
  mutate(rn = row_number()) %>%
  rename(value = 1) %>% 
  ungroup %>%
  pivot_wider(names_from = grp, values_from = value) %>%
  select(-rn)

varPairs.complete <- unique(lateral.df$day) %>%
  as_tibble %>% 
  rename(Column1 = value) %>% 
  mutate(Column2 = unique(lateral.df$day)) %>% 
  rbind(varPairs)

stats <- as_tibble(varPairs.complete) %>% 
  mutate(pval = NA,
         r2 = NA,
         Z = NA)
i = 0

## For all day-population combinations
for(z in (1:nrow(varPairs.complete))) {
  i = i+1
  coords1 <- gdf$coords[,,gdf$day == varPairs.complete$Column1[z]][,,1:30] ### CEC is column 1
  coords2 <- gdf$coords[,,gdf$day == varPairs.complete$Column2[z]][,,31:60] ### MSK is column 2
  coords.both <- abind(coords1, coords2, along = 3)
  
  pop.test <- as.factor(c(rep("CEC", 30), rep("MSK", 30)))
  gdf.test <- geomorph.data.frame(coords = coords.both, pop = pop.test)
  
  fit <- procD.lm(coords ~ pop,
                  data = gdf.test)
  
  stats$pval[i] <- fit$aov.table$`Pr(>F)`[1]
  stats$r2[i] <- fit$aov.table$Rsq[1]
  stats$Z[i] <- fit$aov.table$Z[1]
}

stats <- stats %>% 
  mutate(qval = p.adjust(pval, method='fdr', n = nrow(.)))


### Full matrix
stats1 <- stats %>% select(Column1, Column2, Z)
Zscores <- stats1 %>% pivot_wider(names_from = Column1, values_from = Z)