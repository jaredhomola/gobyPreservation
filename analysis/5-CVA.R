#################################################################
######       Ethanol experiment morphometric analyses      ######
######################     June  2020     #######################
##################################################################Bailey's setwd


library(geomorph)
library(Morpho)
library(shapes)
library(tidyverse)

## Establish factors
pop = as.factor(c(rep(c("CEC"), 1620), rep(c("MSK"), 1620)))
ind = as.factor(c(rep(1:30, each = 18, times = 3), rep(1:30, each = 18, times = 3)))
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

dorsal.df <- data.frame(two.d.array(adj.shape)) %>% 
  mutate(pop = pop, 
         ind = ind, 
         rep = rep, 
         day = day) %>% 
  group_by(pop, ind, day) %>% 
  summarize_if(is.numeric, funs(mean))

array.dorsal <- arrayspecs(dorsal.df[,4:15], 6, 2)
pop.daySub <- as.factor(c(rep("CEC", 30), rep("MSK", 30)))

## Function to perform test ##
gobyCVA <- function(dayTarget, ...) {
  dat.day <- array.dorsal[,,dorsal.df$day == dayTarget]
  proc.day <- procSym(dat.day)
  CVA(proc.day$orpdata, pop.daySub, rounds = 10000, cv=TRUE)
}

## Loop through applying function to each day in preservation time series ##
results.dorsal <- tibble(1:length(unique(day.vec))) %>% 
  rename(reassignmentRate.dorsal = `1:length(unique(day.vec))`) %>% 
  mutate(day = as.numeric(as.character(unique(day))))
i = 0

for(z in unique(day.vec)) {
  i = i+1
  cva.out <- gobyCVA(z)
  reassignment  <- cva.out$posterior %>% 
    as_tibble(rownames = NA) %>% 
    mutate(assignment = if_else(CEC > 0.5, "CEC", "MSK"),
           assignment = as.factor(assignment),
           actual = pop.daySub) %>% 
    mutate(match = if_else(assignment == actual, "yes", "no")) %>% 
    summarize(reassignmentRate.dorsal = sum(match == "yes") / 60)
  
  results.dorsal$reassignmentRate.dorsal[i] <- reassignment$reassignmentRate.dorsal
}




### Lateral ####

dat.side <-
  readland.tps("./Ethanol Experiment/Experiment Photos/Experiment-photos-side-UNBEND/unbentAll.TPS",
               specID = "ID",
               readcurves = FALSE)

dat.side <- dat.side[c(1:17),,]

## Perform Generalized Procrustes transformation
proc.side <- gpagen(dat.side)

gdf.side <- geomorph.data.frame(proc.side,
                                  shape = proc.side$coords) 

fit.simple.allo <- procD.lm(coords ~ log(Csize), 
                            data = gdf.side, print.progress = FALSE)
shape.resid <- arrayspecs(fit.simple.allo$residuals,
                          p=dim(gdf.side$coords)[1], k=dim(gdf.side$coords)[2]) # allometry-adjusted residuals
adj.shape <- shape.resid + array(mshape(gdf.side$coords), dim(shape.resid)) # Size adjusted shapes

side.df <- data.frame(two.d.array(adj.shape)) %>% 
  mutate(pop = pop, 
         ind = ind, 
         rep = rep, 
         day = day) %>% 
  group_by(pop, ind, day) %>% 
  summarize_if(is.numeric, funs(mean))

array.side <- arrayspecs(side.df[,4:37], 17, 2)
pop.daySub <- as.factor(c(rep("CEC", 30), rep("MSK", 30)))

## Function to perform test ##
gobyCVA <- function(dayTarget, ...) {
  dat.day <- array.side[,,side.df$day == dayTarget]
  proc.day <- procSym(dat.day)
  CVA(proc.day$orpdata, pop.daySub, rounds = 10000, cv=TRUE)
}

## Loop through applying function to each day in preservation time series ##
results <- tibble(1:length(unique(day.vec))) %>% 
  rename(reassignmentRate = `1:length(unique(day.vec))`) %>% 
  mutate(day = as.numeric(as.character(unique(day))))
i = 0

for(z in unique(day.vec)) {
  i = i+1
  cva.out <- gobyCVA(z)
  reassignment  <- cva.out$posterior %>% 
    as_tibble(rownames = NA) %>% 
    mutate(assignment = if_else(CEC > 0.5, "CEC", "MSK"),
           assignment = as.factor(assignment),
           actual = pop.daySub) %>% 
    mutate(match = if_else(assignment == actual, "yes", "no")) %>% 
    summarize(reassignmentRate = sum(match == "yes") / 60)
  
  results$reassignmentRate[i] <- reassignment$reassignmentRate
}

results.all <- left_join(results, results.dorsal) %>% 
  pivot_longer(cols = c("reassignmentRate", "reassignmentRate.dorsal"),
               names_to = "view",
               values_to = "reassignmentRate") %>% 
  mutate(view = recode(view,
                       reassignmentRate = "Lateral",
                       reassignmentRate.dorsal = "Dorsal"))

results.all %>% 
  ggplot(aes(x = day, y = reassignmentRate)) +
  geom_point() +
  geom_smooth(method = 'lm', color = "black") +
  xlab("Days post-preservation") +
  ylab("CVA reassignment rate") +
  facet_grid(rows = vars(view)) +
  theme_bw(base_size = 18)


ggsave("./Ethanol Experiment/plots/CVA.pdf",
           width = 9, height = 5, units = "in")

