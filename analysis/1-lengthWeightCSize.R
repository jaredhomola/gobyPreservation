#################################################################
######       Ethanol experiment morphometric analyses      ######
######################     June  2020     #######################
#################################################################

###Load data packages###
library(nlme)
library(lmtest)
library(geomorph)
library(segmented)
library(gobyPreservation)
library(tidyverse)
data("gobyPreservation")

#establish factors
pop = as.factor(c(rep(c("CEC"), 1620), rep(c("MSK"), 1620)))
ind = as.factor(c(rep(1:30, each = 18, times = 3), rep(1:30, each = 18, times = 3)))
rep = as.factor(c(rep(1:3, each = 540), rep(1:3, each = 540)))
day = as.factor(c(rep(rep(c(0, 1, 2, 7, 14, 21, 28, 35, 42, 49, 56, 63, 70, 77, 84, 105, 126, 154), times = 60), times = 3)))



### Set up dataframe, including adding CSize
## Dorsal CSize calculation
proc.dorsal <- gpagen(dat.dorsal)
# Take mean of 3 reps
cSizes.dorsal <- data.frame(keyName=names(proc.dorsal$Csize), value=proc.dorsal$Csize, row.names=NULL) %>% 
  cbind(pop) %>% cbind(ind) %>% cbind(rep) %>% cbind(day) %>% 
  as_tibble() %>% 
  mutate(dayMeasure = day,
         ind = as.numeric(ind)) %>% 
  group_by(pop, ind, dayMeasure) %>% 
  summarize(meanCsizeDorsal = mean(value))

## Side CSize calculation
proc.side <- gpagen(dat.lateral)

# Take mean of 3 reps
cSizes.side <- data.frame(keyName=names(proc.side$Csize), value=proc.side$Csize, row.names=NULL) %>% 
  cbind(pop) %>% cbind(ind) %>% cbind(rep) %>% cbind(day) %>% 
  as_tibble() %>% 
  mutate(dayMeasure = day,
         ind = as.numeric(ind)) %>% 
  group_by(pop, ind, dayMeasure) %>% 
  summarize(meanCsizeSide = mean(value))


### Assembly overall dataframe
dat <- gobyData %>% 
  separate(sampleName,c("time", "pop", "ind")) %>% 
  as_tibble() %>% 
  mutate(pop = recode(pop, "MSK02" = "MSK", "CEC01" = "CEC")) %>% 
  mutate(pop = as.factor(pop),
         ind = as.numeric(ind),
         dayMeasure = as.factor(as.character(dayMeasure)),
         stdlen = forkLength) %>% 
  left_join(cSizes.dorsal) %>% 
  left_join(cSizes.side) %>% 
  mutate(dayMeasure = as.numeric(dayMeasure)) %>% 
  mutate(ind = ifelse(pop == "MSK", 
                      ind + 30, 
                      ind))

### Create Regression plot ####
### Mass
ggplot(dat, aes(dayMeasure, mass, group = pop)) +
  geom_point(aes(shape = pop), alpha = 0.2) +
  geom_smooth(method = "auto", aes(lty = pop)) + #The "span" variable will change the way smoothing is performed.
  ylab("Mass (g)") + #"ylab" lets you rename the y axis label
  xlab("Days post-preservation") +
  scale_color_manual(labels = c("Cedar Creek", "Lake Michigan"), values = c("blue", "red")) + #Change legend names
  guides(color = guide_legend(override.aes = list(fill=NA))) + #Remove gray background from legend symbols
  theme_classic(base_size = 18) + #Changes overall font sizes
  theme(legend.title = element_blank(),  #Removes legend title
        legend.text=element_text(size=18),
        legend.position = "bottom")  #Increases font size for legend


### Standard length
ggplot(dat, aes(dayMeasure, forkLength, group = pop)) +
  geom_point(aes(color = pop), alpha = 0.2) +
  geom_smooth(method = "auto" , aes(color = pop)) +
  xlab("Days post-preservation") +
  ylab("Standard length") +
  scale_color_manual(labels = c("Cedar Creek", "Lake Michigan"), values = c("blue", "red")) +
  guides(color = guide_legend(override.aes = list(fill=NA))) + #Remove gray background from legend symbols
  theme_classic(base_size = 18) + #Changes overall font sizes
  theme(legend.title = element_blank(),  #Removes legend title
        legend.text=element_text(size=18),
        legend.position = "bottom")

### Mass
ggplot(dat, aes(dayMeasure, mass, group = pop)) +
  geom_point(aes(color = pop), alpha = 0.2) +
  geom_smooth(method = "auto", aes(color = pop)) + #The "span" variable will change the way smoothing is performed.
  ylab("Mass (g)") + #"ylab" lets you rename the y axis label
  xlab("Days post-preservation") +
  scale_color_manual(labels = c("Cedar Creek", "Lake Michigan"), values = c("blue", "red")) + #Change legend names
  guides(color = guide_legend(override.aes = list(fill=NA))) + #Remove gray background from legend symbols
  theme_classic(base_size = 18) + #Changes overall font sizes
  theme(legend.title = element_blank(),  #Removes legend title
        legend.text=element_text(size=18),
        legend.position = "bottom")  #Increases font size for legend

#### Generalized least squares modeling to account for autocorrelation structure ####
#### GLS model for mass ####
mass.mod.11 <-
  gls(log(mass) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 1, q = 1),
    method = "ML")

mass.mod.12 <-
  gls(log(mass) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 1, q = 2),
    method = "ML")

mass.mod.13 <-
  gls(log(mass) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 1, q = 3),
    method = "ML")

mass.mod.14 <-
  gls(log(mass) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 1, q = 4),
    method = "ML")

mass.mod.15 <-
  gls(log(mass) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 1, q = 5),
    method = "ML")

mass.mod.21 <-
  gls(log(mass) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 2, q = 1),
    method = "ML")

mass.mod.22 <-
  gls(log(mass) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 2, q = 2),
    method = "ML")

mass.mod.23 <-
  gls(log(mass) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 2, q = 3),
    method = "ML")

mass.mod.24 <-
  gls(log(mass) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 2, q = 4),
    method = "ML")

mass.mod.25 <-
  gls(log(mass) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 2, q = 5),
    method = "ML")

mass.mod.31 <-
  gls(log(mass) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 3, q = 1),
    method = "ML")

mass.mod.32 <-
  gls(log(mass) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 3, q = 2),
    method = "ML")

mass.mod.33 <-
  gls(log(mass) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 3, q = 3),
    method = "ML")

mass.mod.34 <-
  gls(log(mass) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 3, q = 4),
    method = "ML")

mass.mod.35 <-
  gls(log(mass) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 3, q = 5),
    method = "ML")

mass.mod.41 <-
  gls(log(mass) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 4, q = 1),
    method = "ML")

mass.mod.42 <-
  gls(log(mass) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 4, q = 2),
    method = "ML")

mass.mod.43 <-
  gls(log(mass) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 4, q = 3),
    method = "ML")

mass.mod.44 <-
  gls(log(mass) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 4, q = 4),
    method = "ML")

mass.mod.45 <-
  gls(log(mass) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 4, q = 5),
    method = "ML")

mass.mod.51 <-
  gls(log(mass) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 5, q = 1),
    method = "ML")

mass.mod.52 <-
  gls(log(mass) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 5, q = 2),
    method = "ML")

mass.mod.53 <-
  gls(log(mass) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 5, q = 3),
    method = "ML")

mass.mod.54 <-
  gls(log(mass) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 5, q = 4),
    method = "ML")

mass.mod.55 <-
  gls(log(mass) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 5, q = 5),
    method = "ML")


## Test other error autocorrelation models
anova(mass.mod.11,
      mass.mod.12,
      mass.mod.13,
      mass.mod.14,
      mass.mod.15,
      mass.mod.21,
      mass.mod.22,
      mass.mod.23,
      #mass.mod.24, Did not converge
      mass.mod.25,
      #mass.mod.31, Did not converge
      mass.mod.32,
      mass.mod.33,
      #mass.mod.34, Did not converge
      mass.mod.35,
      mass.mod.41,
      mass.mod.42,
      mass.mod.43,
      mass.mod.44,
      #mass.mod.45, Did not converge
      mass.mod.51,
      mass.mod.52,
      mass.mod.53,  ## Most supported model
      mass.mod.54
      #mass.mod.55 Did not converge
      )

summary(mass.mod.53)


#### GLS model for length ####
stdlen.mod.11 <-
  gls(log(stdlen) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 1, q = 1),
    method = "ML")

stdlen.mod.12 <-
  gls(log(stdlen) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 1, q = 2),
    method = "ML")

stdlen.mod.13 <-
  gls(log(stdlen) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 1, q = 3),
    method = "ML")

stdlen.mod.14 <-
  gls(log(stdlen) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 1, q = 4),
    method = "ML")

stdlen.mod.15 <-
  gls(log(stdlen) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 1, q = 5),
    method = "ML")

stdlen.mod.21 <-
  gls(log(stdlen) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 2, q = 1),
    method = "ML")

stdlen.mod.22 <-
  gls(log(stdlen) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 2, q = 2),
    method = "ML")

stdlen.mod.23 <-
  gls(log(stdlen) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 2, q = 3),
    method = "ML")

stdlen.mod.24 <-
  gls(log(stdlen) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 2, q = 4),
    method = "ML")

stdlen.mod.25 <-
  gls(log(stdlen) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 2, q = 5),
    method = "ML")

stdlen.mod.31 <-
  gls(log(stdlen) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 3, q = 1),
    method = "ML")

stdlen.mod.32 <-
  gls(log(stdlen) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 3, q = 2),
    method = "ML")

stdlen.mod.33 <-
  gls(log(stdlen) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 3, q = 3),
    method = "ML")

stdlen.mod.34 <-
  gls(log(stdlen) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 3, q = 4),
    method = "ML")

stdlen.mod.35 <-
  gls(log(stdlen) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 3, q = 5),
    method = "ML")

stdlen.mod.41 <-
  gls(log(stdlen) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 4, q = 1),
    method = "ML")

stdlen.mod.42 <-
  gls(log(stdlen) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 4, q = 2),
    method = "ML")

stdlen.mod.43 <-
  gls(log(stdlen) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 4, q = 3),
    method = "ML")

stdlen.mod.44 <-
  gls(log(stdlen) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 4, q = 4),
    method = "ML")

stdlen.mod.45 <-
  gls(log(stdlen) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 4, q = 5),
    method = "ML")

stdlen.mod.51 <-
  gls(log(stdlen) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 5, q = 1),
    method = "ML")

stdlen.mod.52 <-
  gls(log(stdlen) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 5, q = 2),
    method = "ML")

stdlen.mod.53 <-
  gls(log(stdlen) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 5, q = 3),
    method = "ML")

stdlen.mod.54 <-
  gls(log(stdlen) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 5, q = 4),
    method = "ML")

stdlen.mod.55 <-
  gls(log(stdlen) ~ dayMeasure + pop + dayMeasure*pop,
    data = dat,
    correlation = corARMA(form = ~1|ind, p = 5, q = 5),
    method = "ML")


## Test other error autocorrelation models
anova(stdlen.mod.11,
      stdlen.mod.12,
      stdlen.mod.13, ## Most supported model
      #stdlen.mod.14, Did not converge
      #stdlen.mod.15, Did not converge
      stdlen.mod.21, 
      stdlen.mod.22,
      stdlen.mod.23,
      #stdlen.mod.24, Did not converge
      stdlen.mod.25,
      stdlen.mod.31,
      stdlen.mod.32,
      stdlen.mod.33,
      #stdlen.mod.34, Did not converge
      stdlen.mod.35,
      stdlen.mod.41,
      stdlen.mod.42,
      stdlen.mod.43,
      stdlen.mod.44,
      #stdlen.mod.45, Did not converge
      stdlen.mod.51,
      stdlen.mod.52,
      stdlen.mod.53,  
      stdlen.mod.54,
      stdlen.mod.55)

summary(stdlen.mod.13)


#### GLS model for meanCsizeDorsal ####
meanCsizeDorsal.mod.11 <-
  gls(log(meanCsizeDorsal) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 1, q = 1),
      method = "ML")

meanCsizeDorsal.mod.12 <-
  gls(log(meanCsizeDorsal) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 1, q = 2),
      method = "ML")

meanCsizeDorsal.mod.13 <-
  gls(log(meanCsizeDorsal) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 1, q = 3),
      method = "ML")

meanCsizeDorsal.mod.14 <-
  gls(log(meanCsizeDorsal) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 1, q = 4),
      method = "ML")

meanCsizeDorsal.mod.15 <-
  gls(log(meanCsizeDorsal) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 1, q = 5),
      method = "ML")

meanCsizeDorsal.mod.21 <-
  gls(log(meanCsizeDorsal) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 2, q = 1),
      method = "ML")

meanCsizeDorsal.mod.22 <-
  gls(log(meanCsizeDorsal) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 2, q = 2),
      method = "ML")

meanCsizeDorsal.mod.23 <-
  gls(log(meanCsizeDorsal) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 2, q = 3),
      method = "ML")

meanCsizeDorsal.mod.24 <-
  gls(log(meanCsizeDorsal) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 2, q = 4),
      method = "ML")

meanCsizeDorsal.mod.25 <-
  gls(log(meanCsizeDorsal) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 2, q = 5),
      method = "ML")

meanCsizeDorsal.mod.31 <-
  gls(log(meanCsizeDorsal) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 3, q = 1),
      method = "ML")

meanCsizeDorsal.mod.32 <-
  gls(log(meanCsizeDorsal) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 3, q = 2),
      method = "ML")

meanCsizeDorsal.mod.33 <-
  gls(log(meanCsizeDorsal) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 3, q = 3),
      method = "ML")

meanCsizeDorsal.mod.34 <-
  gls(log(meanCsizeDorsal) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 3, q = 4),
      method = "ML")

meanCsizeDorsal.mod.35 <-
  gls(log(meanCsizeDorsal) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 3, q = 5),
      method = "ML")

meanCsizeDorsal.mod.41 <-
  gls(log(meanCsizeDorsal) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 4, q = 1),
      method = "ML")

meanCsizeDorsal.mod.42 <-
  gls(log(meanCsizeDorsal) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 4, q = 2),
      method = "ML")

meanCsizeDorsal.mod.43 <-
  gls(log(meanCsizeDorsal) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 4, q = 3),
      method = "ML")

meanCsizeDorsal.mod.44 <-
  gls(log(meanCsizeDorsal) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 4, q = 4),
      method = "ML")

meanCsizeDorsal.mod.45 <-
  gls(log(meanCsizeDorsal) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 4, q = 5),
      method = "ML")

meanCsizeDorsal.mod.51 <-
  gls(log(meanCsizeDorsal) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 5, q = 1),
      method = "ML")

meanCsizeDorsal.mod.52 <-
  gls(log(meanCsizeDorsal) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 5, q = 2),
      method = "ML")

meanCsizeDorsal.mod.53 <-
  gls(log(meanCsizeDorsal) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 5, q = 3),
      method = "ML")

meanCsizeDorsal.mod.54 <-
  gls(log(meanCsizeDorsal) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 5, q = 4),
      method = "ML")

meanCsizeDorsal.mod.55 <-
  gls(log(meanCsizeDorsal) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 5, q = 5),
      method = "ML")


## Test other error autocorrelation models
anova(meanCsizeDorsal.mod.11,
      meanCsizeDorsal.mod.12,
      #meanCsizeDorsal.mod.13, Did not converge
      meanCsizeDorsal.mod.14,
      meanCsizeDorsal.mod.15,
      #meanCsizeDorsal.mod.21, Did not converge
      meanCsizeDorsal.mod.22,
      meanCsizeDorsal.mod.23,
      meanCsizeDorsal.mod.24,
      meanCsizeDorsal.mod.25,
      #meanCsizeDorsal.mod.31, Did not converge
      #meanCsizeDorsal.mod.32, Did not converge
      meanCsizeDorsal.mod.33,
      #meanCsizeDorsal.mod.34, Did not converge
      #meanCsizeDorsal.mod.35, Did not converge
      meanCsizeDorsal.mod.41,
      #meanCsizeDorsal.mod.42,Did not converge
      meanCsizeDorsal.mod.43,
      meanCsizeDorsal.mod.44,
      #meanCsizeDorsal.mod.45, Did not converge
      meanCsizeDorsal.mod.51,
      meanCsizeDorsal.mod.52,  ## Most supported model
      meanCsizeDorsal.mod.53,
      meanCsizeDorsal.mod.54
     #meanCsizeDorsal.mod.55
)

summary(meanCsizeDorsal.mod.52)





#### GLS model for meanCsizeSide ####
meanCsizeSide.mod.11 <-
  gls(log(meanCsizeSide) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 1, q = 1),
      method = "ML")

meanCsizeSide.mod.12 <-
  gls(log(meanCsizeSide) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 1, q = 2),
      method = "ML")

meanCsizeSide.mod.13 <-
  gls(log(meanCsizeSide) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 1, q = 3),
      method = "ML")

meanCsizeSide.mod.14 <-
  gls(log(meanCsizeSide) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 1, q = 4),
      method = "ML")

meanCsizeSide.mod.15 <-
  gls(log(meanCsizeSide) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 1, q = 5),
      method = "ML")

meanCsizeSide.mod.21 <-
  gls(log(meanCsizeSide) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 2, q = 1),
      method = "ML")

meanCsizeSide.mod.22 <-
  gls(log(meanCsizeSide) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 2, q = 2),
      method = "ML")

meanCsizeSide.mod.23 <-
  gls(log(meanCsizeSide) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 2, q = 3),
      method = "ML")

meanCsizeSide.mod.24 <-
  gls(log(meanCsizeSide) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 2, q = 4),
      method = "ML")

meanCsizeSide.mod.25 <-
  gls(log(meanCsizeSide) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 2, q = 5),
      method = "ML")

meanCsizeSide.mod.31 <-
  gls(log(meanCsizeSide) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 3, q = 1),
      method = "ML")

meanCsizeSide.mod.32 <-
  gls(log(meanCsizeSide) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 3, q = 2),
      method = "ML")

meanCsizeSide.mod.33 <-
  gls(log(meanCsizeSide) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 3, q = 3),
      method = "ML")

meanCsizeSide.mod.34 <-
  gls(log(meanCsizeSide) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 3, q = 4),
      method = "ML")

meanCsizeSide.mod.35 <-
  gls(log(meanCsizeSide) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 3, q = 5),
      method = "ML")

meanCsizeSide.mod.41 <-
  gls(log(meanCsizeSide) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 4, q = 1),
      method = "ML")

meanCsizeSide.mod.42 <-
  gls(log(meanCsizeSide) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 4, q = 2),
      method = "ML")

meanCsizeSide.mod.43 <-
  gls(log(meanCsizeSide) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 4, q = 3),
      method = "ML")

meanCsizeSide.mod.44 <-
  gls(log(meanCsizeSide) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 4, q = 4),
      method = "ML")

meanCsizeSide.mod.45 <-
  gls(log(meanCsizeSide) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 4, q = 5),
      method = "ML")

meanCsizeSide.mod.51 <-
  gls(log(meanCsizeSide) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 5, q = 1),
      method = "ML")

meanCsizeSide.mod.52 <-
  gls(log(meanCsizeSide) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 5, q = 2),
      method = "ML")

meanCsizeSide.mod.53 <-
  gls(log(meanCsizeSide) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 5, q = 3),
      method = "ML")

meanCsizeSide.mod.54 <-
  gls(log(meanCsizeSide) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 5, q = 4),
      method = "ML")

meanCsizeSide.mod.55 <-
  gls(log(meanCsizeSide) ~ dayMeasure + pop + dayMeasure*pop,
      data = dat,
      correlation = corARMA(form = ~1|ind, p = 5, q = 5),
      method = "ML")


## Test other error autocorrelation models
anova(meanCsizeSide.mod.11,
      meanCsizeSide.mod.12,
      meanCsizeSide.mod.13,
      meanCsizeSide.mod.14,
      meanCsizeSide.mod.15,
      meanCsizeSide.mod.21,
      meanCsizeSide.mod.22,
      meanCsizeSide.mod.23,
      meanCsizeSide.mod.24,  ## Most supported model
      meanCsizeSide.mod.25,
      meanCsizeSide.mod.31,
      meanCsizeSide.mod.32,
      meanCsizeSide.mod.33,
      #meanCsizeSide.mod.34,
      #meanCsizeSide.mod.35,
      meanCsizeSide.mod.41,
      meanCsizeSide.mod.42,
      meanCsizeSide.mod.43,
      meanCsizeSide.mod.44,
      #meanCsizeSide.mod.45,
      meanCsizeSide.mod.51,
      meanCsizeSide.mod.52, 
      meanCsizeSide.mod.53,
      meanCsizeSide.mod.54
      #meanCsizeSide.mod.55
)

summary(meanCsizeSide.mod.24)

#### Segmented regression ####
## Mass
mass.lm <- lm(log(mass) ~ dayMeasure + pop + dayMeasure*pop, data = dat) 
mass.seg <- segmented(mass.lm, 
                    seg.Z = ~ dayMeasure, 
                    npsi = 1)
summary(mass.seg)

# Identify breakpoints
mass.seg$psi

## Standard length
sl.lm <- lm(log(stdlen) ~ dayMeasure + pop + dayMeasure*pop, data = dat) 
sl.seg <- segmented(sl.lm, 
                      seg.Z = ~ dayMeasure, 
                      npsi = 1)
summary(sl.seg)
# Identify breakpoints
sl.seg$psi

## Dorsal centroid size
meanCsizedorsal.lm <- lm(log(meanCsizeDorsal) ~ dayMeasure + pop + dayMeasure*pop, data = dat) 
meanCsizedorsal.seg <- segmented(meanCsizedorsal.lm, 
                    seg.Z = ~ dayMeasure, 
                    npsi = 1)
summary(meanCsizedorsal.seg)
# Identify breakpoints
meanCsizedorsal.seg$psi


## Side centroid size
meanCsizeSide.lm <- lm(log(meanCsizeSide) ~ dayMeasure + pop + dayMeasure*pop, data = dat) 
meanCsizeSide.seg <- segmented(meanCsizeSide.lm, 
                                 seg.Z = ~ dayMeasure, 
                                 npsi = 1)
summary(meanCsizeSide.seg)
# Identify breakpoints
meanCsizeSide.seg$psi


