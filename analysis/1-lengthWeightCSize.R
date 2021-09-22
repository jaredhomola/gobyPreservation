#################################################################
######       Ethanol experiment morphometric analyses      ######
######################     Sept  2021     #######################
#################################################################

###Load data packages###
library(nlme)
library(lmtest)
library(geomorph)
library(segmented)
library(car)
library(gobyPreservation)
library(tidyverse)
data("gobyPreservation")

### Test for normality
# Raw data
shapiro.test(gobyData$mass)
shapiro.test(gobyData$stdlen)
shapiro.test(gobyData$meanCsizeDorsal)
shapiro.test(gobyData$meanCsizeSide)

# Transformed data
shapiro.test(log(gobyData$mass))
shapiro.test(log(gobyData$stdlen))
shapiro.test(log(gobyData$meanCsizeDorsal))
shapiro.test(log(gobyData$meanCsizeSide))

### Plot the data
datPlot <- gobyData %>% 
  select(pop, ind, dayMeasure, mass, stdlen, meanCsizeDorsal, meanCsizeSide) %>% 
  pivot_longer(cols = mass:meanCsizeSide,
               names_to = "variable",
               values_to = "value") %>% 
  mutate(variable = recode(variable,
                           mass = "Mass",
                           stdlen = "Standard_Length",
                           meanCsizeDorsal = "Centroid_Size_-_Dorsal",
                           meanCsizeSide = "Centroid_Size_-_Lateral"))

datPlot2 <- datPlot %>%
  group_by(variable, ind) %>%
  mutate(difference = value - lag(value, default = 0)) %>% 
  filter(dayMeasure != 0) %>% 
  ungroup()

datPlot2 %>% 
  ggplot(aes(x = as.factor(dayMeasure), y = difference, fill = factor(pop))) +
  geom_boxplot(outlier.alpha = 0.2) +
  scale_fill_grey(start = 0.5, end = 1) +
  xlab("Days post-preservation") +
  ylab("Difference from previous time point") +
  facet_grid(rows = vars(variable),
             scales = "free_y") +
  theme_bw(base_size = 18) + 
  theme(legend.title = element_blank(),
        legend.text=element_text(size=18),
        legend.position = "none")


#### Generalized least squares modeling to account for autocorrelation structure ####
#### GLS model for mass ####
modVar <- crossing(var1 = 0:5, var2 = 0:5) %>% 
  rename(p = 1, 
         q = 2)

i = 1
mass.res.df <- data.frame(p = numeric(), 
                     q = numeric(), 
                     aic = numeric())

for (i in 1:nrow(modVar)) {
  tryCatch({
    mod <- gls(log(mass) ~ dayMeasure + pop + dayMeasure*pop,
               data = gobyData,
               correlation = corARMA(form = ~1|ind, 
                                     p = as.numeric(modVar[i, 1]), 
                                     q = as.numeric(modVar[i, 2])),
               method = "ML")
    
    mass.res.df[i,1] <- as.numeric(modVar[i, 1])
    mass.res.df[i,2] <- as.numeric(modVar[i, 2])
    mass.res.df[i,3] <- AIC(mod)
  }, error=function(e){})
  
  i = i+1
}

mass.res.df %>% arrange(aic) %>% head()

mass.mod <-
  gls(log(mass) ~ dayMeasure + pop + dayMeasure*pop,
      data = gobyData,
      correlation = corARMA(form = ~1|ind, p = 4, q = 3),
      method = "ML")

summary(mass.mod)
Anova(mass.mod)

acf(residuals(mass.mod))


#### GLS model for length ####
i = 1
stdlen.res.df <- data.frame(p = numeric(), 
                     q = numeric(), 
                     aic = numeric())

for (i in 1:nrow(modVar)) {
  
  tryCatch({
    mod <- gls(log(stdlen) ~ dayMeasure + pop + dayMeasure*pop,
               data = gobyData,
               correlation = corARMA(form = ~1|ind, 
                                     p = as.numeric(modVar[i, 1]), 
                                     q = as.numeric(modVar[i, 2])),
               method = "ML")
    
    stdlen.res.df[i,1] <- as.numeric(modVar[i, 1])
    stdlen.res.df[i,2] <- as.numeric(modVar[i, 2])
    stdlen.res.df[i,3] <- AIC(mod)
  }, error=function(e){})
  
  i = i+1
}

stdlen.res.df %>% arrange(aic) %>% head()

stdlen.mod <-
  gls(log(stdlen) ~ dayMeasure + pop + dayMeasure*pop,
      data = gobyData,
      correlation = corARMA(form = ~1|ind, p = 2, q = 1),
      method = "ML")

summary(stdlen.mod)
Anova(stdlen.mod)

acf(residuals(stdlen.mod))


#### GLS model for meanCsizeDorsal ####
i = 1
meanCsizeDorsal.res.df <- data.frame(p = numeric(), 
                     q = numeric(), 
                     aic = numeric())

for (i in 1:nrow(modVar)) {
  
  tryCatch({
    mod <- gls(log(meanCsizeDorsal) ~ dayMeasure + pop + dayMeasure*pop,
               data = gobyData,
               correlation = corARMA(form = ~1|ind, 
                                     p = as.numeric(modVar[i, 1]), 
                                     q = as.numeric(modVar[i, 2])),
               method = "ML")
    
    meanCsizeDorsal.res.df[i,1] <- as.numeric(modVar[i, 1])
    meanCsizeDorsal.res.df[i,2] <- as.numeric(modVar[i, 2])
    meanCsizeDorsal.res.df[i,3] <- AIC(mod)
  }, error=function(e){})
  
  i = i+1
}

meanCsizeDorsal.res.df %>% arrange(aic) %>% head()

meanCsizeDorsal.mod <-
  gls(log(meanCsizeDorsal) ~ dayMeasure + pop + dayMeasure*pop,
      data = gobyData,
      correlation = corARMA(form = ~1|ind, p = 2, q = 1),
      method = "ML")

summary(meanCsizeDorsal.mod)
Anova(meanCsizeDorsal.mod)

acf(residuals(meanCsizeDorsal.mod))


#### GLS model for meanCsizeSide ####
i = 1
meanCsizeSide.res.df <- data.frame(p = numeric(), 
                     q = numeric(), 
                     aic = numeric())

for (i in 1:nrow(modVar)) {
  
  tryCatch({
    mod <- gls(log(meanCsizeSide) ~ dayMeasure + pop + dayMeasure*pop,
               data = gobyData,
               correlation = corARMA(form = ~1|ind, 
                                     p = as.numeric(modVar[i, 1]), 
                                     q = as.numeric(modVar[i, 2])),
               method = "ML")
    
    meanCsizeSide.res.df[i,1] <- as.numeric(modVar[i, 1])
    meanCsizeSide.res.df[i,2] <- as.numeric(modVar[i, 2])
    meanCsizeSide.res.df[i,3] <- AIC(mod)
  }, error=function(e){})
  
  i = i+1
}

meanCsizeSide.res.df %>% arrange(aic) %>% head()

meanCsizeSide.mod <-
  gls(log(meanCsizeSide) ~ dayMeasure + pop + dayMeasure*pop,
      data = gobyData,
      correlation = corARMA(form = ~1|ind, p = 2, q = 4),
      method = "ML")

summary(meanCsizeSide.mod)
Anova(meanCsizeSide.mod)

acf(residuals(meanCsizeSide.mod))

#### Segmented regression ####
## Mass
mass.lm <- lm(log(mass) ~ dayMeasure + pop + dayMeasure*pop, data = gobyData) 
mass.seg <- segmented(mass.lm, 
                    seg.Z = ~ dayMeasure, 
                    npsi = 1)
summary(mass.seg)

# Identify breakpoints
mass.seg$psi

## Standard length
sl.lm <- lm(log(stdlen) ~ dayMeasure + pop + dayMeasure*pop, data = gobyData) 
sl.seg <- segmented(sl.lm, 
                      seg.Z = ~ dayMeasure, 
                      npsi = 1)
summary(sl.seg)
# Identify breakpoints
sl.seg$psi

## Dorsal centroid size
meanCsizedorsal.lm <- lm(log(meanCsizeDorsal) ~ dayMeasure + pop + dayMeasure*pop, data = gobyData) 
meanCsizedorsal.seg <- segmented(meanCsizedorsal.lm, 
                    seg.Z = ~ dayMeasure, 
                    npsi = 1)
summary(meanCsizedorsal.seg)
# Identify breakpoints
meanCsizedorsal.seg$psi


## Side centroid size
meanCsizeSide.lm <- lm(log(meanCsizeSide) ~ dayMeasure + pop + dayMeasure*pop, data = gobyData) 
meanCsizeSide.seg <- segmented(meanCsizeSide.lm, 
                                 seg.Z = ~ dayMeasure, 
                                 npsi = 1)
summary(meanCsizeSide.seg)
# Identify breakpoints
meanCsizeSide.seg$psi


