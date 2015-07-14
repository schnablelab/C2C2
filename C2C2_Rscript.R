## Last Update: 7/14/2015

## ----------------------------------------------------------- Load required packages -------------------------------------------------- ##

library(magrittr)
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(mgcv)
library(pracma)
library(Hmisc)

## ---------------------------------------------------------- Read in required data sets ------------------------------------------------ ##

## read in five data sets:
##     r1: Carbon emission and greenhouse data generated from the first run (2/21/14 - 4/17/14).
##     r2: Carbon emission and greenhouse data generated from the second run (3/28/14 - 5/22/14).
##         Note that in both runs, measurements were taken every two hours for the first two days, and no more than once a day thereafter
##     GH_r1: greenhouse temperature (was measured every 15 minutes) corresponding to the first run (2/21/14 - 4/17/14).              
##     GH_r2: greenhouse temperature (was measured every 15 minutes) corresponding to the second run (2/21/14 - 4/17/14).              
##     LigninData: measures of lignin composition, lignin content and input carbon of each genotype in both runs.

r1 <- read.csv("Run1.csv", header = TRUE) 
r2 <- read.csv("Run2.csv", header = TRUE)
GH_r1 <- read.csv("C2C2_GH_data_Run1.csv", header = TRUE) 
GH_r2 <- read.csv("C2C2_GH_data_Run2.csv", header = TRUE) 
LigninData <- read.csv("LigninData_Update.csv", header = TRUE)


## ----------------------------------------------------------------- Data analysis ------------------------------------------------------- ##

## --------------------------------------------------------------------------------------------------------------------------------------- ##
## ------------------ Section: Estimation of total carbon emission from each mesocosm in Supplementary Materials.------------------------- ##
## --------------------------------------------------------------------------------------------------------------------------------------- ##


## For the four data sets, create sine and cosine terms for the variable "Day", which will be required in the model fitting.
r1 <- r1 %>% mutate(s = sin(2*pi*Day), c = cos(2*pi*Day))
r2 <- r2 %>% mutate(s = sin(2*pi*Day), c = cos(2*pi*Day))
GH_r1 <- GH_r1 %>% mutate(s = sin(2*pi*Day), c = cos(2*pi*Day))
GH_r2 <- GH_r2 %>% mutate(s = sin(2*pi*Day), c = cos(2*pi*Day))


run1_result <- NULL ## run1_results will store the predicted CO2 emission for the first run

for(i in 1:19){ ## go through genotype 1 to 19
  
  ## iGeno: carbon emission and greenhouse data for the ith genotype
  iGeno <-  r1 %>% filter(Genotype == Genotype[i]) 

  ## Step 1: model the relationship between CO2 emission and soil temperature using gam function (formula (1) in Supplementary Materials)
  Flux.gam <- gam(log(Exp_Flux) ~ s(Day, bs = "cr", k=18) + Soil_temp,  select=TRUE, data = iGeno)
  
  ## Step 2: model the relationship between soil temperature and air temperature using gam function (formula (2) in Supplementary Materials)
  Soil.gam = gam(Soil_temp~s+c+s(Day, bs="cr")+s(GH_temp_lag1, bs="cr"), select=TRUE, data = iGeno)
  
  ## Step 3: given estimates of Soil.gam and greenhouse air temperatures recorded every 15 minutes (contained in data GH_r1), 
  ##        predict the soil temperature every 15 minutes from the start of the first run to the end of day 55.
  iGeno_result <- GH_r1 %>% mutate(Genotype = as.character(iGeno$Genotype[1]), Soil_temp = predict(Soil.gam, newdata = GH_r1))

  ## Step 4: given estimates of Flux.gam in Step 1 and those soil temperature estimates obtained in Step 3, 
  ##        predict the CO2 flux every 15 min from the start of the first run to the end of day 55.
  iGeno_result <- iGeno_result %>% mutate(Flux_pred = exp(predict(Flux.gam, newdata = iGeno_result)))
  
  run1_result <- rbind(run1_result, iGeno_result) ## add the predicted CO2 flux for ith genotype to run1_result
}


run2_result <- NULL ## run2_results will store the predicted CO2 emission for the second run

for(i in 1:19){ ## go through genotype 1 to 19
  
  ## iGeno: carbon emission and greenhouse data for the ith genotype
  iGeno <-  r2 %>% filter(Genotype == Genotype[i]) 
  
  ## Step 1: model the relationship between CO2 emission and soil temperature using gam function (formula (1) in Supplementary Materials)
  Flux.gam <- gam(log(Exp_Flux) ~ s(Day, bs = "cr", k=18) + Soil_temp,  select=TRUE, data = iGeno)
  
  ## Step 2: model the relationship between soil temperature and air temperature using gam function (formula (2) in Supplementary Materials)
  Soil.gam = gam(Soil_temp~s+c+s(Day, bs="cr")+s(GH_temp_lag1, bs="cr"), select=TRUE, data = iGeno)
  
  ## Step 3: given estimates of Soil.gam and greenhouse air temperatures recorded every 15 minutes (contained in data GH_r2), 
  ##        predict the soil temperature every 15 minutes from the start of the second run to the end of day 55.
  iGeno_result <- GH_r2 %>% mutate(Genotype = as.character(iGeno$Genotype[1]), Soil_temp = predict(Soil.gam, newdata = GH_r2))
  
  ## Step 4: given estimates of Flux.gam in Step 1 and those soil temperature estimates obtained in Step 3, 
  ##        predict the CO2 flux every 15 min from the start of the second run to the end of day 55.
  iGeno_result <- iGeno_result %>% mutate(Flux_pred = exp(predict(Flux.gam, newdata = iGeno_result)))
  
  run2_result <- rbind(run2_result, iGeno_result) ## add the predicted CO2 flux for ith genotype to run2_result
}


## For each genotype, subtract the estimated CO2 flux values by the corresponding estimated CO2 flux of the blank control within each run.
run1_blank <- run1_result %>% filter(Genotype == "Blank") 
run1_result_adjusted <- run1_result %>% group_by(Genotype) %>% mutate(Flux = Flux_pred - run1_blank$Flux_pred, Flux = pmax(0, Flux))
run2_blank <- run2_result %>% filter(Genotype == "Blank") 
run2_result_adjusted <- run2_result %>% group_by(Genotype) %>% mutate(Flux = Flux_pred - run2_blank$Flux_pred, Flux = pmax(0, Flux))

## Calculate the area under each CO2 emission curve (AUC) for both runs using equation (3) in Supplementary Materials.
AUC_r1 <- run1_result_adjusted %>%
            group_by(Genotype) %>%
            summarise(TotalFlux = trapz(Day,Flux)) %>%
            mutate(C_released = TotalFlux*24*60*60*12*10^-6, Run = 1)
AUC_r2 <- run2_result_adjusted %>%
            group_by(Genotype) %>%
            summarise(TotalFlux = trapz(Day,Flux)) %>%
            mutate(C_released = TotalFlux*24*60*60*12*10^-6, Run = 2)
AUC <- data.frame(rbind(AUC_r1, AUC_r2)) ## combine AUCs for the two runs.
AUC <- AUC %>% filter(Genotype != "Blank") ## the AUC for Blank should be zero, remove Blank from AUC data set.



## --------------------------------------------------------------------------------------------------------------------------------------- ##
## -------------------- Section: Testing genotype difference in total carbon emission in Supplementary Materials. ------------------------ ##
## --------------------------------------------------------------------------------------------------------------------------------------- ##

## Fit a linear model with additive run and genotype effects to AUC data.
AUC.lm <- lm(C_released ~ factor(Run) + Genotype, data = AUC)
anova(AUC.lm)
## Tukey-adjusted multiple comparisons .
aov.out = aov(AUC.lm)
TukeyHSD(aov.out, which=c("Genotype"), conf.level=.95)



## --------------------------------------------------------------------------------------------------------------------------------------- ##
## ---------------------------------------- Section: Lignin analysis in Supplementary Materials. ----------------------------------------- ##
## --------------------------------------------------------------------------------------------------------------------------------------- ##

## merge AUC and LigninData by genotype and run.
dat <- merge(AUC, LigninData, by = c("Genotype", "Run")) 
Lig_r1 <- dat %>% arrange(Run, Genotype) %>% filter(Run == 1)
Lig_r2 <- dat %>% arrange(Run, Genotype) %>% filter(Run == 2)

## Lig_avg (Table 1 in Main Text): estimated total carbon emission and lignin values of 18 genotypes averaged across two runs.
Lig_avg <- (Lig_r1[, -c(1:2)] + Lig_r2[, -c(1:2)])/2
Lig_avg <- Lig_avg %>% mutate(SG_rat = S/G, Genotype = Lig_r1$Genotype)

## create scatterplots of total carbon emission vs lignin content.
col.names <- c("C_released", "G", "S", "SG_rat", "SplusG", "avg_Klason", "C_content")
pairs(Lig_r1[, col.names])
pairs(Lig_r2[, col.names])
pairs(Lig_avg[, col.names])


## the correlation test between total carbon emission and lignin measurements.
rcorr(as.matrix(Lig_r1[, col.names]))
rcorr(as.matrix(Lig_r2[, col.names]))
rcorr(as.matrix(Lig_avg[, col.names])) ## produce Table S1 in Supplementary Materials



## --------------------------------------------------------------------------------------------------------------------------------------- ##
## ----------------------------------------------------------------- Figures ------------------------------------------------------------- ##
## --------------------------------------------------------------------------------------------------------------------------------------- ##

## Fig.1 in Main Text: plot the estimated CO2 flux curves for three genotypes (CML247, Oh7B and blank) in the first run.
estimated.data <- run1_result %>% filter(Genotype == "CML247" | Genotype == "Blank" | Genotype == "Oh7B")
raw.data <- r1 %>% filter(Genotype == "CML247" | Genotype == "Blank" | Genotype == "Oh7B")
custompal <- c("grey", "red", "black")

## The main figure
Fig.1 <- ggplot(data = estimated.data, aes(Day, Flux_pred, color = Genotype), xlab = "Day", ylab = "CO2 emission") + 
         geom_line() +
         geom_point(data = raw.data, aes(Day, Exp_Flux), size = 3) +
         scale_y_continuous("CO2 Emission (umol/m2/s)") +
         scale_x_continuous(breaks = seq(0, 55, by = 10)) +
         theme_bw() +
         theme(panel.border = element_blank(), 
                  axis.line = element_line(),
           panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank(), 
                 legend.key = element_blank()) + scale_colour_manual(values=custompal)
Fig.1 

## The first inset
Fig.1.Inset.1 <- ggplot(data = estimated.data, aes(Day, Flux_pred, color = Genotype), xlab = "Day", ylab = "CO2 emission") + 
         geom_line() +
         geom_point(data = raw.data, aes(Day, Exp_Flux), size = 3) +
         scale_y_continuous("CO2 Emission (umol/m2/s)", limit=c(0,15)) +
         scale_x_continuous(limit=c(10,20), breaks=seq(10,20, by=2)) +
         theme_bw() +
         theme(panel.border = element_blank(), 
                  axis.line = element_line(),
           panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank(), 
                 legend.key = element_blank()) + scale_colour_manual(values=custompal)
Fig.1.Inset.1

## The second inset
Fig.1.Inset.2 <- ggplot(data = estimated.data, aes(Day, Flux_pred, color = Genotype), xlab = "Day", ylab = "CO2 emission") + 
         geom_line() +
         geom_point(data = raw.data, aes(Day, Exp_Flux), size = 3) +
         scale_y_continuous("CO2 Emission (umol/m2/s)", limit=c(0,15)) +
         scale_x_continuous(limit=c(40,50), breaks=seq(40,50, by=2)) +
         theme_bw() +
         theme(panel.border = element_blank(), 
                  axis.line = element_line(),
           panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank(), 
                 legend.key = element_blank()) + scale_colour_manual(values=custompal)
Fig.1.Inset.2


## Fig.S1 in Supplementary Materials: plot the estimated CO2 flux curves for 19 genotypes in the first run. 
Fig.S1 <- ggplot(data = run1_result, aes(Day, Flux_pred, color = Genotype)) + 
          geom_line() +
          geom_point(data = r1, aes(Day, Exp_Flux)) +
          scale_y_continuous("Run 1 CO2 Emission (umol/m2/s)", limit = c(0, 220)) +
          scale_x_continuous(breaks = seq(0, 55, by = 10)) +
          theme_bw() +
          theme(panel.border = element_blank(), 
                   axis.line = element_line(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
                  legend.key = element_blank()) 
Fig.S1 

## Fig.S2 in Supplementary Materials: plot the estimated CO2 flux curves for 19 genotypes in the second run. 
Fig.S2 <- ggplot(data = run2_result, aes(Day, Flux_pred, color = Genotype)) + 
          geom_line() +
          geom_point(data = r2, aes(Day, Exp_Flux)) +
          scale_y_continuous("Run 2 CO2 Emission (umol/m2/s)", limit = c(0, 220)) +
          scale_x_continuous(breaks=seq(0, 55, by = 10)) +
          theme_bw() +
          theme(panel.border = element_blank(), 
                   axis.line = element_line(),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
                  legend.key = element_blank()) 
Fig.S2 
