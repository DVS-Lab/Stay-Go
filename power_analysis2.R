library(plyr)                                   
library(dplyr)
library(tidyverse)
library(readr)      
library(haven)

library(lme4)
library(lmerTest)
library(simr)

# Daniel Sazhin
# 01/26/2023
# Neuroeconomics Lab
# Temple University

# This code does the power analysis for the StayGo2 Project. We do this by 
# importing data from StayGo1 and complete a power analysis for the variables 
# we already know. 

setwd('A:/Data/Stay-Go/')

# Import Mike McCloskeys AUDIT data.

AUDIT_all <- read_sav('A:/Data/Stay-Go/ISTART_Master_Merged.sav') # Use AUDIT_tot
AUDIT <- data.frame(AUDIT_all$AUDIT_tot)

# Import the data.

data_all <- list.files(path = 'A:/Data/Stay-Go/Participants',  # Identify all CSV files
                       pattern = "*.csv", full.names = TRUE) %>% 
  lapply(read_csv) %>%                              # Store all files in list
  bind_rows                                         # Combine data sets into one data set 
data_all                                            # Print data to RStudio console

# Make an ID for each participants.

# Generate an AUDIT dataset


subset <- AUDIT[complete.cases(nrow(AUDIT), 1000), ]
subset <- AUDIT[sample(nrow(AUDIT), 123), ]

AUDIT_test <- data.frame()
ID <- data.frame()

N = 123

for (ii in 1:N) {

AUDIT_process = rep(subset[ii],55)
AUDIT_process = matrix(AUDIT_process)

subject = rep(ii,55)
subject= matrix(subject)

ID = rbind(ID,subject) 
AUDIT_test = rbind(AUDIT_test,AUDIT_process)
}

colnames(ID) <- c("ID")
colnames(AUDIT_test) <- c("AUDIT_test")


# Concatenate the IDs to be tied to each participant

Full_Set = cbind(data_all, ID, AUDIT_test)

Full_Set_Double = rbind(Full_Set, Full_Set) # N = 246

Full_Set_Triple = rbind(Full_Set, Full_Set_Double) # N = 369

#random_slope_model = fitlme(mixed_effects,'turns ~ value * growth + (value | ID) + (growth | ID)')

# Generating a simple HLM model. This one tests using AUDIT and Final Value as the main regressors, with ID as the random variable.

# GLMER (Poission distribution for dependent variable. 80 percent power.)
# Pilot phase
# Kreft, I. G. G. 1996. Are multilevel techniques necessary? An overview, including simulation studies. Los Angeles, CA: California State University at Los Angeles.
# Kreft (1996) suggested 30 Level 2 groups and 30 Level 1 observations per group. Simulation studies indicate that the balance swings toward more than 30 groups and fewer than 30 observations per group for hypotheses about the effects of Level 2 variables (e.g., Maas & Hox, 2005; Snijders & Bosker, 1993). For a two-level model, Hox (2010) recommends at least 20 observations for 50 groups to test cross-level interactions, and at least 10 observations for 100 groups to test random effects. He also presents detailed power analyses based on both sample sizes and anticipated effect sizes. Spybrook, Raudenbush, Congdon, and Martínez (2011) provide software for power analysis for specific MLM situations.




m1 <- lmer(Full_Set$Turn_Left ~ Full_Set$AUDIT_test * Full_Set$Final_Value + (1 | Full_Set$ID), 
           data = Full_Set,
           REML = F)


# Testing Power for VAR2
powerSim(m1, 
         nsim=1000, 
         test = fcompare(Full_Set$Turn_Left ~ Full_Set$AUDIT_test + Full_Set$Final_Value:Full_Set$AUDIT_test))

# Power = 59.40%

# Testing Power for VAR1-VAR2 Interaction 
powerSim(m1, 
         nsim=1000, 
         test = fcompare(Full_Set$Turn_Left ~ Full_Set$AUDIT_test + Full_Set$Final_Value))

# Power = 49.80%

# Double the N.

m1 <- lmer(Full_Set_Double$Turn_Left ~ Full_Set_Double$AUDIT_test * Full_Set_Double$Final_Value + (1 | Full_Set_Double$ID), 
           data = Full_Set_Double,
           REML = F)

powerSim(m1, 
         nsim=1000, 
         test = fcompare(Full_Set_Double$Turn_Left ~ Full_Set_Double$AUDIT_test + Full_Set_Double$Final_Value:Full_Set_Double$AUDIT_test))

# Power = 66.50%

powerSim(m1, 
         nsim=1000, 
         test = fcompare(Full_Set_Double$Turn_Left ~ Full_Set_Double$AUDIT_test + Full_Set_Double$Final_Value))

# Power = 48.60%

m1a <- glmer(Full_Set$Turn_Left ~ Full_Set$AUDIT_test + (1 | Full_Set$ID), data = Full_Set, family = poisson(link = "log"))

# Testing Power for VAR2
powerSim(m1a, 
         nsim=1000, 
         test = fcompare(Full_Set$Turn_Left ~ Full_Set$AUDIT_test + Full_Set$Final_Value:Full_Set$AUDIT_test))

# Power = 59.40%

# Testing Power for VAR1-VAR2 Interaction 
powerSim(m1a, 
         nsim=1000, 
         test = fcompare(Full_Set$Turn_Left ~ Full_Set$AUDIT_test + Full_Set$Final_Value))

# Power = 49.80%


# Model 2: AUDIT and Alpha only 

m2 <- lmer(Full_Set$Turn_Left ~ Full_Set$AUDIT_test * Full_Set$Alpha + (1 | Full_Set$ID), 
           data = Full_Set,
           REML = F)

powerSim(m2, 
         nsim=1000, 
         test = fcompare(Full_Set$Turn_Left ~ Full_Set$AUDIT_test + Full_Set$Alpha:Full_Set$AUDIT_test))

# Power = 100.0% 

# Testing Power for VAR1-VAR2 Interaction
powerSim(m2, 
         nsim=1000, 
         test = fcompare(Full_Set$Turn_Left ~ Full_Set$AUDIT_test + Full_Set$Alpha))

# Power = 53.20% 

# Model 3: Final Values and Alpha only 

m3 <- lmer(Full_Set$Turn_Left ~ Full_Set$Final_Value * Full_Set$Alpha + (1 | Full_Set$ID), 
           data = Full_Set,
           REML = F)

powerSim(m3, 
         nsim=1000, 
         test = fcompare(Full_Set$Turn_Left ~ Full_Set$Final_Value + Full_Set$Alpha:Full_Set$Final_Value))

# Power = 74.10%

# Testing Power for VAR1-VAR2 Interaction
powerSim(m3, 
         nsim=1000, 
         test = fcompare(Full_Set$Turn_Left ~ Full_Set$Final_Value + Full_Set$Alpha))

# Power = 49.50%
