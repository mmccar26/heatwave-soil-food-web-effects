###########################################
# Predator-Prey Analysis for the code no Sminthuridae contimination
###########################################

# Load packages:
library(dplyr)
library(tidyverse)
library(lme4)
library(ggplot2)
library(Rmisc)
library(forcats)
library(jtools)
library(broom.mixed)
library(wacolors)

## IMPORT DATASETS FOR AMBIENT AND HEATWAVE TEMP TREATMENTS
# Relative pathname
ambient_temp <- file.path(".", "Data", "ambient_temp_final.csv")
heatwave_temp<- file.path(".", "Data", "heatwave_temp_final.csv")
print(ambient_temp)

## LOAD DATA WITH RIGHT VARIABLES
# Ambient
ambient <- read_csv(ambient_temp)%>%
  mutate(predator = staphylinidae+staphylinidae_larval, Temperature = "Ambient")%>%
  select(-collembola_overhead, -larval_overhead, -oribatida, -staphylinidae, -staphylinidae_larval, -staphylinidae_overhead)%>%
  drop_na()

# Heatwave
heatwave <- read_csv(heatwave_temp)%>%
  filter(sminthuridae <= 20)%>%
  mutate(predator = staphylinidae+staphylinidae_larval, Temperature = "Heatwave")%>%
  select(-date_processed, -sminthuridae, -oribatida, -collembola_overhead, -staphylinidae, -staphylinidae_larval, -staphylinidae_overhead, -larval_overhead,
         -mesostig, -prostig)%>%
  drop_na()

# Combine heatwave and ambient datasets
all_data <- rbind(ambient, heatwave)

# Subset data for Staphylinid predators only
rove<-
  all_data%>%
  filter(predation_trt == "CB")

## LINEAR MODELS
#  Test of collembola: prey abundance across all treatments
mod1<-aov(log(collembola+1)~Temperature*predation_trt*harvest_week*site, data = all_data)
anova(mod1) # Results printed in Table 2 of main text
# Temperature and site affect how predators influence collembola abundance (interactive effects 0.002)
plot(mod1)
# Rerun analysis to just look at prey abundance in ambient temp trt only
mod3<-aov(log(collembola+1)~predation_trt*harvest_week*site, data = ambient)
anova(mod3) # Results printed in Supplemental Index: Table S3
# And again for prey abundance in heatwave trt
mod4<-aov(log(collembola+1)~predation_trt*harvest_week*site, data = heatwave)
anova(mod4) # Results printed in Supplemental Index: Table S4

# Test of rove beetles: rove beetle abundance across all treatments when present
mod2<-aov(log(predator+1)~Temperature*harvest_week*site, data = rove)
anova(mod2) # Results printed in Table 2 of main text

## RUN POST-HOC ANALYSES TO TEASE OUT INTERACTIONS FROM MODELS
tukey1 <- TukeyHSD(mod1)
tukey2 <- TukeyHSD(mod2)

# Save Tukey results to create tables in Supplemental Appendix
# Table S1
tukey1_data <- lapply(names(tukey1), function(name) {
  df <- as.data.frame(tukey1[[name]])
  df$Comparison <- rownames(df)
  df$Factor <- name
  return(df)
})
tukey1_all <- do.call(rbind, tukey1_data)
tukey1_all <- tukey1_all[, c("Factor", "Comparison", "diff", "lwr", "upr", "p adj")]
write.csv(tukey1_all, file = "tukey1_results.csv", row.names = FALSE)
# Table S2
tukey2_data <- lapply(names(tukey2), function(name) {
  df <- as.data.frame(tukey2[[name]])
  df$Comparison <- rownames(df)
  df$Factor <- name
  return(df)
})
tukey2_all <- do.call(rbind, tukey2_data)
tukey2_all <- tukey2_all[, c("Factor", "Comparison", "diff", "lwr", "upr", "p adj")]
write.csv(tukey2_all, file = "tukey2_results.csv", row.names = FALSE)

