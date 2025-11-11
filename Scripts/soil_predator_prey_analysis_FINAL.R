###########################################
# Predator-Prey Analysis and Visualization - Ambient & Heatwave Temperature Treatments
# The analysis entails running linear models for treatment over time, post-hoc analysis
# on significant results, and data visualization for figures.
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
  drop_na()

# Heatwave
heatwave <- read_csv(heatwave_temp)%>%
  mutate(predator = staphylinidae+staphylinidae_larval, Temperature = "Heatwave")%>%
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

## PREP DATA TO CREATE FIGURES
# Filter data to only include microcosms where predators present + count only predator abundance
rove_relabel<-
  all_data%>%
  filter(predation_trt == "CB")%>%
  select(-collembola)%>%
  mutate(predation_trt = dplyr::recode(predation_trt, "CB" = "predators"))

# Filter data to only count prey abundance in microcosms and relabel predation treatment
# to indicate presence or absence of predation pressure
collem_relabel<-
  all_data%>%
  select(-predator)%>%
  mutate(predation_trt = dplyr::recode(predation_trt, "C" = "collembola_only", "CB" = "collem_predator"))

# Merge datasets
clean_data<-rbind.fill(collem_relabel, rove_relabel)

# Fill in NAs with zeroes
clean_data[is.na(clean_data)] <- 0

final_figure_data<-
  clean_data%>%
  mutate(Total = collembola+predator)

# Find mean collembola/staphylinidae count for treatments as a function of time
all_figure_summary <- summarySE(final_figure_data, measurevar = "Total", groupvars = c("site", "harvest_week", "predation_trt", "Temperature"))

## CREATE LINE CHART SHOWING AVG # OF SOIL ARTHROPODS PER WEEK OF EXPERIMENT - AMBIENT TEMP TRT
# Find mean collembola/staphylinidae count for ambient temp trt only
ambient_figure<-
  all_figure_summary%>%
  filter(Temperature == "Ambient")

# Add "Week 0" to ambient data to show initial microcosm prey + predator counts in figure to show starting conditions
## Start with agricultural site:
week0.ambient.ag <- data.frame(c("Intensive Agriculture","Intensive Agriculture","Intensive Agriculture"), 
                               c(0,0,0),
                               c("collem_predator","collembola_only","predators"),
                               c("Ambient","Ambient","Ambient"), 
                               c(3,3,3),
                               c(20,20,5),
                               c(0,0,0),
                               c(0,0,0),
                               c(0,0,0))
names(week0.ambient.ag) <- c("site","harvest_week","predation_trt","Temperature","N","Total","sd","se","ci")
# Add other sites:
week0.ambient.urb <- data.frame(c("Urbanized","Urbanized","Urbanized"), 
                                c(0,0,0),
                                c("collem_predator","collembola_only","predators"),
                                c("Ambient","Ambient","Ambient"), 
                                c(3,3,3),
                                c(20,20,5),
                                c(0,0,0),
                                c(0,0,0),
                                c(0,0,0))
names(week0.ambient.urb) <- c("site","harvest_week","predation_trt","Temperature","N","Total","sd","se","ci")
week0.ambient.grass <- data.frame(c("Grassland","Grassland","Grassland"), 
                                  c(0,0,0),
                                  c("collem_predator","collembola_only","predators"),
                                  c("Ambient","Ambient","Ambient"), 
                                  c(3,3,3),
                                  c(20,20,5),
                                  c(0,0,0),
                                  c(0,0,0),
                                  c(0,0,0))
names(week0.ambient.grass) <- c("site","harvest_week","predation_trt","Temperature","N","Total","sd","se","ci")
week0.ambient <- bind_rows(week0.ambient.ag,week0.ambient.urb,week0.ambient.grass)
ambient_figure.allweeks <- rbind(ambient_figure, week0.ambient)

# Create figure for ambient temperature treatment
ambient_figures <- ambient_figure.allweeks %>% 
  mutate(Treatment = fct_relevel(predation_trt, 
                                 "collembola_only", "collem_predator", "predators"))%>%
  mutate(Location = fct_relevel(site, "Grassland", "Intensive Agriculture", "Urbanized"))%>%
  ggplot(aes(x=harvest_week, y=Total, group = Treatment)) +
  geom_line(aes(color = Treatment, group = Treatment)) +
  geom_point(aes(color = factor(Treatment)), size = 2.5) +
  geom_errorbar(aes(ymin = Total - se, ymax = Total + se, color = Treatment), size = 1, width = 0) + # adds error bars
  xlab("Microcosm harvest week") +
  ylab("Average abundance \n (arthropods/microcosm)") +
  scale_y_continuous(limits = c(0, 90))+
  # wasn't showing one of the error bars because the y-axis limit was too small
  scale_color_manual(values=c("#7BAEA0", "#386276", "#D9B96E"))+
  theme(axis.text.x = element_text(colour= "black", face = "bold", size = 10),
        axis.text.y = element_text(colour= "black", face = "bold", size = 10),
        axis.line = element_line(colour = "black", size = .3),
        axis.line.x = element_line(colour = "black", size = .3),
        axis.ticks.x = element_line(colour = "black", size = 1),
        axis.ticks.y = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size = .3),
        axis.title=element_text(size=12),
        strip.text = element_text(size = 12, color = 'black', face = "bold"),
        strip.background = element_rect(fill = "grey"),
        legend.position = "none",
        panel.background = element_rect(fill= "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = 0.5))+ 
  scale_x_continuous(breaks = seq(1, 3, 1)) + 
  facet_grid(~Location)

ambient_figures

## CREATE LINE CHART SHOWING AVG # OF SOIL ARTHROPODS PER WEEK OF EXPERIMENT - AMBIENT TEMP TRT
heatwave_figure<-
  all_figure_summary%>%
  filter(Temperature == "Heatwave")

# Add "Week 0" to heatwave data
## Start with agricultural data:
week0.heatwave.ag <- data.frame(c("Intensive Agriculture","Intensive Agriculture","Intensive Agriculture"), 
                                c(0,0,0),
                                c("collem_predator","collembola_only","predators"),
                                c("Heatwave","Heatwave","Heatwave"), 
                                c(3,3,3),
                                c(20,20,5),
                                c(0,0,0),
                                c(0,0,0),
                                c(0,0,0))
names(week0.heatwave.ag) <- c("site","harvest_week","predation_trt","Temperature","N","Total","sd","se","ci")
week0.heatwave.urb <- data.frame(c("Urbanized","Urbanized","Urbanized"), 
                                 c(0,0,0),
                                 c("collem_predator","collembola_only","predators"),
                                 c("Heatwave","Heatwave","Heatwave"), 
                                 c(3,3,3),
                                 c(20,20,5),
                                 c(0,0,0),
                                 c(0,0,0),
                                 c(0,0,0))
names(week0.heatwave.urb) <- c("site","harvest_week","predation_trt","Temperature","N","Total","sd","se","ci")
week0.heatwave.grass <- data.frame(c("Grassland","Grassland","Grassland"), 
                                   c(0,0,0),
                                   c("collem_predator","collembola_only","predators"),
                                   c("Heatwave","Heatwave","Heatwave"), 
                                   c(3,3,3),
                                   c(20,20,5),
                                   c(0,0,0),
                                   c(0,0,0),
                                   c(0,0,0))
names(week0.heatwave.grass) <- c("site","harvest_week","predation_trt","Temperature","N","Total","sd","se","ci")
week0.heatwave <- bind_rows(week0.heatwave.ag,week0.heatwave.urb,week0.heatwave.grass)
heatwave_figure.allweeks <- rbind(heatwave_figure, week0.heatwave)

# Create figure for heatwave temperature treatment
heatwave_figures <- heatwave_figure.allweeks %>% 
  mutate(Treatment = fct_relevel(predation_trt, 
                                 "collembola_only", "collem_predator", "predators"))%>%
  mutate(Location = fct_relevel(site, "Grassland", "Intensive Agriculture", "Urbanized"))%>%
  ggplot(aes(x=harvest_week, y=Total, group = Treatment)) +
  geom_line(aes(color = Treatment, group = Treatment)) +
  geom_point(aes(color = factor(Treatment)), size = 2.5) +
  geom_errorbar(aes(ymin = Total - se, ymax = Total + se, color = Treatment), size = 1, width = 0) + # adds error bars
  xlab("Microcosm harvest week") +
  ylab("Average abundance \n (arthropods/microcosm)") +
  scale_y_continuous(limits = c(0, 90))+
  scale_color_manual(values=c("#7BAEA0", "#386276", "#D9B96E"))+
  theme(axis.text.x = element_text(colour= "black", face = "bold", size = 10),
        axis.text.y = element_text(colour= "black", face = "bold", size = 10),
        axis.line = element_line(colour = "black", size = .3),
        axis.line.x = element_line(colour = "black", size =.3),
        axis.ticks.x = element_line(colour = "black", size = 1),
        axis.ticks.y = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size =.3),
        axis.title=element_text(size=12),
        strip.text = element_text(size = 12, color = 'black', face = "bold"),
        strip.background = element_rect(fill = "grey"),
        legend.position = "none",
        panel.background = element_rect(fill= "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = 0.5))+ 
  scale_x_continuous(breaks = seq(1, 3, 1)) + 
  facet_grid(~Location)

heatwave_figures

## Save figures
# Ambient
ggsave(file = "ambient_figures.jpg", plot = ambient_figures, width=7.5, height=3, units="in", dpi=1000)
# Heatwave
ggsave(file = "heatwave_figures.jpg", plot = heatwave_figures, width=7.5, height=3, units="in", dpi=1000)

#-------------------------------------------------------------------------------------------

## ANALYSIS OF FINAL ARTHROPOD ABUNDANCES @ END OF EXPERIMENT (WEEK 3)
# Organize data
week3<-all_data %>% 
  filter(harvest_week == 3)

###Linear models of the last week of the experiment
mod5<-aov(log(collembola+1)~predation_trt*site*Temperature, data = week3)
anova(mod5) # Results printed in Supplemental Index: Table S3
# And again for prey abundance in heatwave trt

###Linear models of the last week of the experiment for predatprs
mod6<-aov(log(predator+1)~site*Temperature, data = week3)
anova(mod6) # Results printed in Supplemental Index: Table S3
# And again for prey abundance in heatwave trt

#Posthoc test
tukey5 <- TukeyHSD(mod5)
tukey6 <- TukeyHSD(mod6)

# Summary statistics for prey
sum_week3<-summarySE(data = week3, measurevar = "collembola", groupvars = c("predation_trt"))# Summary statistics

# Summary statistics for predator
sum_week3_predator<-summarySE(data = week3, measurevar = "predator", groupvars = c("Temperature", "site"))


## Visualize data for final arthropod abundance in bar graph
# Colorblind-friendly color palette
coast_palette <- pal_vector("coast",6)

week3_figure2 <- ggplot(sum_week3, aes(x=predation_trt, y=collembola, fill=Temperature)) +
  geom_bar(stat="identity", position = "dodge") +
  geom_errorbar(aes(ymin=collembola-se, ymax=collembola+se), width=0.2, position = position_dodge(width=0.9)) +
  ggtitle("Avg overall Collembola abundance at week 3") + 
    # aka effects on predation x temp x predation_trt interaction on collembola abundance at conclusion of experiment
  xlab("Predation pressure") +
  ylab("Average abundance \n(Collembola/microcosm)") +
  scale_x_discrete(labels = c("Absent", "Present")) +
  scale_y_continuous(limits = c(0, 90)) +
  scale_fill_manual(values = c("#CC7810","#EFC519"))+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        axis.text.x = element_text(colour= "black", face = "bold", size = 8),
        axis.text.y = element_text(colour= "black", face = "bold", size = 8),
        axis.line = element_line(colour = "black", size = .3),
        axis.line.x = element_line(colour = "black", size = .3),
        axis.ticks.x = element_line(colour = "black", size = 1),
        axis.ticks.y = element_line(colour = "black", size = 1),
        axis.line.y = element_line(colour = "black", size = .3),
        axis.title=element_text(size=12),
        strip.text = element_text(size = 10, color = 'black', face = "bold"),
        strip.background = element_rect(fill = "grey"),
        legend.position = "none",
        #legend.title = element_text(color = "black", size = 10),
        #legend.text = element_text(color = "black", size = 8),
        panel.background = element_rect(fill= "white"),
        panel.border = element_rect(fill = NA, colour = "black", size = 0.5)) +
  facet_wrap(~ site)

week3_figure2

# Save bar graph
ggsave(file = "week3_figure2.jpg", plot = week3_figure2, width=5, height=3, units="in", dpi=1000)



