library(ggplot2)
library(readxl)
library(reshape2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(ggalluvial)

setwd("/Users/roucogar/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/1_FACS-sorter/Flowjo_analysis_nov2023/")

#### 58 at different stages datasets ####
#E105
E105_normalised <- read.csv("E105/58vsmutantclones_E105_normalised.csv")
E105_normalised <- E105_normalised[E105_normalised$clone == "58",]
E105_normalised$stage <- "E105"
all_stages58_normalised <- E105_normalised

E105_FL_averages_data <- read.csv("E105/58vsmutantclones_E105_average_normalised_FL.csv")
E105_FL_averages_data <- E105_FL_averages_data[E105_FL_averages_data$clone == "58",]
E105_FL_averages_data$stage <- "E105"
all_stages58_normalised_avg_FL <- E105_FL_averages_data

E105_HL_averages_data <- read.csv("E105/58vsmutantclones_E105_average_normalised_HL.csv")
E105_HL_averages_data <- E105_HL_averages_data[E105_HL_averages_data$clone == "58",]
E105_HL_averages_data$stage <- "E105"
all_stages58_normalised_avg_HL <- E105_HL_averages_data

rm(E105_normalised, E105_FL_averages_data,E105_HL_averages_data)

#E110
E110_normalised <- read.csv("E110/58vsmutantclones_E110_normalised.csv")
E110_normalised <- E110_normalised[E110_normalised$clone == "58",]
E110_normalised$stage <- "E110"
all_stages58_normalised <- rbind(all_stages58_normalised,E110_normalised)

E110_FL_averages_data <- read.csv("E110/58vsmutantclones_E110_average_normalised_FL.csv")
E110_FL_averages_data <- E110_FL_averages_data[E110_FL_averages_data$clone == "58",]
E110_FL_averages_data$stage <- "E110"
all_stages58_normalised_avg_FL <- rbind(all_stages58_normalised_avg_FL,E110_FL_averages_data)

E110_HL_averages_data <- read.csv("E110/58vsmutantclones_E110_average_normalised_HL.csv")
E110_HL_averages_data <- E110_HL_averages_data[E110_HL_averages_data$clone == "58",]
E110_HL_averages_data$stage <- "E110"
all_stages58_normalised_avg_HL <- rbind(all_stages58_normalised_avg_HL,E110_HL_averages_data)

rm(E110_normalised, E110_FL_averages_data,E110_HL_averages_data)

#E115
E115_normalised <- read.csv("E115/58vsmutantclones_E115_normalised.csv")
E115_normalised <- E115_normalised[E115_normalised$clone == "58",]
E115_normalised$stage <- "E115"
all_stages58_normalised <- rbind(all_stages58_normalised,E115_normalised)

E115_FL_averages_data <- read.csv("E115/58vsmutantclones_E115_average_normalised_FL.csv")
E115_FL_averages_data <- E115_FL_averages_data[E115_FL_averages_data$clone == "58",]
E115_FL_averages_data$stage <- "E115"
all_stages58_normalised_avg_FL <- rbind(all_stages58_normalised_avg_FL,E115_FL_averages_data)

E115_HL_averages_data <- read.csv("E115/58vsmutantclones_E115_average_normalised_HL.csv")
E115_HL_averages_data <- E115_HL_averages_data[E115_HL_averages_data$clone == "58",]
E115_HL_averages_data$stage <- "E115"
all_stages58_normalised_avg_HL <- rbind(all_stages58_normalised_avg_HL,E115_HL_averages_data)

rm(E115_normalised, E115_FL_averages_data,E115_HL_averages_data)

#E125 + E130
E125_E130_normalised <- read.csv("E125_E130/58vsmutantclones_E125_E130_normalised.csv")
E125_E130_normalised <- E125_E130_normalised[17:25,c(1:14,16,17,15)]
head(E125_E130_normalised)
all_stages58_normalised <- rbind(all_stages58_normalised,E125_E130_normalised)

E125_FL_averages_data <- read.csv("E125_E130/58vsmutantclones_E125_average_normalised_FL.csv")
E125_FL_averages_data$stage <- "E125"
all_stages58_normalised_avg_FL <- rbind(all_stages58_normalised_avg_FL,E125_FL_averages_data)

E125_HL_averages_data <- read.csv("E125_E130/58vsmutantclones_E125_average_normalised_HL.csv")
E125_HL_averages_data$stage <- "E125"
all_stages58_normalised_avg_HL <- rbind(all_stages58_normalised_avg_HL,E125_HL_averages_data)

rm (E125_E130_normalised, E125_FL_averages_data,E125_HL_averages_data)

E130_FL_averages_data <- read.csv("E125_E130/58vsmutantclones_E130_average_normalised_FL.csv")
E130_FL_averages_data$stage <- "E130"
all_stages58_normalised_avg_FL <- rbind(all_stages58_normalised_avg_FL,E130_FL_averages_data)

E130_HL_averages_data <- read.csv("E125_E130/58vsmutantclones_E130_average_normalised_HL.csv")
E130_HL_averages_data$stage <- "E130"
all_stages58_normalised_avg_HL <- rbind(all_stages58_normalised_avg_HL,E130_HL_averages_data)

rm(E130_FL_averages_data,E130_HL_averages_data)

#E135
E135_normalised <- read.csv("E135/58vsmutantclones_E135_normalised.csv")
E135_normalised <- E135_normalised[E135_normalised$clone == "58",]
E135_normalised$stage <- "E135"
all_stages58_normalised <- rbind(all_stages58_normalised,E135_normalised)

E135_FL_averages_data <- read.csv("E135/58vsmutantclones_E135_average_normalised_FL.csv")
E135_FL_averages_data <- E135_FL_averages_data[E135_FL_averages_data$clone == "58",]
E135_FL_averages_data$stage <- "E135"
all_stages58_normalised_avg_FL <- rbind(all_stages58_normalised_avg_FL,E135_FL_averages_data)

E135_HL_averages_data <- read.csv("E135/58vsmutantclones_E135_average_normalised_HL.csv")
E135_HL_averages_data <- E135_HL_averages_data[E135_HL_averages_data$clone == "58",]
E135_HL_averages_data$stage <- "E135"
all_stages58_normalised_avg_HL <- rbind(all_stages58_normalised_avg_HL,E135_HL_averages_data)

rm(E135_normalised, E135_FL_averages_data,E135_HL_averages_data)

#E140
E140_normalised <- read.csv("E140/58vsmutantclones_E140_normalised.csv")
E140_normalised <- E140_normalised[E140_normalised$clone == "58",]
E140_normalised$stage <- "E140"
E140_normalised <-E140_normalised %>% rename(EYFP. = EYFP)
all_stages58_normalised <- rbind(all_stages58_normalised,E140_normalised)

E140_FL_averages_data <- read.csv("E140/58vsmutantclones_E140_average_normalised_FL.csv")
E140_FL_averages_data <- E140_FL_averages_data[E140_FL_averages_data$clone == "58",]
E140_FL_averages_data$stage <- "E140"
E140_FL_averages_data <-E140_FL_averages_data %>% rename(EYFP. = EYFP)
all_stages58_normalised_avg_FL <- rbind(all_stages58_normalised_avg_FL,E140_FL_averages_data)

E140_HL_averages_data <- read.csv("E140/58vsmutantclones_E140_average_normalised_HL.csv")
E140_HL_averages_data <- E140_HL_averages_data[E140_HL_averages_data$clone == "58",]
E140_HL_averages_data <-E140_HL_averages_data %>% rename(EYFP. = EYFP)
E140_HL_averages_data$stage <- "E140"
all_stages58_normalised_avg_HL <- rbind(all_stages58_normalised_avg_HL,E140_HL_averages_data)

rm(E140_normalised, E140_FL_averages_data,E140_HL_averages_data)

#E145
E145_normalised <- read.csv("E145_removing_dataset29032023/58vsmutantclones_E145_normalised_without_dataset29032023.csv")
E145_normalised <- E145_normalised[E145_normalised$clone == "58",]
E145_normalised$stage <- "E145"
all_stages58_normalised <- rbind(all_stages58_normalised,E145_normalised)

E145_FL_averages_data <- read.csv("E145_removing_dataset29032023/58vsmutantclones_E145_average_normalised_FL_without_dataset29032023.csv")
E145_FL_averages_data <- E145_FL_averages_data[E145_FL_averages_data$clone == "58",]
E145_FL_averages_data$stage <- "E145"
all_stages58_normalised_avg_FL <- rbind(all_stages58_normalised_avg_FL,E145_FL_averages_data)

E145_HL_averages_data <- read.csv("E145_removing_dataset29032023/58vsmutantclones_E145_average_normalised_HL_without_dataset29032023.csv")
E145_HL_averages_data <- E145_HL_averages_data[E145_HL_averages_data$clone == "58",]
E145_HL_averages_data$stage <- "E145"
all_stages58_normalised_avg_HL <- rbind(all_stages58_normalised_avg_HL,E145_HL_averages_data)

rm(E145_normalised, E145_FL_averages_data,E145_HL_averages_data)

#save summary tables
write.csv(all_stages58_normalised,"all_stages58_normalised.csv",row.names = FALSE)
write.csv(all_stages58_normalised_avg_FL,"all_stages58_normalised_avg_FL.csv",row.names = FALSE)
write.csv(all_stages58_normalised_avg_HL,"all_stages58_normalised_avg_HL.csv",row.names = FALSE)

#### 58 at different stages plots + selected stages ####
## plots ##
#barplot all populations
all_stages58_normalised_avg_FL$limbtype <- "FL"
all_stages58_normalised_avg_HL$limbtype <- "HL"

all_stages58_normalised_avg <- rbind(all_stages58_normalised_avg_FL,all_stages58_normalised_avg_HL)

all_stages58_normalised_avg_long <- all_stages58_normalised_avg %>% 
  gather(FACs_class, average_prop, -clone, -stage,-limbtype)

ggplot(all_stages58_normalised_avg_long, aes(fill= FACs_class, x=stage,y=average_prop, label = round(average_prop, 1))) + 
  geom_col() + facet_grid(.~limbtype)+
  geom_text(stat = "stratum", aes(stratum = FACs_class))
  
ggsave("Average_Proportion_clone58_FACsclass_from_citoflex_by_stage_and_limbtype.pdf", width = 12, height = 8)

#Figure3C
#barplot selected stages
all_stages58_normalised_avg_FL <- read.csv("/Users/roucogar/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/1_FACS-sorter/Flowjo_analysis_nov2023/58_all_stages/all_stages58_normalised_avg_FL.csv")
all_stages58_normalised_avg_HL <- read.csv("/Users/roucogar/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/1_FACS-sorter/Flowjo_analysis_nov2023/58_all_stages/all_stages58_normalised_avg_HL.csv")

all_stages58_normalised_avg_FL$limbtype <- "FL"
all_stages58_normalised_avg_HL$limbtype <- "HL"

all_stages58_normalised_avg <- rbind(all_stages58_normalised_avg_FL,all_stages58_normalised_avg_HL)

exclude_stages <- c('E110','E130','E140')
all_stages58_normalised_avg <- subset(all_stages58_normalised_avg,!(stage %in% exclude_stages))

all_stages58_normalised_avg_long <- all_stages58_normalised_avg %>% 
  gather(FACs_class, average_prop, -clone, -stage,-limbtype)

ggplot(all_stages58_normalised_avg_long, aes(fill= FACs_class, x=stage,y=average_prop, label = round(average_prop, 1))) + 
  geom_col() + facet_grid(.~limbtype)+
  geom_text(stat = "stratum", aes(stratum = FACs_class))

ggsave("Average_Proportion_clone58_selected_stages_FACsclass_from_citoflex_by_stage_and_limbtype.pdf", width = 12, height = 8)

#save summary tables
write.csv(all_stages58_normalised_avg,"all_stages58_normalised_selected_stages.csv",row.names = FALSE)

### 251 all stages datasets ####
#E105
E105_normalised <- read.csv("E105/58vsmutantclones_E105_normalised.csv")
E105_normalised <- E105_normalised[E105_normalised$clone == "251",]
E105_normalised$stage <- "E105"
all_stages251_normalised <- E105_normalised

E105_FL_averages_data <- read.csv("E105/58vsmutantclones_E105_average_normalised_FL.csv")
E105_FL_averages_data <- E105_FL_averages_data[E105_FL_averages_data$clone == "251",]
E105_FL_averages_data$stage <- "E105"
all_stages251_normalised_avg_FL <- E105_FL_averages_data

E105_HL_averages_data <- read.csv("E105/58vsmutantclones_E105_average_normalised_HL.csv")
E105_HL_averages_data <- E105_HL_averages_data[E105_HL_averages_data$clone == "251",]
E105_HL_averages_data$stage <- "E105"
all_stages251_normalised_avg_HL <- E105_HL_averages_data

rm(E105_normalised, E105_FL_averages_data,E105_HL_averages_data)

#E110
E110_normalised <- read.csv("E110/58vsmutantclones_E110_normalised.csv")
E110_normalised <- E110_normalised[E110_normalised$clone == "251",]
E110_normalised$stage <- "E110"
all_stages251_normalised <- rbind(all_stages251_normalised,E110_normalised)

E110_FL_averages_data <- read.csv("E110/58vsmutantclones_E110_average_normalised_FL.csv")
E110_FL_averages_data <- E110_FL_averages_data[E110_FL_averages_data$clone == "251",]
E110_FL_averages_data$stage <- "E110"
all_stages251_normalised_avg_FL <- rbind(all_stages251_normalised_avg_FL,E110_FL_averages_data)

E110_HL_averages_data <- read.csv("E110/58vsmutantclones_E110_average_normalised_HL.csv")
E110_HL_averages_data <- E110_HL_averages_data[E110_HL_averages_data$clone == "251",]
E110_HL_averages_data$stage <- "E110"
all_stages251_normalised_avg_HL <- rbind(all_stages251_normalised_avg_HL,E110_HL_averages_data)

rm(E110_normalised, E110_FL_averages_data,E110_HL_averages_data)

#E115
E115_normalised <- read.csv("E115/58vsmutantclones_E115_normalised.csv")
E115_normalised <- E115_normalised[E115_normalised$clone == "251",]
E115_normalised$stage <- "E115"
all_stages251_normalised <- rbind(all_stages251_normalised,E115_normalised)

E115_FL_averages_data <- read.csv("E115/58vsmutantclones_E115_average_normalised_FL.csv")
E115_FL_averages_data <- E115_FL_averages_data[E115_FL_averages_data$clone == "251",]
E115_FL_averages_data$stage <- "E115"
all_stages251_normalised_avg_FL <- rbind(all_stages251_normalised_avg_FL,E115_FL_averages_data)

E115_HL_averages_data <- read.csv("E115/58vsmutantclones_E115_average_normalised_HL.csv")
E115_HL_averages_data <- E115_HL_averages_data[E115_HL_averages_data$clone == "251",]
E115_HL_averages_data$stage <- "E115"
all_stages251_normalised_avg_HL <- rbind(all_stages251_normalised_avg_HL,E115_HL_averages_data)

rm(E115_normalised, E115_FL_averages_data,E115_HL_averages_data)

#E135
E135_normalised <- read.csv("E135/58vsmutantclones_E135_normalised.csv")
E135_normalised <- E135_normalised[E135_normalised$clone == "251",]
E135_normalised$stage <- "E135"
all_stages251_normalised <- rbind(all_stages251_normalised,E135_normalised)

E135_FL_averages_data <- read.csv("E135/58vsmutantclones_E135_average_normalised_FL.csv")
E135_FL_averages_data <- E135_FL_averages_data[E135_FL_averages_data$clone == "251",]
E135_FL_averages_data$stage <- "E135"
all_stages251_normalised_avg_FL <- rbind(all_stages251_normalised_avg_FL,E135_FL_averages_data)

E135_HL_averages_data <- read.csv("E135/58vsmutantclones_E135_average_normalised_HL.csv")
E135_HL_averages_data <- E135_HL_averages_data[E135_HL_averages_data$clone == "251",]
E135_HL_averages_data$stage <- "E135"
all_stages251_normalised_avg_HL <- rbind(all_stages251_normalised_avg_HL,E135_HL_averages_data)

rm(E135_normalised, E135_FL_averages_data,E135_HL_averages_data)

#E140
E140_normalised <- read.csv("E140/58vsmutantclones_E140_normalised.csv")
E140_normalised <- E140_normalised[E140_normalised$clone == "251",]
E140_normalised$stage <- "E140"
E140_normalised <-E140_normalised %>% rename(EYFP. = EYFP)
all_stages251_normalised <- rbind(all_stages251_normalised,E140_normalised)

E140_FL_averages_data <- read.csv("E140/58vsmutantclones_E140_average_normalised_FL.csv")
E140_FL_averages_data <- E140_FL_averages_data[E140_FL_averages_data$clone == "251",]
E140_FL_averages_data$stage <- "E140"
E140_FL_averages_data <-E140_FL_averages_data %>% rename(EYFP. = EYFP)
all_stages251_normalised_avg_FL <- rbind(all_stages251_normalised_avg_FL,E140_FL_averages_data)

E140_HL_averages_data <- read.csv("E140/58vsmutantclones_E140_average_normalised_HL.csv")
E140_HL_averages_data <- E140_HL_averages_data[E140_HL_averages_data$clone == "251",]
E140_HL_averages_data <-E140_HL_averages_data %>% rename(EYFP. = EYFP)
E140_HL_averages_data$stage <- "E140"
all_stages251_normalised_avg_HL <- rbind(all_stages251_normalised_avg_HL,E140_HL_averages_data)

rm(E140_normalised, E140_FL_averages_data,E140_HL_averages_data)

#E145
E145_normalised <- read.csv("E145_removing_dataset29032023/58vsmutantclones_E145_normalised_without_dataset29032023.csv")
E145_normalised <- E145_normalised[E145_normalised$clone == "251",]
E145_normalised$stage <- "E145"
all_stages251_normalised <- rbind(all_stages251_normalised,E145_normalised)

E145_FL_averages_data <- read.csv("E145_removing_dataset29032023/58vsmutantclones_E145_average_normalised_FL_without_dataset29032023.csv")
E145_FL_averages_data <- E145_FL_averages_data[E145_FL_averages_data$clone == "251",]
E145_FL_averages_data$stage <- "E145"
all_stages251_normalised_avg_FL <- rbind(all_stages251_normalised_avg_FL,E145_FL_averages_data)

E145_HL_averages_data <- read.csv("E145_removing_dataset29032023/58vsmutantclones_E145_average_normalised_HL_without_dataset29032023.csv")
E145_HL_averages_data <- E145_HL_averages_data[E145_HL_averages_data$clone == "251",]
E145_HL_averages_data$stage <- "E145"
all_stages251_normalised_avg_HL <- rbind(all_stages251_normalised_avg_HL,E145_HL_averages_data)

rm(E145_normalised, E145_FL_averages_data,E145_HL_averages_data)

#save summary tables
write.csv(all_stages251_normalised,"all_stages251_normalised.csv",row.names = FALSE)
write.csv(all_stages251_normalised_avg_FL,"all_stages251_normalised_avg_FL.csv",row.names = FALSE)
write.csv(all_stages251_normalised_avg_HL,"all_stages251_normalised_avg_HL.csv",row.names = FALSE)

### 251 and 58 selected stages ####
all_stages58_normalised <- read.csv("/Users/roucogar/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/1_FACS-sorter/Flowjo_analysis_nov2023/58_all_stages/all_stages58_normalised.csv")
exclude_stages <- c('E110','E125','E130','E140')
all_stages58_normalised3 <- subset(all_stages58_normalised,!(stage %in% exclude_stages))
all_stages58_normalised3

all_stages251_normalised <- read.csv("/Users/roucogar/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/1_FACS-sorter/Flowjo_analysis_nov2023/251_all_stages/all_stages251_normalised.csv")
exclude_stages2 <- c('E110','E140')
all_stages251_normalised2 <- subset(all_stages251_normalised,!(stage %in% exclude_stages2))
all_stages251_normalised2

all_stages58and251 <- rbind(all_stages58_normalised3,all_stages251_normalised2)
all_stages58and251$clone <- factor(all_stages58and251$clone, levels = c('58', '251'))

#DP
p1 <-ggplot(all_stages58and251, aes(x = clone, y = DP.., fill = clone)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  geom_jitter(shape = 16, position = position_jitter(0.06), color = "black") +
  ylim(0, 50) +
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 11, face = "bold")) +
  stat_summary(fun = mean, geom = "crossbar", linetype = "dashed", linewidth = 0.2, width = 0.5) +
  facet_grid(tissue ~ stage) +
  labs(y = "Percentage of Double Positive", x = "Stage") +
  scale_fill_manual(values = c('58' = '#0CB702', '251' = '#00A9FF'))

my_comparisons1 <- list(c("58","251"))

p1 + stat_compare_means(comparisons = my_comparisons1, method = "t.test") +
  ggtitle("Percentage of Normalised Double Positive cells across different tissues and genotypes, t-test")

ggsave("Percentage_of_Norm_Double_Positive_in58vs251_cells_selected_stages_across_different_stages_and_limbtype.pdf", width = 12, height = 8)

#EYFP
p2 <-ggplot(all_stages58and251, aes(x = clone, y = EYFP., fill = clone)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  geom_jitter(shape = 16, position = position_jitter(0.06), color = "black") +
  ylim(0, 50) +
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 11, face = "bold")) +
  stat_summary(fun = mean, geom = "crossbar", linetype = "dashed", linewidth = 0.2, width = 0.5) +
  facet_grid(tissue ~ stage) +
  labs(y = "Percentage of EYFP", x = "Stage") +
  scale_fill_manual(values = c('58' = '#0CB702', '251' = '#00A9FF'))

my_comparisons1 <- list(c("58","251"))

p2 + stat_compare_means(comparisons = my_comparisons1, method = "t.test") +
  ggtitle("Percentage of Normalised EYFP cells across different tissues and genotypes, t-test")

ggsave("Percentage_of_Norm_EYFP_in58vs251_cells_selected_stages_across_different_stages_and_limbtype.pdf", width = 12, height = 8)

#dsmch
p3 <-ggplot(all_stages58and251, aes(x = clone, y = dsmch., fill = clone)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  geom_jitter(shape = 16, position = position_jitter(0.06), color = "black") +
  ylim(0, 40) +
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 11, face = "bold")) +
  stat_summary(fun = mean, geom = "crossbar", linetype = "dashed", linewidth = 0.2, width = 0.5) +
  facet_grid(tissue ~ stage) +
  labs(y = "Percentage of dsmCherry", x = "Stage") +
  scale_fill_manual(values = c('58' = '#0CB702', '251' = '#00A9FF'))

my_comparisons1 <- list(c("58","251"))

p3 + stat_compare_means(comparisons = my_comparisons1, method = "t.test") +
  ggtitle("Percentage of Normalised dsmCherry cells across different tissues and genotypes, t-test")

ggsave("Percentage_of_Norm_dsmcherry_in58vs251_cells_selected_stages_across_different_stages_and_limbtype.pdf", width = 12, height = 8)

#LowEYFP
p4 <-ggplot(all_stages58and251, aes(x = clone, y = Low_EYFP, fill = clone)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  geom_jitter(shape = 16, position = position_jitter(0.06), color = "black") +
  ylim(0, 20) +
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 11, face = "bold")) +
  stat_summary(fun = mean, geom = "crossbar", linetype = "dashed", linewidth = 0.2, width = 0.5) +
  facet_grid(tissue ~ stage) +
  labs(y = "Percentage of LowEYFP", x = "Stage") +
  scale_fill_manual(values = c('58' = '#0CB702', '251' = '#00A9FF'))

my_comparisons1 <- list(c("58","251"))

p4 + stat_compare_means(comparisons = my_comparisons1, method = "t.test") +
  ggtitle("Percentage of Normalised LowEYFP cells across different tissues and genotypes, t-test")

ggsave("Percentage_of_Norm_LowEYFP_in58vs251_cells_selected_stages__across_different_stages_and_limbtype.pdf", width = 12, height = 8)

#Interm
p5 <-ggplot(all_stages58and251, aes(x = clone, y = Interm, fill = clone)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  geom_jitter(shape = 16, position = position_jitter(0.06), color = "black") +
  ylim(0, 20) +
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 11, face = "bold")) +
  stat_summary(fun = mean, geom = "crossbar", linetype = "dashed", linewidth = 0.2, width = 0.5) +
  facet_grid(tissue ~ stage) +
  labs(y = "Percentage of Interm", x = "Stage") +
  scale_fill_manual(values = c('58' = '#0CB702', '251' = '#00A9FF'))

my_comparisons1 <- list(c("58","251"))

p5 + stat_compare_means(comparisons = my_comparisons1, method = "t.test") +
  ggtitle("Percentage of Normalised Interm cells across different tissues and genotypes, t-test")

ggsave("Percentage_of_Norm_Interm_in58vs251_cells_selected_stages__across_different_stages_and_limbtype.pdf", width = 12, height = 8)

#NEG
p6 <-ggplot(all_stages58and251, aes(x = clone, y = NEG., fill = clone)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  geom_jitter(shape = 16, position = position_jitter(0.06), color = "black") +
  ylim(0, 110) +
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 11, face = "bold")) +
  stat_summary(fun = mean, geom = "crossbar", linetype = "dashed", linewidth = 0.2, width = 0.5) +
  facet_grid(tissue ~ stage) +
  labs(y = "Percentage of NEG", x = "Stage") +
  scale_fill_manual(values = c('58' = '#0CB702', '251' = '#00A9FF'))

my_comparisons1 <- list(c("58","251"))

p6 + stat_compare_means(comparisons = my_comparisons1, method = "t.test") +
  ggtitle("Percentage of Normalised NEG cells across different tissues and genotypes, t-test")

ggsave("Percentage_of_Norm_NEG_in58vs251_cells_selected_stages_across_different_stages_and_limbtype.pdf", width = 12, height = 8)

### 242 all stages ####
#E105
E105_normalised <- read.csv("E105/58vsmutantclones_E105_normalised.csv")
E105_normalised <- E105_normalised[E105_normalised$clone == "242",]
E105_normalised$stage <- "E105"
all_stages242_normalised <- E105_normalised

E105_FL_averages_data <- read.csv("E105/58vsmutantclones_E105_average_normalised_FL.csv")
E105_FL_averages_data <- E105_FL_averages_data[E105_FL_averages_data$clone == "242",]
E105_FL_averages_data$stage <- "E105"
all_stages242_normalised_avg_FL <- E105_FL_averages_data

E105_HL_averages_data <- read.csv("E105/58vsmutantclones_E105_average_normalised_HL.csv")
E105_HL_averages_data <- E105_HL_averages_data[E105_HL_averages_data$clone == "242",]
E105_HL_averages_data$stage <- "E105"
all_stages242_normalised_avg_HL <- E105_HL_averages_data

rm(E105_normalised, E105_FL_averages_data,E105_HL_averages_data)

#E110
E110_normalised <- read.csv("E110/58vsmutantclones_E110_normalised.csv")
E110_normalised <- E110_normalised[E110_normalised$clone == "242",]
E110_normalised$stage <- "E110"
all_stages242_normalised <- rbind(all_stages242_normalised,E110_normalised)

E110_FL_averages_data <- read.csv("E110/58vsmutantclones_E110_average_normalised_FL.csv")
E110_FL_averages_data <- E110_FL_averages_data[E110_FL_averages_data$clone == "242",]
E110_FL_averages_data$stage <- "E110"
all_stages242_normalised_avg_FL <- rbind(all_stages242_normalised_avg_FL,E110_FL_averages_data)

E110_HL_averages_data <- read.csv("E110/58vsmutantclones_E110_average_normalised_HL.csv")
E110_HL_averages_data <- E110_HL_averages_data[E110_HL_averages_data$clone == "242",]
E110_HL_averages_data$stage <- "E110"
all_stages242_normalised_avg_HL <- rbind(all_stages242_normalised_avg_HL,E110_HL_averages_data)

rm(E110_normalised, E110_FL_averages_data,E110_HL_averages_data)

#E115
E115_normalised <- read.csv("E115/58vsmutantclones_E115_normalised.csv")
E115_normalised <- E115_normalised[E115_normalised$clone == "242",]
E115_normalised$stage <- "E115"
all_stages242_normalised <- rbind(all_stages242_normalised,E115_normalised)

E115_FL_averages_data <- read.csv("E115/58vsmutantclones_E115_average_normalised_FL.csv")
E115_FL_averages_data <- E115_FL_averages_data[E115_FL_averages_data$clone == "242",]
E115_FL_averages_data$stage <- "E115"
all_stages242_normalised_avg_FL <- rbind(all_stages242_normalised_avg_FL,E115_FL_averages_data)

E115_HL_averages_data <- read.csv("E115/58vsmutantclones_E115_average_normalised_HL.csv")
E115_HL_averages_data <- E115_HL_averages_data[E115_HL_averages_data$clone == "242",]
E115_HL_averages_data$stage <- "E115"
all_stages242_normalised_avg_HL <- rbind(all_stages242_normalised_avg_HL,E115_HL_averages_data)

rm(E115_normalised, E115_FL_averages_data,E115_HL_averages_data)

#E135
E135_normalised <- read.csv("E135/58vsmutantclones_E135_normalised.csv")
E135_normalised <- E135_normalised[E135_normalised$clone == "242",]
E135_normalised$stage <- "E135"
all_stages242_normalised <- rbind(all_stages242_normalised,E135_normalised)

E135_FL_averages_data <- read.csv("E135/58vsmutantclones_E135_average_normalised_FL.csv")
E135_FL_averages_data <- E135_FL_averages_data[E135_FL_averages_data$clone == "242",]
E135_FL_averages_data$stage <- "E135"
all_stages242_normalised_avg_FL <- rbind(all_stages242_normalised_avg_FL,E135_FL_averages_data)

E135_HL_averages_data <- read.csv("E135/58vsmutantclones_E135_average_normalised_HL.csv")
E135_HL_averages_data <- E135_HL_averages_data[E135_HL_averages_data$clone == "242",]
E135_HL_averages_data$stage <- "E135"
all_stages242_normalised_avg_HL <- rbind(all_stages242_normalised_avg_HL,E135_HL_averages_data)

rm(E135_normalised, E135_FL_averages_data,E135_HL_averages_data)

#save summary tables
write.csv(all_stages242_normalised,"all_stages242_normalised.csv",row.names = FALSE)
write.csv(all_stages242_normalised_avg_FL,"all_stages242_normalised_avg_FL.csv",row.names = FALSE)
write.csv(all_stages242_normalised_avg_HL,"all_stages242_normalised_avg_HL.csv",row.names = FALSE)

### 242 and 58 selected stages ####
all_stages58_normalised <- read.csv("/Users/roucogar/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/1_FACS-sorter/Flowjo_analysis_nov2023/58_all_stages/all_stages58_normalised.csv")
exclude_stages <- c('E110','E125','E130','E140','E145')
all_stages58_normalised3 <- subset(all_stages58_normalised,!(stage %in% exclude_stages))
all_stages58_normalised3

all_stages242_normalised <- read.csv("/Users/roucogar/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/1_FACS-sorter/Flowjo_analysis_nov2023/242_all_stages/all_stages242_normalised.csv")
exclude_stages2 <- c('E110')
all_stages242_normalised2 <- subset(all_stages242_normalised,!(stage %in% exclude_stages2))
all_stages242_normalised2

all_stages58and242 <- rbind(all_stages58_normalised3,all_stages242_normalised2)
all_stages58and242$clone <- factor(all_stages58and242$clone, levels = c('58', '242'))

#DP
p1 <-ggplot(all_stages58and242, aes(x = clone, y = DP.., fill = clone)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  geom_jitter(shape = 16, position = position_jitter(0.06), color = "black") +
  ylim(0, 50) +
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 11, face = "bold")) +
  stat_summary(fun = mean, geom = "crossbar", linetype = "dashed", linewidth = 0.2, width = 0.5) +
  facet_grid(tissue ~ stage) +
  labs(y = "Percentage of Double Positive", x = "Stage") +
  scale_fill_manual(values = c('58' = '#0CB702', '242' = '#00A9FF'))

my_comparisons1 <- list(c("58","242"))

p1 + stat_compare_means(comparisons = my_comparisons1, method = "t.test") +
  ggtitle("Percentage of Normalised Double Positive cells across different tissues and genotypes, t-test")

ggsave("Percentage_of_Norm_Double_Positive_in58vs242_cells_selected_stages_across_different_stages_and_limbtype.pdf", width = 12, height = 8)

#EYFP
p2 <-ggplot(all_stages58and242, aes(x = clone, y = EYFP., fill = clone)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  geom_jitter(shape = 16, position = position_jitter(0.06), color = "black") +
  ylim(0, 50) +
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 11, face = "bold")) +
  stat_summary(fun = mean, geom = "crossbar", linetype = "dashed", linewidth = 0.2, width = 0.5) +
  facet_grid(tissue ~ stage) +
  labs(y = "Percentage of EYFP", x = "Stage") +
  scale_fill_manual(values = c('58' = '#0CB702', '242' = '#00A9FF'))

my_comparisons1 <- list(c("58","242"))

p2 + stat_compare_means(comparisons = my_comparisons1, method = "t.test") +
  ggtitle("Percentage of Normalised EYFP cells across different tissues and genotypes, t-test")

ggsave("Percentage_of_Norm_EYFP_in58vs242_cells_selected_stages_across_different_stages_and_limbtype.pdf", width = 12, height = 8)

#dsmch
p3 <-ggplot(all_stages58and242, aes(x = clone, y = dsmch., fill = clone)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  geom_jitter(shape = 16, position = position_jitter(0.06), color = "black") +
  ylim(0, 40) +
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 11, face = "bold")) +
  stat_summary(fun = mean, geom = "crossbar", linetype = "dashed", linewidth = 0.2, width = 0.5) +
  facet_grid(tissue ~ stage) +
  labs(y = "Percentage of dsmCherry", x = "Stage") +
  scale_fill_manual(values = c('58' = '#0CB702', '242' = '#00A9FF'))

my_comparisons1 <- list(c("58","242"))

p3 + stat_compare_means(comparisons = my_comparisons1, method = "t.test") +
  ggtitle("Percentage of Normalised dsmCherry cells across different tissues and genotypes, t-test")

ggsave("Percentage_of_Norm_dsmcherry_in58vs242_cells_selected_stages_across_different_stages_and_limbtype.pdf", width = 12, height = 8)

#LowEYFP
p4 <-ggplot(all_stages58and242, aes(x = clone, y = Low_EYFP, fill = clone)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  geom_jitter(shape = 16, position = position_jitter(0.06), color = "black") +
  ylim(0, 20) +
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 11, face = "bold")) +
  stat_summary(fun = mean, geom = "crossbar", linetype = "dashed", linewidth = 0.2, width = 0.5) +
  facet_grid(tissue ~ stage) +
  labs(y = "Percentage of LowEYFP", x = "Stage") +
  scale_fill_manual(values = c('58' = '#0CB702', '242' = '#00A9FF'))

my_comparisons1 <- list(c("58","242"))

p4 + stat_compare_means(comparisons = my_comparisons1, method = "t.test") +
  ggtitle("Percentage of Normalised LowEYFP cells across different tissues and genotypes, t-test")

ggsave("Percentage_of_Norm_LowEYFP_in58vs242_cells_selected_stages__across_different_stages_and_limbtype.pdf", width = 12, height = 8)

#Interm
p5 <-ggplot(all_stages58and242, aes(x = clone, y = Interm, fill = clone)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  geom_jitter(shape = 16, position = position_jitter(0.06), color = "black") +
  ylim(0, 20) +
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 11, face = "bold")) +
  stat_summary(fun = mean, geom = "crossbar", linetype = "dashed", linewidth = 0.2, width = 0.5) +
  facet_grid(tissue ~ stage) +
  labs(y = "Percentage of Interm", x = "Stage") +
  scale_fill_manual(values = c('58' = '#0CB702', '242' = '#00A9FF'))

my_comparisons1 <- list(c("58","242"))

p5 + stat_compare_means(comparisons = my_comparisons1, method = "t.test") +
  ggtitle("Percentage of Normalised Interm cells across different tissues and genotypes, t-test")

ggsave("Percentage_of_Norm_Interm_in58vs242_cells_selected_stages__across_different_stages_and_limbtype.pdf", width = 12, height = 8)

#NEG
p6 <-ggplot(all_stages58and242, aes(x = clone, y = NEG., fill = clone)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  geom_jitter(shape = 16, position = position_jitter(0.06), color = "black") +
  ylim(0, 110) +
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 11, face = "bold")) +
  stat_summary(fun = mean, geom = "crossbar", linetype = "dashed", linewidth = 0.2, width = 0.5) +
  facet_grid(tissue ~ stage) +
  labs(y = "Percentage of NEG", x = "Stage") +
  scale_fill_manual(values = c('58' = '#0CB702', '242' = '#00A9FF'))

my_comparisons1 <- list(c("58","242"))

p6 + stat_compare_means(comparisons = my_comparisons1, method = "t.test") +
  ggtitle("Percentage of Normalised NEG cells across different tissues and genotypes, t-test")

ggsave("Percentage_of_Norm_NEG_in58vs242_cells_selected_stages_across_different_stages_and_limbtype.pdf", width = 12, height = 8)

### 198 all stages ####
#E115
E115_normalised <- read.csv("E115/58vsmutantclones_E115_normalised.csv")
E115_normalised <- E115_normalised[E115_normalised$clone == "198",]
E115_normalised$stage <- "E115"
all_stages198_normalised <- E115_normalised

E115_FL_averages_data <- read.csv("E115/58vsmutantclones_E115_average_normalised_FL.csv")
E115_FL_averages_data <- E115_FL_averages_data[E115_FL_averages_data$clone == "198",]
E115_FL_averages_data$stage <- "E115"
all_stages198_normalised_avg_FL <- E115_FL_averages_data

E115_HL_averages_data <- read.csv("E115/58vsmutantclones_E115_average_normalised_HL.csv")
E115_HL_averages_data <- E115_HL_averages_data[E115_HL_averages_data$clone == "198",]
E115_HL_averages_data$stage <- "E115"
all_stages198_normalised_avg_HL <- E115_HL_averages_data

rm(E115_normalised, E115_FL_averages_data,E115_HL_averages_data)

#E135
E135_normalised <- read.csv("E135/58vsmutantclones_E135_normalised.csv")
E135_normalised <- E135_normalised[E135_normalised$clone == "198",]
E135_normalised$stage <- "E135"
all_stages198_normalised <- rbind(all_stages198_normalised,E135_normalised)

E135_FL_averages_data <- read.csv("E135/58vsmutantclones_E135_average_normalised_FL.csv")
E135_FL_averages_data <- E135_FL_averages_data[E135_FL_averages_data$clone == "198",]
E135_FL_averages_data$stage <- "E135"
all_stages198_normalised_avg_FL <- rbind(all_stages198_normalised_avg_FL,E135_FL_averages_data)

E135_HL_averages_data <- read.csv("E135/58vsmutantclones_E135_average_normalised_HL.csv")
E135_HL_averages_data <- E135_HL_averages_data[E135_HL_averages_data$clone == "198",]
E135_HL_averages_data$stage <- "E135"
all_stages198_normalised_avg_HL <- rbind(all_stages198_normalised_avg_HL,E135_HL_averages_data)

rm(E135_normalised, E135_FL_averages_data,E135_HL_averages_data)

#E140
E140_normalised <- read.csv("E140/58vsmutantclones_E140_normalised.csv")
E140_normalised <- E140_normalised[E140_normalised$clone == "198",]
E140_normalised$stage <- "E140"
E140_normalised <-E140_normalised %>% rename(EYFP. = EYFP)
all_stages198_normalised <- rbind(all_stages198_normalised,E140_normalised)

E140_FL_averages_data <- read.csv("E140/58vsmutantclones_E140_average_normalised_FL.csv")
E140_FL_averages_data <- E140_FL_averages_data[E140_FL_averages_data$clone == "198",]
E140_FL_averages_data$stage <- "E140"
E140_FL_averages_data <-E140_FL_averages_data %>% rename(EYFP. = EYFP)
all_stages198_normalised_avg_FL <- rbind(all_stages198_normalised_avg_FL,E140_FL_averages_data)

E140_HL_averages_data <- read.csv("E140/58vsmutantclones_E140_average_normalised_HL.csv")
E140_HL_averages_data <- E140_HL_averages_data[E140_HL_averages_data$clone == "198",]
E140_HL_averages_data <-E140_HL_averages_data %>% rename(EYFP. = EYFP)
E140_HL_averages_data$stage <- "E140"
all_stages198_normalised_avg_HL <- rbind(all_stages198_normalised_avg_HL,E140_HL_averages_data)

rm(E140_normalised, E140_FL_averages_data,E140_HL_averages_data)

#E145
E145_normalised <- read.csv("E145_removing_dataset29032023/58vsmutantclones_E145_normalised_without_dataset29032023.csv")
E145_normalised <- E145_normalised[E145_normalised$clone == "198",]
E145_normalised$stage <- "E145"
all_stages198_normalised <- rbind(all_stages198_normalised,E145_normalised)

E145_FL_averages_data <- read.csv("E145_removing_dataset29032023/58vsmutantclones_E145_average_normalised_FL_without_dataset29032023.csv")
E145_FL_averages_data <- E145_FL_averages_data[E145_FL_averages_data$clone == "198",]
E145_FL_averages_data$stage <- "E145"
all_stages198_normalised_avg_FL <- rbind(all_stages198_normalised_avg_FL,E145_FL_averages_data)

E145_HL_averages_data <- read.csv("E145_removing_dataset29032023/58vsmutantclones_E145_average_normalised_HL_without_dataset29032023.csv")
E145_HL_averages_data <- E145_HL_averages_data[E145_HL_averages_data$clone == "198",]
E145_HL_averages_data$stage <- "E145"
all_stages198_normalised_avg_HL <- rbind(all_stages198_normalised_avg_HL,E145_HL_averages_data)

rm(E145_normalised, E145_FL_averages_data,E145_HL_averages_data)

#save summary tables
write.csv(all_stages198_normalised,"all_stages198_normalised.csv",row.names = FALSE)
write.csv(all_stages198_normalised_avg_FL,"all_stages198_normalised_avg_FL.csv",row.names = FALSE)
write.csv(all_stages198_normalised_avg_HL,"all_stages198_normalised_avg_HL.csv",row.names = FALSE)

### 198 and 58 selected stages ####
all_stages58_normalised <- read.csv("/Users/roucogar/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/1_FACS-sorter/Flowjo_analysis_nov2023/58_all_stages/all_stages58_normalised.csv")
exclude_stages <- c('E105','E110','E125','E130','E140')
all_stages58_normalised3 <- subset(all_stages58_normalised,!(stage %in% exclude_stages))
all_stages58_normalised3

all_stages198_normalised <- read.csv("/Users/roucogar/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/1_FACS-sorter/Flowjo_analysis_nov2023/198_all_stages/all_stages198_normalised.csv")
exclude_stages2 <- c('E140')
all_stages198_normalised2 <- subset(all_stages198_normalised,!(stage %in% exclude_stages2))
all_stages198_normalised2

all_stages58and198 <- rbind(all_stages58_normalised3,all_stages198_normalised2)
all_stages58and198$clone <- factor(all_stages58and198$clone, levels = c('58', '198'))

#DP
p1 <-ggplot(all_stages58and198, aes(x = clone, y = DP.., fill = clone)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  geom_jitter(shape = 16, position = position_jitter(0.06), color = "black") +
  ylim(0, 50) +
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 11, face = "bold")) +
  stat_summary(fun = mean, geom = "crossbar", linetype = "dashed", linewidth = 0.2, width = 0.5) +
  facet_grid(tissue ~ stage) +
  labs(y = "Percentage of Double Positive", x = "Stage") +
  scale_fill_manual(values = c('58' = '#0CB702', '198' = '#00A9FF'))

my_comparisons1 <- list(c("58","198"))

p1 + stat_compare_means(comparisons = my_comparisons1, method = "t.test") +
  ggtitle("Percentage of Normalised Double Positive cells across different tissues and genotypes, t-test")

ggsave("Percentage_of_Norm_Double_Positive_in58vs198_cells_selected_stages_across_different_stages_and_limbtype.pdf", width = 12, height = 8)

#EYFP
p2 <-ggplot(all_stages58and198, aes(x = clone, y = EYFP., fill = clone)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  geom_jitter(shape = 16, position = position_jitter(0.06), color = "black") +
  ylim(0, 50) +
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 11, face = "bold")) +
  stat_summary(fun = mean, geom = "crossbar", linetype = "dashed", linewidth = 0.2, width = 0.5) +
  facet_grid(tissue ~ stage) +
  labs(y = "Percentage of EYFP", x = "Stage") +
  scale_fill_manual(values = c('58' = '#0CB702', '198' = '#00A9FF'))

my_comparisons1 <- list(c("58","198"))

p2 + stat_compare_means(comparisons = my_comparisons1, method = "t.test") +
  ggtitle("Percentage of Normalised EYFP cells across different tissues and genotypes, t-test")

ggsave("Percentage_of_Norm_EYFP_in58vs198_cells_selected_stages_across_different_stages_and_limbtype.pdf", width = 12, height = 8)

#dsmch
p3 <-ggplot(all_stages58and198, aes(x = clone, y = dsmch., fill = clone)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  geom_jitter(shape = 16, position = position_jitter(0.06), color = "black") +
  ylim(0, 40) +
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 11, face = "bold")) +
  stat_summary(fun = mean, geom = "crossbar", linetype = "dashed", linewidth = 0.2, width = 0.5) +
  facet_grid(tissue ~ stage) +
  labs(y = "Percentage of dsmCherry", x = "Stage") +
  scale_fill_manual(values = c('58' = '#0CB702', '198' = '#00A9FF'))

my_comparisons1 <- list(c("58","198"))

p3 + stat_compare_means(comparisons = my_comparisons1, method = "t.test") +
  ggtitle("Percentage of Normalised dsmCherry cells across different tissues and genotypes, t-test")

ggsave("Percentage_of_Norm_dsmcherry_in58vs198_cells_selected_stages_across_different_stages_and_limbtype.pdf", width = 12, height = 8)

#LowEYFP
p4 <-ggplot(all_stages58and198, aes(x = clone, y = Low_EYFP, fill = clone)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  geom_jitter(shape = 16, position = position_jitter(0.06), color = "black") +
  ylim(0, 20) +
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 11, face = "bold")) +
  stat_summary(fun = mean, geom = "crossbar", linetype = "dashed", linewidth = 0.2, width = 0.5) +
  facet_grid(tissue ~ stage) +
  labs(y = "Percentage of LowEYFP", x = "Stage") +
  scale_fill_manual(values = c('58' = '#0CB702', '198' = '#00A9FF'))

my_comparisons1 <- list(c("58","198"))

p4 + stat_compare_means(comparisons = my_comparisons1, method = "t.test") +
  ggtitle("Percentage of Normalised LowEYFP cells across different tissues and genotypes, t-test")

ggsave("Percentage_of_Norm_LowEYFP_in58vs198_cells_selected_stages__across_different_stages_and_limbtype.pdf", width = 12, height = 8)

#Interm
p5 <-ggplot(all_stages58and198, aes(x = clone, y = Interm, fill = clone)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  geom_jitter(shape = 16, position = position_jitter(0.06), color = "black") +
  ylim(0, 20) +
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 11, face = "bold")) +
  stat_summary(fun = mean, geom = "crossbar", linetype = "dashed", linewidth = 0.2, width = 0.5) +
  facet_grid(tissue ~ stage) +
  labs(y = "Percentage of Interm", x = "Stage") +
  scale_fill_manual(values = c('58' = '#0CB702', '198' = '#00A9FF'))

my_comparisons1 <- list(c("58","198"))

p5 + stat_compare_means(comparisons = my_comparisons1, method = "t.test") +
  ggtitle("Percentage of Normalised Interm cells across different tissues and genotypes, t-test")

ggsave("Percentage_of_Norm_Interm_in58vs198_cells_selected_stages__across_different_stages_and_limbtype.pdf", width = 12, height = 8)

#NEG
p6 <-ggplot(all_stages58and198, aes(x = clone, y = NEG., fill = clone)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "black") +
  geom_jitter(shape = 16, position = position_jitter(0.06), color = "black") +
  ylim(0, 110) +
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 11, face = "bold")) +
  stat_summary(fun = mean, geom = "crossbar", linetype = "dashed", linewidth = 0.2, width = 0.5) +
  facet_grid(tissue ~ stage) +
  labs(y = "Percentage of NEG", x = "Stage") +
  scale_fill_manual(values = c('58' = '#0CB702', '198' = '#00A9FF'))

my_comparisons1 <- list(c("58","198"))

p6 + stat_compare_means(comparisons = my_comparisons1, method = "t.test") +
  ggtitle("Percentage of Normalised NEG cells across different tissues and genotypes, t-test")

ggsave("Percentage_of_Norm_NEG_in58vs198_cells_selected_stages_across_different_stages_and_limbtype.pdf", width = 12, height = 8)




