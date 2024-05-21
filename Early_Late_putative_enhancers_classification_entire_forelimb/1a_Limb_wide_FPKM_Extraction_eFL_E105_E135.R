library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyverse)

#### Loading eFL samples ####

DESeq2Analysis_PerGenoFL <- read.delim("/Users/roucogar/Library/CloudStorage/Dropbox/Trajectories_Rouco_et_al/Limb_wide_early-late/outputs/mergedTables/AllCufflinks_Simplified_norm.txt")
head(DESeq2Analysis_PerGenoFL)

head(DESeq2Analysis_PerGenoFL[,c(2,4:7)])

#### Loading list of genes to filter ####

gene_list <- read_lines("/Users/roucogar/Library/CloudStorage/Dropbox/Trajectories_Rouco_et_al/Limb_wide_early-late/Curated_Gene_list.txt")

df2 <- filter(DESeq2Analysis_PerGenoFL, gene_short_name %in% gene_list)

df_eFL2 <-df2[,c(2,4:7)]

df_eFL2 <-melt(df_eFL2)

df_eFL2$stage <- df_eFL2$variable
df_eFL2$stage <- gsub("[0-9]$","", df_eFL2$stage)
df_eFL2$stage <- gsub("_rep","", df_eFL2$stage)

df_eFL2$avg_value <- ave(df_eFL2$value,df_eFL2$gene_short_name,df_eFL2$stage)

##### keeping average expression #####
GOI_df2 <- df_eFL2[,c(1,4:5)]
head(GOI_df2)

unique_GOI_df2 <- GOI_df2 %>%  distinct(gene_short_name, stage, .keep_all = TRUE)

GOI_df3 <- dcast (unique_GOI_df2, gene_short_name ~ stage, value.var = "avg_value")

GOI_df3$ratio <- GOI_df3$FPKM_FL_E105/GOI_df3$FPKM_FL_E135
head(GOI_df3)

write.csv(GOI_df3, file = "/Users/roucogar/Library/CloudStorage/Dropbox/Trajectories_Rouco_et_al/Limb_wide_early-late/Curated_Gene_list_avg_ratio.csv", row.names = FALSE)


#### filtering genes of interest ###
## criteria 1 fpkm >5 at E10.5 and E13.5
summary(GOI_df3)

filtered_GOI_1 <- GOI_df3[GOI_df3$FPKM_FL_E105 > 5,]
head(filtered_GOI_1)
summary(filtered_GOI_1)

filtered_GOI_2 <- filtered_GOI_1[filtered_GOI_1$FPKM_FL_E135 > 5,]
head(filtered_GOI_2)
summary(filtered_GOI_2)

## criteria 2 ratio < 3x & > 0.333333

filtered_GOI_3 <- filtered_GOI_2[filtered_GOI_2$ratio < 3,]
head(filtered_GOI_3)
summary(filtered_GOI_3)

filtered_GOI <- filtered_GOI_3[filtered_GOI_3$ratio > 0.333333,]
head(filtered_GOI)
summary(filtered_GOI)

write.csv(filtered_GOI, file = "/Users/roucogar/Library/CloudStorage/Dropbox/Trajectories_Rouco_et_al/Limb_wide_early-late/Filtered_Curated_Gene_list_avg_ratio.csv", row.names = FALSE)


