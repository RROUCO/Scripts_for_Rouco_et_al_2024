library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(gplots)
library(RColorBrewer)

setwd("Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/7_RNA-seq/Final_analysis/R_analysis_Heatmap/")

#### Loading FL samples ####
DESeq2Analysis_PerGenoFL <- read.delim("~/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/7_RNA-seq/Final_analysis/RNA-seq_58_Forelimb/outputs/mergedTables/AllCufflinks_Simplified_norm.txt")
head(DESeq2Analysis_PerGenoFL)

#start editing the dataframe to keep data columns and rows that you're interested in
#To work with FPKM

#Forelimb
head(DESeq2Analysis_PerGenoFL[,c(2,4:43)])
df_FL <-DESeq2Analysis_PerGenoFL[,c(2,4:43)]

df_FL <-melt(df_FL)
df_FL$condition <- gsub("[0-9]*$","", df_FL$variable)
df_FL$condition <- gsub("_RR$","", df_FL$condition)

df_FL$facs_stage <- gsub("FPKM_","", df_FL$condition)
df_FL$facs_stage <- gsub("[0-9]*$","", df_FL$facs_stage)

df_FL$condition <- gsub("[0-9]$","", df_FL$condition)
df_FL$condition <- gsub("FPKM_","", df_FL$condition)

df_FL$stage <- gsub("FL58_","", df_FL$condition)
df_FL$stage <- gsub("_[A-z]*","", df_FL$stage)

df_FL$facs_cell <- gsub("FL58_E","", df_FL$condition)
df_FL$facs_cell <- gsub("[0-9]*_","", df_FL$facs_cell)

df_FL$avg_value <- ave(df_FL$value,df_FL$gene_short_name,df_FL$facs_stage)

df_FL2<- df_FL[df_FL$stage != "E110",]

#### Forelimb heatmap for figure 3D ####
#Keep only for figure 3D
df_FL_facs <- df_FL2[df_FL2$facs_cell %in% c("DP","EYFP","NEG"),]

#genes for figure 3D
#gene_list2 <- c("Shox2","CRE","dsmCherry","EYFP","Hoxd13","Sall3","Sall4","Meis1",
#                "Meis2","Hoxa11","Dcn","Lum","Kera","Osr1","Acan","Col3a1","Etv4","Shh",
#                "Irx1","Col9a3","Irx3","Irx5","Msx1")

gene_list2 <- c("Shox2","CRE","EYFP","Hoxd13","Hoxa11","Sall3","Sall4","Irx3","Msx1","Etv4","Shh",
                "Lum","Osr1","Col3a1","Acan","Col9a3","Irx1")

GOI_df2 <- filter(df_FL_facs, gene_short_name %in% gene_list2)
GOI_df2 <- data.frame(GOI_df2)

heatmap_data2 <- GOI_df2[,c(1,5,8)]
heatmap_data2 <- heatmap_data2 %>%  distinct(gene_short_name, facs_stage, .keep_all = TRUE)
heatmap_data2 <- reshape(heatmap_data2, idvar = "gene_short_name", timevar = "facs_stage", direction = "wide")
rownames(heatmap_data2) <- heatmap_data2$gene_short_name
heatmap_data2 <- heatmap_data2[, -1]

# Normalize the data (z-score normalization)
#You may want to z-score normalize your FPKM values to ensure that the heatmap represents relative expression levels across samples.
#This is only to have an idea to the values that are going to be plot in the heatmap
#there is no need to feed this dataframe to the heatmap2 function
normalized_data2 <- as.data.frame(t(apply(heatmap_data2, 1, scale)))
colnames(normalized_data2) <- colnames(heatmap_data2)

# Define the custom order of columns
custom_order3 <- c("avg_value.FL58_E105_DP","avg_value.FL58_E115_DP","avg_value.FL58_E125_DP","avg_value.FL58_E135_DP",
                   "avg_value.FL58_E105_EYFP","avg_value.FL58_E115_EYFP","avg_value.FL58_E125_EYFP","avg_value.FL58_E135_EYFP",
                   "avg_value.FL58_E105_NEG","avg_value.FL58_E115_NEG","avg_value.FL58_E125_NEG","avg_value.FL58_E135_NEG")

# Reorder the numeric matrix based on the custom order
heatmap_data2 <- heatmap_data2[, match(custom_order3, colnames(heatmap_data2))]

#Define the custom order of rows
custom_order4 <- c("Shox2","CRE","EYFP","Hoxd13","Hoxa11","Sall3","Sall4","Irx3","Msx1","Etv4","Shh",
                   "Lum","Osr1","Col3a1","Acan","Col9a3","Irx1")

# Reorder rows based on the custom order
heatmap_data2 <- heatmap_data2[match(custom_order4, rownames(heatmap_data2)),]

#color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
color =colorRampPalette(brewer.pal(9, "BuPu"))(100)
#cex=character expansion
# Create the heatmap
heatmap.2(
  as.matrix(heatmap_data2),  # Convert the data to a matrix
  Rowv = FALSE,
  Colv = FALSE,
  col = color,
  scale = "row",
  dendrogram = "none",
  ColSideColors=rep(c("#be1622","#F1D31E","#000000"), each=4),
  trace = "none",
  key = TRUE,
  cexRow = 1,
  cexCol = 1,
  main = "Scale by Row FPKM Heatmap",
  margins = c(10, 10)
)


### Loading HL Samples ####
DESeq2Analysis_PerGeno <- read.delim("~/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/7_RNA-seq/Final_analysis/RNA-seq_58_Hindlimb/outputs/mergedTables/AllCufflinks_Simplified_norm.txt")
head(DESeq2Analysis_PerGeno)

#start editing the dataframe to keep data columns and rows that you're interested in
#To work with FPKM
head(DESeq2Analysis_PerGeno[,c(2,4:37)])
df_HL <-DESeq2Analysis_PerGeno[,c(2,4:37)]

df_HL <-melt(df_HL)
df_HL$condition <- gsub("[0-9]*$","", df_HL$variable)
df_HL$condition <- gsub("_RR$","", df_HL$condition)

df_HL$facs_stage <- gsub("FPKM_","", df_HL$condition)
df_HL$facs_stage <- gsub("[0-9]*$","", df_HL$facs_stage)

df_HL$condition <- gsub("[0-9]$","", df_HL$condition)
df_HL$condition <- gsub("FPKM_","", df_HL$condition)

df_HL$stage <- gsub("HL58_","", df_HL$condition)
df_HL$stage <- gsub("_[A-z]*","", df_HL$stage)

df_HL$facs_cell <- gsub("HL58_E","", df_HL$condition)
df_HL$facs_cell <- gsub("[0-9]*_","", df_HL$facs_cell)

df_HL$avg_value <- ave(df_HL$value,df_HL$gene_short_name,df_HL$facs_stage)

#### Hindlimb Heatmap for Figure 3B ####
#Keep only E11.5 for figure 3B
df_HL_E115 <- df_HL[df_HL$stage == "E115",]

#genes for figure 3B
gene_list1 <- c("Shox2","CRE","dsmCherry","EYFP","Hoxd13","Hoxa11","Prrx1","Msx1",
                "Twist1","Hbb-y","Krt5","Krt14","Wnt6","Cdh5","Myod1","Col2a1")

GOI_df1 <- filter(df_HL_E115, gene_short_name %in% gene_list1)
GOI_df1 <- data.frame(GOI_df1)

heatmap_data <- GOI_df1[,c(1,5,8)]
heatmap_data <- heatmap_data %>%  distinct(gene_short_name, facs_stage, .keep_all = TRUE)
heatmap_data <- reshape(heatmap_data, idvar = "gene_short_name", timevar = "facs_stage", direction = "wide")
rownames(heatmap_data) <- heatmap_data$gene_short_name
heatmap_data <- heatmap_data[, -1]

# Normalize the data (z-score normalization)
#You may want to z-score normalize your FPKM values to ensure that the heatmap represents relative expression levels across samples.
#This is only to have an idea to the values that are going to be plot in the heatmap
#there is no need to feed this dataframe to the heatmap2 function
normalized_data <- as.data.frame(t(apply(heatmap_data, 1, scale)))
colnames(normalized_data) <- colnames(heatmap_data)

# Define the custom order of columns
custom_order <- c("avg_value.HL58_E115_dsmch","avg_value.HL58_E115_LowEYFP","avg_value.HL58_E115_DP",
                  "avg_value.HL58_E115_Interm","avg_value.HL58_E115_EYFP","avg_value.HL58_E115_NEG")

# Reorder the numeric matrix based on the custom order
heatmap_data <- heatmap_data[, match(custom_order, colnames(heatmap_data))]

#Define the custom order of rows
custom_order2 <- c("Shox2","CRE","dsmCherry","EYFP","Hoxd13","Hoxa11","Prrx1",
                   "Twist1","Msx1","Col2a1","Hbb-y","Krt5","Krt14","Wnt6","Cdh5","Myod1")

# Reorder rows based on the custom order
heatmap_data <- heatmap_data[match(custom_order2, rownames(heatmap_data)),]

#color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
color =colorRampPalette(brewer.pal(9, "BuPu"))(100)
#cex=character expansion
# Create the heatmap
heatmap.2(
  as.matrix(heatmap_data),  # Convert the data to a matrix
  Rowv = FALSE,
  Colv = FALSE,
  col = color,
  scale = "row",
  dendrogram = "none",
  trace = "none",
  key = TRUE,
  cexRow = 1,
  cexCol = 1,
  main = "Scale by Row FPKM Heatmap",
  margins = c(10, 10)
)

#### Hindlimb heatmap for figure 3D ####
#Keep only for figure 3D
df_HL_facs <- df_HL[df_HL$facs_cell %in% c("DP","EYFP","NEG"),]

#genes for figure 3D
#gene_list2 <- c("Shox2","CRE","dsmCherry","EYFP","Hoxd13","Sall3","Sall4","Meis1",
#                "Meis2","Hoxa11","Dcn","Lum","Kera","Osr1","Acan","Col3a1","Etv4","Shh",
#                "Irx1","Col9a3","Irx3","Irx5","Msx1")

gene_list2 <- c("Shox2","CRE","EYFP","Hoxd13","Hoxa11","Sall3","Sall4","Irx3","Msx1","Etv4","Shh",
                "Lum","Osr1","Col3a1","Acan","Col9a3","Irx1")

GOI_df2 <- filter(df_HL_facs, gene_short_name %in% gene_list2)
GOI_df2 <- data.frame(GOI_df2)

heatmap_data2 <- GOI_df2[,c(1,5,8)]
heatmap_data2 <- heatmap_data2 %>%  distinct(gene_short_name, facs_stage, .keep_all = TRUE)
heatmap_data2 <- reshape(heatmap_data2, idvar = "gene_short_name", timevar = "facs_stage", direction = "wide")
rownames(heatmap_data2) <- heatmap_data2$gene_short_name
heatmap_data2 <- heatmap_data2[, -1]

# Normalize the data (z-score normalization)
#You may want to z-score normalize your FPKM values to ensure that the heatmap represents relative expression levels across samples.
#This is only to have an idea to the values that are going to be plot in the heatmap
#there is no need to feed this dataframe to the heatmap2 function
normalized_data2 <- as.data.frame(t(apply(heatmap_data2, 1, scale)))
colnames(normalized_data2) <- colnames(heatmap_data2)

# Define the custom order of columns
custom_order3 <- c("avg_value.HL58_E105_DP","avg_value.HL58_E115_DP","avg_value.HL58_E125_DP","avg_value.HL58_E135_DP",
                   "avg_value.HL58_E105_EYFP","avg_value.HL58_E115_EYFP","avg_value.HL58_E125_EYFP","avg_value.HL58_E135_EYFP",
                   "avg_value.HL58_E105_NEG","avg_value.HL58_E115_NEG","avg_value.HL58_E125_NEG","avg_value.HL58_E135_NEG")

# Reorder the numeric matrix based on the custom order
heatmap_data2 <- heatmap_data2[, match(custom_order3, colnames(heatmap_data2))]

#Define the custom order of rows
custom_order4 <- c("Shox2","CRE","EYFP","Hoxd13","Hoxa11","Sall3","Sall4","Irx3","Msx1","Etv4","Shh",
                    "Lum","Osr1","Col3a1","Acan","Col9a3","Irx1")

# Reorder rows based on the custom order
heatmap_data2 <- heatmap_data2[match(custom_order4, rownames(heatmap_data2)),]

#color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
color =colorRampPalette(brewer.pal(9, "BuPu"))(100)
#cex=character expansion
# Create the heatmap
heatmap.2(
  as.matrix(heatmap_data2),  # Convert the data to a matrix
  Rowv = FALSE,
  Colv = FALSE,
  col = color,
  scale = "row",
  dendrogram = "none",
  ColSideColors=rep(c("#be1622","#F1D31E","#000000"), each=4),
  trace = "none",
  key = TRUE,
  cexRow = 1,
  cexCol = 1,
  main = "Scale by Row FPKM Heatmap",
  margins = c(10, 10)
)

### Supplementary table S3 only HL ####
#gene list1
list <- c("Shox2","CRE","dsmCherry","EYFP","Hoxd13","Hoxa11","Sall3","Sall4","Irx3","Msx1","Etv4","Shh","Col2a1",
          "Lum","Osr1","Col3a1","Acan","Col9a3","Irx1", "Prrx1","Twist1","Hbb-y","Krt5","Krt14","Wnt6","Cdh5","Myod1" )

GOI_df <- filter(DESeq2Analysis_PerGeno, gene_short_name %in% list)
GOI_df <- data.frame(GOI_df)
head(GOI_df)

write.csv(GOI_df,"~/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/7_RNA-seq/Final_analysis/R_analysis_Heatmap/Selected_Genes_of_Interest_FPKM.csv", row.names = FALSE)

colnames(GOI_df)
GOI_df2 <- GOI_df[,c(1,2,3,4,5,12,13,22,23,30,31,6,14,15,24,25,32,33,7,8,16,17,26,27,34,35,9,19,20,21,10,11,18,28,29,34,35)]
write.csv(GOI_df,"~/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/7_RNA-seq/Final_analysis/R_analysis_Heatmap/Selected_Genes_of_Interest_FPKM_reorderbycelltype.csv", row.names = FALSE)

