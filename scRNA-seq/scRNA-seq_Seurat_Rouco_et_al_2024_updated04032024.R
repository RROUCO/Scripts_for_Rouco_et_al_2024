library(Seurat)
library(dplyr)
library(ggplot2)
library(DoubletFinder)
library(tidyverse)
library(ggraph)
library(clustree)
library(velocyto.R)
library(SeuratDisk)
library(SeuratWrappers)
library(Nebulosa)
library(scCustomize)
library(viridis)
library(magrittr)
library(scales)
library(monocle3)

setwd("/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/")

#### A.LOADING THE SAMPLES DATA #####
#### A.1 -> HL_E105_rep1 getting the cellranger data and the velocity data ####
#First step load rawdata coming from CellRanger into a Seurat object 
#This Seurat Object includes only cells with more 200 genes detected and the genes present in at least 3 cells
HL_E105_rep1.data <- Read10X(data.dir = "/Volumes/RAQUEL-LAB3/Single_cell_58_data/CellRanger_58/cr_HL_E105_58_REP1/outs/filtered_feature_bc_matrix/")

## since the previous files include some extra labelling attached to the cellname I'm gonna remove this extra tag
#data
colnames(HL_E105_rep1.data)
colnames(HL_E105_rep1.data) <- gsub("-1$","",colnames(HL_E105_rep1.data))
colnames(HL_E105_rep1.data)

#Velocyto.data
HL_E105_rep1.velo <- ReadVelocity(file = "/Volumes/RAQUEL-LAB3/Single_cell_58_data/CellRanger_58/cr_HL_E105_58_REP1/velocyto/cr_HL_E105_58_REP1.loom")

## since the previous files include some extra labelling attached to the cellname I'm gonna remove this extra tag
#velo
colnames(HL_E105_rep1.velo$spliced)
colnames(HL_E105_rep1.velo$spliced) <- gsub("cr_HL_E105_58_REP1:","",colnames(HL_E105_rep1.velo$spliced))
colnames(HL_E105_rep1.velo$unspliced) <- gsub("cr_HL_E105_58_REP1:","",colnames(HL_E105_rep1.velo$unspliced))
colnames(HL_E105_rep1.velo$ambiguous) <- gsub("cr_HL_E105_58_REP1:","",colnames(HL_E105_rep1.velo$ambiguous))
colnames(HL_E105_rep1.velo$spliced)

colnames(HL_E105_rep1.velo$spliced) <- gsub("x$","",colnames(HL_E105_rep1.velo$spliced))
colnames(HL_E105_rep1.velo$unspliced) <- gsub("x$","",colnames(HL_E105_rep1.velo$unspliced))
colnames(HL_E105_rep1.velo$ambiguous) <- gsub("x$","",colnames(HL_E105_rep1.velo$ambiguous))
colnames(HL_E105_rep1.velo$spliced)

#test
A<-colnames(HL_E105_rep1.data)
B<-colnames(HL_E105_rep1.velo$ambiguous)
## Test to confirm if the cell names (Colnames) are equal in .data and .velo ##
#if I use all.equal(A,B) without sorting I got this message "3947 string mismatches" 
all.equal(A,B)
#but if I sort B, before running all.equal(A,B) then the answer is TRUE
B <- sort(B)
all.equal(A,B)
#also I've check the similarities by checkin the items in A that are not in B
# it can be done with %in%
A[!(A %in% B)]
# it can be also done with setdiff
setdiff(A,B)
#items in B that are not in A
B[!(B %in% A)]
setdiff(B,A)

## Checking if there are differences in the genes between the file .data (a) and .velo (b) ##
b <- rownames(HL_E105_rep1.velo$spliced)
a <- rownames(HL_E105_rep1.data)
a <- sort(a)
b <- sort(b)
length(a)
length(b)

## Merging both Matrix in one list to generate the Seurat object ##
HL_E105_rep1 <- HL_E105_rep1.velo
HL_E105_rep1$RNA <- HL_E105_rep1.data

## checking the data inside the Matrix ##
HL_E105_rep1$spliced[c("Pitx1", "Shox2", "Hoxa13","mt-Co2"),1:30]
HL_E105_rep1$RNA[c("Pitx1", "Shox2", "Hoxa13","mt-Co2"),1:30]

#create Seurat object
HL_E105_rep1 <- as.Seurat(x = HL_E105_rep1)
head(HL_E105_rep1@meta.data,5)
tail(HL_E105_rep1@meta.data,5)
HL_E105_rep1

DefaultAssay(object = HL_E105_rep1) <- "RNA"
HL_E105_rep1

rm(a,b,A,B,HL_E105_rep1.data,HL_E105_rep1.velo)

# Rename orig.ident for future analysis #
levels(HL_E105_rep1$orig.ident) <- "HL_E105_rep1"
head(HL_E105_rep1@meta.data,5)
#Rename levels also
levels(HL_E105_rep1)
HL_E105_rep1 <- RenameIdents(HL_E105_rep1, "SeuratProject" = "HL_E105_rep1")
levels(HL_E105_rep1)

#calculate % of reads belonging to Mitochondrial genes
HL_E105_rep1[["percent.mt"]] <- PercentageFeatureSet(HL_E105_rep1, pattern = "^mt-",assay = "RNA")
head(HL_E105_rep1@meta.data,5)
summary(HL_E105_rep1@meta.data)

#count cells in the seurat object
t1 <- table(HL_E105_rep1@meta.data[["orig.ident"]])
capture.output(bind_rows(t1), file=paste0("numbers_total.txt"))
t1

### A.1.2 QUALITY CONTROL 1 (QC_1) ###
FeatureScatter(HL_E105_rep1, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(HL_E105_rep1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
VlnPlot(HL_E105_rep1, features = "nCount_RNA", pt.size=0)
VlnPlot(HL_E105_rep1, features = "nFeature_RNA", pt.size=0)

pdf(file=paste0("HL_58_E105_QC_rep1_comparison_before_filtering.pdf"))
VlnPlot(HL_E105_rep1, features = "nCount_RNA", pt.size=0)+xlab("orig.ident")+ylab("# UMI")+labs(title="# UMI")
VlnPlot(HL_E105_rep1, features = "nFeature_RNA", pt.size=0)+xlab("orig.ident")+ylab("# genes")+labs(title="# genes")
VlnPlot(HL_E105_rep1, features = "percent.mt", pt.size=0)+xlab("orig.ident")+ylab("% mito.reads")+labs(title="% mito.reads")
FeatureScatter(HL_E105_rep1, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size=0.2)+xlab("# UMI")+ylab("# percent.mt")
FeatureScatter(HL_E105_rep1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.2)+xlab("# UMI")+ylab("# genes")
dev.off()

# save as seurat object 
saveRDS(HL_E105_rep1, file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E105_rep1_unfiltered.rds")

### A.1. filtering FINAL --> MT > 0.5 < 5 & nfeatures >200 <7500 ###
HL_E105_rep1 <-readRDS(file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E105_rep1_unfiltered.rds")
HL_E105_rep1 <- subset(HL_E105_rep1, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt >0.5 & percent.mt < 5)
FeatureScatter(HL_E105_rep1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(HL_E105_rep1, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size=0.2)+xlab("# UMI")+ylab("# percent.mt")
table(HL_E105_rep1@meta.data[["orig.ident"]])

pdf(file=paste0("HL_58_E105_QC_rep1_comparison_after_filtering_mt0.5-5_feature200-7500.pdf"))
VlnPlot(HL_E105_rep1, features = "nCount_RNA", pt.size=0)+xlab("orig.ident")+ylab("# UMI")+labs(title="# UMI")
VlnPlot(HL_E105_rep1, features = "nFeature_RNA", pt.size=0)+xlab("orig.ident")+ylab("# genes")+labs(title="# genes")
VlnPlot(HL_E105_rep1, features = "percent.mt", pt.size=0)+xlab("orig.ident")+ylab("% mito.reads")+labs(title="% mito.reads")
FeatureScatter(HL_E105_rep1, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size=0.2)+xlab("# UMI")+ylab("# percent.mt")
FeatureScatter(HL_E105_rep1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.2)+xlab("# UMI")+ylab("# genes")
dev.off()

# save as seurat object 
saveRDS(HL_E105_rep1, file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E105_rep1_filtered.rds")

#### A.2 -> HL_E105_rep2 ####
#First step load rawdata coming from CellRanger into a Seurat object 
#This Seurat Object includes only cells with more 200 genes detected and the genes present in at least 3 cells
HL_E105_rep2.data <- Read10X(data.dir = "/Volumes/RAQUEL-LAB3/Single_cell_58_data/CellRanger_58/cr_HL_E105_58_REP2/outs/filtered_feature_bc_matrix/")

## since the previous files include some extra labelling attached to the cellname I'm gonna remove this extra tag
#data
colnames(HL_E105_rep2.data)
colnames(HL_E105_rep2.data) <- gsub("-1$","",colnames(HL_E105_rep2.data))
colnames(HL_E105_rep2.data)

#Velocyto.data
HL_E105_rep2.velo <- ReadVelocity(file = "/Volumes/RAQUEL-LAB3/Single_cell_58_data/CellRanger_58/cr_HL_E105_58_REP2/velocyto/cr_HL_E105_58_REP2.loom")

## since the previous files include some extra labelling attached to the cellname I'm gonna remove this extra tag
#velo
colnames(HL_E105_rep2.velo$spliced)
colnames(HL_E105_rep2.velo$spliced) <- gsub("cr_HL_E105_58_REP2:","",colnames(HL_E105_rep2.velo$spliced))
colnames(HL_E105_rep2.velo$unspliced) <- gsub("cr_HL_E105_58_REP2:","",colnames(HL_E105_rep2.velo$unspliced))
colnames(HL_E105_rep2.velo$ambiguous) <- gsub("cr_HL_E105_58_REP2:","",colnames(HL_E105_rep2.velo$ambiguous))
colnames(HL_E105_rep2.velo$spliced)

colnames(HL_E105_rep2.velo$spliced) <- gsub("x$","",colnames(HL_E105_rep2.velo$spliced))
colnames(HL_E105_rep2.velo$unspliced) <- gsub("x$","",colnames(HL_E105_rep2.velo$unspliced))
colnames(HL_E105_rep2.velo$ambiguous) <- gsub("x$","",colnames(HL_E105_rep2.velo$ambiguous))
colnames(HL_E105_rep2.velo$spliced)

#test
A<-colnames(HL_E105_rep2.data)
B<-colnames(HL_E105_rep2.velo$ambiguous)
## Test to confirm if the cell names (Colnames) are equal in .data and .velo ##
#if I use all.equal(A,B) without sorting I got this message "5251 string mismatches" 
all.equal(A,B)
#but if I sort B, before running all.equal(A,B) then the answer is TRUE
B <- sort(B)
all.equal(A,B)
#also I've check the similarities by checkin the items in A that are not in B
# it can be done with %in%
A[!(A %in% B)]
# it can be also done with setdiff
setdiff(A,B)
#items in B that are not in A
B[!(B %in% A)]
setdiff(B,A)

## Checking if there are differences in the genes between the file .data (a) and .velo (b) ##
b <- rownames(HL_E105_rep2.velo$spliced)
a <- rownames(HL_E105_rep2.data)
a <- sort(a)
b <- sort(b)
length(a)
length(b)

## Merging both Matrix in one list to generate the Seurat object ##
HL_E105_rep2 <- HL_E105_rep2.velo
HL_E105_rep2$RNA <- HL_E105_rep2.data

## checking the data inside the Matrix ##
HL_E105_rep2$spliced[c("Pitx1", "Shox2", "Hoxa13","mt-Co2"),1:30]
HL_E105_rep2$RNA[c("Pitx1", "Shox2", "Hoxa13","mt-Co2"),1:30]

#create Seurat object
HL_E105_rep2 <- as.Seurat(x = HL_E105_rep2)
head(HL_E105_rep2@meta.data,5)
tail(HL_E105_rep2@meta.data,5)
HL_E105_rep2

DefaultAssay(object = HL_E105_rep2) <- "RNA"
HL_E105_rep2

rm(a,b,A,B,HL_E105_rep2.data,HL_E105_rep2.velo)

# Rename orig.ident for future analysis #
levels(HL_E105_rep2$orig.ident) <- "HL_E105_rep2"
head(HL_E105_rep2@meta.data,5)
#Rename levels also
levels(HL_E105_rep2)
HL_E105_rep2 <- RenameIdents(HL_E105_rep2, "SeuratProject" = "HL_E105_rep2")
levels(HL_E105_rep2)

#calculate % of reads belonging to Mitochondrial genes
HL_E105_rep2[["percent.mt"]] <- PercentageFeatureSet(HL_E105_rep2, pattern = "^mt-",assay = "RNA")
head(HL_E105_rep2@meta.data,5)
summary(HL_E105_rep2@meta.data)

#count cells in the seurat object
t1 <- table(HL_E105_rep2@meta.data[["orig.ident"]])
capture.output(bind_rows(t1), file=paste0("numbers_total.txt"))
t1

### A.1.2 QUALITY CONTROL 1 (QC_1) ###
FeatureScatter(HL_E105_rep2, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(HL_E105_rep2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
VlnPlot(HL_E105_rep2, features = "nCount_RNA", pt.size=0)
VlnPlot(HL_E105_rep2, features = "nFeature_RNA", pt.size=0)

pdf(file=paste0("HL_58_E105_QC_rep2_comparison_before_filtering.pdf"))
VlnPlot(HL_E105_rep2, features = "nCount_RNA", pt.size=0)+xlab("orig.ident")+ylab("# UMI")+labs(title="# UMI")
VlnPlot(HL_E105_rep2, features = "nFeature_RNA", pt.size=0)+xlab("orig.ident")+ylab("# genes")+labs(title="# genes")
VlnPlot(HL_E105_rep2, features = "percent.mt", pt.size=0)+xlab("orig.ident")+ylab("% mito.reads")+labs(title="% mito.reads")
FeatureScatter(HL_E105_rep2, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size=0.2)+xlab("# UMI")+ylab("# percent.mt")
FeatureScatter(HL_E105_rep2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.2)+xlab("# UMI")+ylab("# genes")
dev.off()

# save as seurat object 
saveRDS(HL_E105_rep2, file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E105_rep2_unfiltered.rds")

### A.1. filtering FINAL --> MT > 1 < 5 & nfeatures >200 <7500 ###
HL_E105_rep2 <- subset(HL_E105_rep2, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt >1 & percent.mt < 5)
FeatureScatter(HL_E105_rep2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(HL_E105_rep2, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size=0.2)+xlab("# UMI")+ylab("# percent.mt")
table(HL_E105_rep2@meta.data[["orig.ident"]])

pdf(file=paste0("HL_58_E105_QC_rep2_comparison_after_filtering_mt1-5_feature200-7500.pdf"))
VlnPlot(HL_E105_rep2, features = "nCount_RNA", pt.size=0)+xlab("orig.ident")+ylab("# UMI")+labs(title="# UMI")
VlnPlot(HL_E105_rep2, features = "nFeature_RNA", pt.size=0)+xlab("orig.ident")+ylab("# genes")+labs(title="# genes")
VlnPlot(HL_E105_rep2, features = "percent.mt", pt.size=0)+xlab("orig.ident")+ylab("% mito.reads")+labs(title="% mito.reads")
FeatureScatter(HL_E105_rep2, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size=0.2)+xlab("# UMI")+ylab("# percent.mt")
FeatureScatter(HL_E105_rep2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.2)+xlab("# UMI")+ylab("# genes")
dev.off()

# save as seurat object 
saveRDS(HL_E105_rep2, file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E105_rep2_filtered.rds")

#### A.3 -> HL_E115_rep1 ####
#First step load rawdata coming from CellRanger into a Seurat object 
#This Seurat Object includes only cells with more 200 genes detected and the genes present in at least 3 cells
HL_E115_rep1.data <- Read10X(data.dir = "/Volumes/RAQUEL-LAB3/Single_cell_58_data/CellRanger_58/cr_HL_E115_58_REP1/outs/filtered_feature_bc_matrix/")

## since the previous files include some extra labelling attached to the cellname I'm gonna remove this extra tag
#data
colnames(HL_E115_rep1.data)
colnames(HL_E115_rep1.data) <- gsub("-1$","",colnames(HL_E115_rep1.data))
colnames(HL_E115_rep1.data)

#Velocyto.data
HL_E115_rep1.velo <- ReadVelocity(file = "/Volumes/RAQUEL-LAB3/Single_cell_58_data/CellRanger_58/cr_HL_E115_58_REP1/velocyto/cr_HL_E115_58_REP1.loom")

## since the previous files include some extra labelling attached to the cellname I'm gonna remove this extra tag
#velo
colnames(HL_E115_rep1.velo$spliced)
colnames(HL_E115_rep1.velo$spliced) <- gsub("cr_HL_E115_58_REP1:","",colnames(HL_E115_rep1.velo$spliced))
colnames(HL_E115_rep1.velo$unspliced) <- gsub("cr_HL_E115_58_REP1:","",colnames(HL_E115_rep1.velo$unspliced))
colnames(HL_E115_rep1.velo$ambiguous) <- gsub("cr_HL_E115_58_REP1:","",colnames(HL_E115_rep1.velo$ambiguous))
colnames(HL_E115_rep1.velo$spliced)

colnames(HL_E115_rep1.velo$spliced) <- gsub("x$","",colnames(HL_E115_rep1.velo$spliced))
colnames(HL_E115_rep1.velo$unspliced) <- gsub("x$","",colnames(HL_E115_rep1.velo$unspliced))
colnames(HL_E115_rep1.velo$ambiguous) <- gsub("x$","",colnames(HL_E115_rep1.velo$ambiguous))
colnames(HL_E115_rep1.velo$spliced)

#test
A<-colnames(HL_E115_rep1.data)
B<-colnames(HL_E115_rep1.velo$ambiguous)
## Test to confirm if the cell names (Colnames) are equal in .data and .velo ##
#if I use all.equal(A,B) without sorting I got this message "3947 string mismatches" 
all.equal(A,B)
#but if I sort B, before running all.equal(A,B) then the answer is TRUE
B <- sort(B)
all.equal(A,B)
#also I've check the similarities by checkin the items in A that are not in B
# it can be done with %in%
A[!(A %in% B)]
# it can be also done with setdiff
setdiff(A,B)
#items in B that are not in A
B[!(B %in% A)]
setdiff(B,A)

## Checking if there are differences in the genes between the file .data (a) and .velo (b) ##
b <- rownames(HL_E115_rep1.velo$spliced)
a <- rownames(HL_E115_rep1.data)
a <- sort(a)
b <- sort(b)
length(a)
length(b)

## Merging both Matrix in one list to generate the Seurat object ##
HL_E115_rep1 <- HL_E115_rep1.velo
HL_E115_rep1$RNA <- HL_E115_rep1.data

## checking the data inside the Matrix ##
HL_E115_rep1$spliced[c("Pitx1", "Shox2", "Hoxa13","mt-Co2"),1:30]
HL_E115_rep1$RNA[c("Pitx1", "Shox2", "Hoxa13","mt-Co2"),1:30]

#create Seurat object
HL_E115_rep1 <- as.Seurat(x = HL_E115_rep1)
head(HL_E115_rep1@meta.data,5)
tail(HL_E115_rep1@meta.data,5)
HL_E115_rep1

DefaultAssay(object = HL_E115_rep1) <- "RNA"
HL_E115_rep1

rm(a,b,A,B,HL_E115_rep1.data,HL_E115_rep1.velo)

# Rename orig.ident for future analysis #
levels(HL_E115_rep1$orig.ident) <- "HL_E115_rep1"
head(HL_E115_rep1@meta.data,5)
#Rename levels also
levels(HL_E115_rep1)
HL_E115_rep1 <- RenameIdents(HL_E115_rep1, "SeuratProject" = "HL_E115_rep1")
levels(HL_E115_rep1)

#calculate % of reads belonging to Mitochondrial genes
HL_E115_rep1[["percent.mt"]] <- PercentageFeatureSet(HL_E115_rep1, pattern = "^mt-",assay = "RNA")
head(HL_E115_rep1@meta.data,5)
summary(HL_E115_rep1@meta.data)

#count cells in the seurat object
t1 <- table(HL_E115_rep1@meta.data[["orig.ident"]])
capture.output(bind_rows(t1), file=paste0("numbers_total.txt"))
t1

### A.1.2 QUALITY CONTROL 1 (QC_1) ###
FeatureScatter(HL_E115_rep1, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(HL_E115_rep1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
VlnPlot(HL_E115_rep1, features = "nCount_RNA", pt.size=0)
VlnPlot(HL_E115_rep1, features = "nFeature_RNA", pt.size=0)

pdf(file=paste0("HL_58_E115_QC_rep1_comparison_before_filtering.pdf"))
VlnPlot(HL_E115_rep1, features = "nCount_RNA", pt.size=0)+xlab("orig.ident")+ylab("# UMI")+labs(title="# UMI")
VlnPlot(HL_E115_rep1, features = "nFeature_RNA", pt.size=0)+xlab("orig.ident")+ylab("# genes")+labs(title="# genes")
VlnPlot(HL_E115_rep1, features = "percent.mt", pt.size=0)+xlab("orig.ident")+ylab("% mito.reads")+labs(title="% mito.reads")
FeatureScatter(HL_E115_rep1, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size=0.2)+xlab("# UMI")+ylab("# percent.mt")
FeatureScatter(HL_E115_rep1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.2)+xlab("# UMI")+ylab("# genes")
dev.off()

# save as seurat object 
saveRDS(HL_E115_rep1, file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E115_rep1_unfiltered.rds")

### A.1. filtering FINAL --> MT > 1 < 5 & nfeatures >200 <7500 ###
HL_E115_rep1 <- readRDS(file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E115_rep1_unfiltered.rds")
HL_E115_rep1 <- subset(HL_E115_rep1, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt >1 & percent.mt < 5)
FeatureScatter(HL_E115_rep1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(HL_E115_rep1, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size=0.2)+xlab("# UMI")+ylab("# percent.mt")
table(HL_E115_rep1@meta.data[["orig.ident"]])

pdf(file=paste0("HL_58_E115_QC_rep1_comparison_after_filtering_mt1-5_feature200-7500.pdf"))
VlnPlot(HL_E115_rep1, features = "nCount_RNA", pt.size=0)+xlab("orig.ident")+ylab("# UMI")+labs(title="# UMI")
VlnPlot(HL_E115_rep1, features = "nFeature_RNA", pt.size=0)+xlab("orig.ident")+ylab("# genes")+labs(title="# genes")
VlnPlot(HL_E115_rep1, features = "percent.mt", pt.size=0)+xlab("orig.ident")+ylab("% mito.reads")+labs(title="% mito.reads")
FeatureScatter(HL_E115_rep1, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size=0.2)+xlab("# UMI")+ylab("# percent.mt")
FeatureScatter(HL_E115_rep1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.2)+xlab("# UMI")+ylab("# genes")
dev.off()

# save as seurat object 
saveRDS(HL_E115_rep1, file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E115_rep1_filtered.rds")

#### A.4 -> HL_E115_rep2 ####
#First step load rawdata coming from CellRanger into a Seurat object 
#This Seurat Object includes only cells with more 200 genes detected and the genes present in at least 3 cells
HL_E115_rep2.data <- Read10X(data.dir = "/Volumes/RAQUEL-LAB3/Single_cell_58_data/CellRanger_58/cr_HL_E115_58_REP2/outs/filtered_feature_bc_matrix/")

## since the previous files include some extra labelling attached to the cellname I'm gonna remove this extra tag
#data
colnames(HL_E115_rep2.data)
colnames(HL_E115_rep2.data) <- gsub("-1$","",colnames(HL_E115_rep2.data))
colnames(HL_E115_rep2.data)

#Velocyto.data
HL_E115_rep2.velo <- ReadVelocity(file = "/Volumes/RAQUEL-LAB3/Single_cell_58_data/CellRanger_58/cr_HL_E115_58_REP2/velocyto/cr_HL_E115_58_REP2.loom")

## since the previous files include some extra labelling attached to the cellname I'm gonna remove this extra tag
#velo
colnames(HL_E115_rep2.velo$spliced)
colnames(HL_E115_rep2.velo$spliced) <- gsub("cr_HL_E115_58_REP2:","",colnames(HL_E115_rep2.velo$spliced))
colnames(HL_E115_rep2.velo$unspliced) <- gsub("cr_HL_E115_58_REP2:","",colnames(HL_E115_rep2.velo$unspliced))
colnames(HL_E115_rep2.velo$ambiguous) <- gsub("cr_HL_E115_58_REP2:","",colnames(HL_E115_rep2.velo$ambiguous))
colnames(HL_E115_rep2.velo$spliced)

colnames(HL_E115_rep2.velo$spliced) <- gsub("x$","",colnames(HL_E115_rep2.velo$spliced))
colnames(HL_E115_rep2.velo$unspliced) <- gsub("x$","",colnames(HL_E115_rep2.velo$unspliced))
colnames(HL_E115_rep2.velo$ambiguous) <- gsub("x$","",colnames(HL_E115_rep2.velo$ambiguous))
colnames(HL_E115_rep2.velo$spliced)

#test
A<-colnames(HL_E115_rep2.data)
B<-colnames(HL_E115_rep2.velo$ambiguous)
## Test to confirm if the cell names (Colnames) are equal in .data and .velo ##
#if I use all.equal(A,B) without sorting I got this message "5251 string mismatches" 
all.equal(A,B)
#but if I sort B, before running all.equal(A,B) then the answer is TRUE
B <- sort(B)
all.equal(A,B)
#also I've check the similarities by checkin the items in A that are not in B
# it can be done with %in%
A[!(A %in% B)]
# it can be also done with setdiff
setdiff(A,B)
#items in B that are not in A
B[!(B %in% A)]
setdiff(B,A)

## Checking if there are differences in the genes between the file .data (a) and .velo (b) ##
b <- rownames(HL_E115_rep2.velo$spliced)
a <- rownames(HL_E115_rep2.data)
a <- sort(a)
b <- sort(b)
length(a)
length(b)

## Merging both Matrix in one list to generate the Seurat object ##
HL_E115_rep2 <- HL_E115_rep2.velo
HL_E115_rep2$RNA <- HL_E115_rep2.data

## checking the data inside the Matrix ##
HL_E115_rep2$spliced[c("Pitx1", "Shox2", "Hoxa13","mt-Co2"),1:30]
HL_E115_rep2$RNA[c("Pitx1", "Shox2", "Hoxa13","mt-Co2"),1:30]

#create Seurat object
HL_E115_rep2 <- as.Seurat(x = HL_E115_rep2)
head(HL_E115_rep2@meta.data,5)
tail(HL_E115_rep2@meta.data,5)
HL_E115_rep2

DefaultAssay(object = HL_E115_rep2) <- "RNA"
HL_E115_rep2

rm(a,b,A,B,HL_E115_rep2.data,HL_E115_rep2.velo)

# Rename orig.ident for future analysis #
levels(HL_E115_rep2$orig.ident) <- "HL_E115_rep2"
head(HL_E115_rep2@meta.data,5)
#Rename levels also
levels(HL_E115_rep2)
HL_E115_rep2 <- RenameIdents(HL_E115_rep2, "SeuratProject" = "HL_E115_rep2")
levels(HL_E115_rep2)

#calculate % of reads belonging to Mitochondrial genes
HL_E115_rep2[["percent.mt"]] <- PercentageFeatureSet(HL_E115_rep2, pattern = "^mt-",assay = "RNA")
head(HL_E115_rep2@meta.data,5)
summary(HL_E115_rep2@meta.data)

#count cells in the seurat object
t1 <- table(HL_E115_rep2@meta.data[["orig.ident"]])
capture.output(bind_rows(t1), file=paste0("numbers_total.txt"))
t1

### A.1.2 QUALITY CONTROL 1 (QC_1) ###
FeatureScatter(HL_E115_rep2, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(HL_E115_rep2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
VlnPlot(HL_E115_rep2, features = "nCount_RNA", pt.size=0)
VlnPlot(HL_E115_rep2, features = "nFeature_RNA", pt.size=0)

pdf(file=paste0("HL_58_E115_QC_rep2_comparison_before_filtering.pdf"))
VlnPlot(HL_E115_rep2, features = "nCount_RNA", pt.size=0)+xlab("orig.ident")+ylab("# UMI")+labs(title="# UMI")
VlnPlot(HL_E115_rep2, features = "nFeature_RNA", pt.size=0)+xlab("orig.ident")+ylab("# genes")+labs(title="# genes")
VlnPlot(HL_E115_rep2, features = "percent.mt", pt.size=0)+xlab("orig.ident")+ylab("% mito.reads")+labs(title="% mito.reads")
FeatureScatter(HL_E115_rep2, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size=0.2)+xlab("# UMI")+ylab("# percent.mt")
FeatureScatter(HL_E115_rep2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.2)+xlab("# UMI")+ylab("# genes")
dev.off()

# save as seurat object 
saveRDS(HL_E115_rep2, file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E115_rep2_unfiltered.rds")

### A.1. filtering FINAL --> MT > 0.5 < 5 & nfeatures >200 <7500 ###
HL_E115_rep2 <-readRDS(file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E115_rep2_unfiltered.rds")
HL_E115_rep2 <- subset(HL_E115_rep2, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt >0.5 & percent.mt < 5)
FeatureScatter(HL_E115_rep2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(HL_E115_rep2, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size=0.2)+xlab("# UMI")+ylab("# percent.mt")
table(HL_E115_rep2@meta.data[["orig.ident"]])

pdf(file=paste0("HL_58_E115_QC_rep2_comparison_after_filtering_mt0.5-5_feature200-7500.pdf"))
VlnPlot(HL_E115_rep2, features = "nCount_RNA", pt.size=0)+xlab("orig.ident")+ylab("# UMI")+labs(title="# UMI")
VlnPlot(HL_E115_rep2, features = "nFeature_RNA", pt.size=0)+xlab("orig.ident")+ylab("# genes")+labs(title="# genes")
VlnPlot(HL_E115_rep2, features = "percent.mt", pt.size=0)+xlab("orig.ident")+ylab("% mito.reads")+labs(title="% mito.reads")
FeatureScatter(HL_E115_rep2, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size=0.2)+xlab("# UMI")+ylab("# percent.mt")
FeatureScatter(HL_E115_rep2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.2)+xlab("# UMI")+ylab("# genes")
dev.off()

# save as seurat object 
saveRDS(HL_E115_rep2, file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E115_rep2_filtered.rds")

#### A.5 -> HL_E125_rep1 ####
#First step load rawdata coming from CellRanger into a Seurat object 
#This Seurat Object includes only cells with more 200 genes detected and the genes present in at least 3 cells
HL_E125_rep1.data <- Read10X(data.dir = "/Volumes/RAQUEL-LAB3/Single_cell_58_data/CellRanger_58/cr_HL_E125_58_REP1/outs/filtered_feature_bc_matrix/")

## since the previous files include some extra labelling attached to the cellname I'm gonna remove this extra tag
#data
colnames(HL_E125_rep1.data)
colnames(HL_E125_rep1.data) <- gsub("-1$","",colnames(HL_E125_rep1.data))
colnames(HL_E125_rep1.data)

#Velocyto.data
HL_E125_rep1.velo <- ReadVelocity(file = "/Volumes/RAQUEL-LAB3/Single_cell_58_data/CellRanger_58/cr_HL_E125_58_REP1/velocyto/cr_HL_E125_58_REP1.loom")

## since the previous files include some extra labelling attached to the cellname I'm gonna remove this extra tag
#velo
colnames(HL_E125_rep1.velo$spliced)
colnames(HL_E125_rep1.velo$spliced) <- gsub("cr_HL_E125_58_REP1:","",colnames(HL_E125_rep1.velo$spliced))
colnames(HL_E125_rep1.velo$unspliced) <- gsub("cr_HL_E125_58_REP1:","",colnames(HL_E125_rep1.velo$unspliced))
colnames(HL_E125_rep1.velo$ambiguous) <- gsub("cr_HL_E125_58_REP1:","",colnames(HL_E125_rep1.velo$ambiguous))
colnames(HL_E125_rep1.velo$spliced)

colnames(HL_E125_rep1.velo$spliced) <- gsub("x$","",colnames(HL_E125_rep1.velo$spliced))
colnames(HL_E125_rep1.velo$unspliced) <- gsub("x$","",colnames(HL_E125_rep1.velo$unspliced))
colnames(HL_E125_rep1.velo$ambiguous) <- gsub("x$","",colnames(HL_E125_rep1.velo$ambiguous))
colnames(HL_E125_rep1.velo$spliced)

#test
A<-colnames(HL_E125_rep1.data)
B<-colnames(HL_E125_rep1.velo$ambiguous)
## Test to confirm if the cell names (Colnames) are equal in .data and .velo ##
#if I use all.equal(A,B) without sorting I got this message "3947 string mismatches" 
all.equal(A,B)
#but if I sort B, before running all.equal(A,B) then the answer is TRUE
B <- sort(B)
all.equal(A,B)
#also I've check the similarities by checkin the items in A that are not in B
# it can be done with %in%
A[!(A %in% B)]
# it can be also done with setdiff
setdiff(A,B)
#items in B that are not in A
B[!(B %in% A)]
setdiff(B,A)

## Checking if there are differences in the genes between the file .data (a) and .velo (b) ##
b <- rownames(HL_E125_rep1.velo$spliced)
a <- rownames(HL_E125_rep1.data)
a <- sort(a)
b <- sort(b)
length(a)
length(b)

## Merging both Matrix in one list to generate the Seurat object ##
HL_E125_rep1 <- HL_E125_rep1.velo
HL_E125_rep1$RNA <- HL_E125_rep1.data

## checking the data inside the Matrix ##
HL_E125_rep1$spliced[c("Pitx1", "Shox2", "Hoxa13","mt-Co2"),1:30]
HL_E125_rep1$RNA[c("Pitx1", "Shox2", "Hoxa13","mt-Co2"),1:30]

#create Seurat object
HL_E125_rep1 <- as.Seurat(x = HL_E125_rep1)
head(HL_E125_rep1@meta.data,5)
tail(HL_E125_rep1@meta.data,5)
HL_E125_rep1

DefaultAssay(object = HL_E125_rep1) <- "RNA"
HL_E125_rep1

rm(a,b,A,B,HL_E125_rep1.data,HL_E125_rep1.velo)

# Rename orig.ident for future analysis #
levels(HL_E125_rep1$orig.ident) <- "HL_E125_rep1"
head(HL_E125_rep1@meta.data,5)
#Rename levels also
levels(HL_E125_rep1)
HL_E125_rep1 <- RenameIdents(HL_E125_rep1, "SeuratProject" = "HL_E125_rep1")
levels(HL_E125_rep1)

#calculate % of reads belonging to Mitochondrial genes
HL_E125_rep1[["percent.mt"]] <- PercentageFeatureSet(HL_E125_rep1, pattern = "^mt-",assay = "RNA")
head(HL_E125_rep1@meta.data,5)
summary(HL_E125_rep1@meta.data)

#count cells in the seurat object
t1 <- table(HL_E125_rep1@meta.data[["orig.ident"]])
capture.output(bind_rows(t1), file=paste0("numbers_total.txt"))
t1

### A.1.2 QUALITY CONTROL 1 (QC_1) ###
FeatureScatter(HL_E125_rep1, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(HL_E125_rep1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
VlnPlot(HL_E125_rep1, features = "nCount_RNA", pt.size=0)
VlnPlot(HL_E125_rep1, features = "nFeature_RNA", pt.size=0)

pdf(file=paste0("HL_58_E125_QC_rep1_comparison_before_filtering.pdf"))
VlnPlot(HL_E125_rep1, features = "nCount_RNA", pt.size=0)+xlab("orig.ident")+ylab("# UMI")+labs(title="# UMI")
VlnPlot(HL_E125_rep1, features = "nFeature_RNA", pt.size=0)+xlab("orig.ident")+ylab("# genes")+labs(title="# genes")
VlnPlot(HL_E125_rep1, features = "percent.mt", pt.size=0)+xlab("orig.ident")+ylab("% mito.reads")+labs(title="% mito.reads")
FeatureScatter(HL_E125_rep1, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size=0.2)+xlab("# UMI")+ylab("# percent.mt")
FeatureScatter(HL_E125_rep1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.2)+xlab("# UMI")+ylab("# genes")
dev.off()

# save as seurat object 
saveRDS(HL_E125_rep1, file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E125_rep1_unfiltered.rds")

### A.1. filtering FINAL --> MT > 0.5 < 5 & nfeatures >200 <7500 ###
HL_E125_rep1 <-readRDS(file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E125_rep1_unfiltered.rds")
HL_E125_rep1 <- subset(HL_E125_rep1, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt >0.5 & percent.mt < 5)
FeatureScatter(HL_E125_rep1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(HL_E125_rep1, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size=0.2)+xlab("# UMI")+ylab("# percent.mt")
table(HL_E125_rep1@meta.data[["orig.ident"]])

pdf(file=paste0("HL_58_E125_QC_rep1_comparison_after_filtering_mt0.5-5_feature200-7500.pdf"))
VlnPlot(HL_E125_rep1, features = "nCount_RNA", pt.size=0)+xlab("orig.ident")+ylab("# UMI")+labs(title="# UMI")
VlnPlot(HL_E125_rep1, features = "nFeature_RNA", pt.size=0)+xlab("orig.ident")+ylab("# genes")+labs(title="# genes")
VlnPlot(HL_E125_rep1, features = "percent.mt", pt.size=0)+xlab("orig.ident")+ylab("% mito.reads")+labs(title="% mito.reads")
FeatureScatter(HL_E125_rep1, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size=0.2)+xlab("# UMI")+ylab("# percent.mt")
FeatureScatter(HL_E125_rep1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.2)+xlab("# UMI")+ylab("# genes")
dev.off()

# save as seurat object 
saveRDS(HL_E125_rep1, file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E125_rep1_filtered.rds")

##### A.6 -> HL_E135_rep1 ####
#First step load rawdata coming from CellRanger into a Seurat object 
#This Seurat Object includes only cells with more 200 genes detected and the genes present in at least 3 cells
HL_E135_rep1.data <- Read10X(data.dir = "/Volumes/RAQUEL-LAB3/Single_cell_58_data/CellRanger_58/cr_HL_E135_58/outs/filtered_feature_bc_matrix/")

## since the previous files include some extra labelling attached to the cellname I'm gonna remove this extra tag
#data
colnames(HL_E135_rep1.data)
colnames(HL_E135_rep1.data) <- gsub("-1$","",colnames(HL_E135_rep1.data))
colnames(HL_E135_rep1.data)

#Velocyto.data
HL_E135_rep1.velo <- ReadVelocity(file = "/Volumes/RAQUEL-LAB3/Single_cell_58_data/CellRanger_58/cr_HL_E135_58/velocyto/cr_HL_E135_58.loom")

## since the previous files include some extra labelling attached to the cellname I'm gonna remove this extra tag
#velo
colnames(HL_E135_rep1.velo$spliced)
colnames(HL_E135_rep1.velo$spliced) <- gsub("cr_HL_E135_58:","",colnames(HL_E135_rep1.velo$spliced))
colnames(HL_E135_rep1.velo$unspliced) <- gsub("cr_HL_E135_58:","",colnames(HL_E135_rep1.velo$unspliced))
colnames(HL_E135_rep1.velo$ambiguous) <- gsub("cr_HL_E135_58:","",colnames(HL_E135_rep1.velo$ambiguous))
colnames(HL_E135_rep1.velo$spliced)

colnames(HL_E135_rep1.velo$spliced) <- gsub("x$","",colnames(HL_E135_rep1.velo$spliced))
colnames(HL_E135_rep1.velo$unspliced) <- gsub("x$","",colnames(HL_E135_rep1.velo$unspliced))
colnames(HL_E135_rep1.velo$ambiguous) <- gsub("x$","",colnames(HL_E135_rep1.velo$ambiguous))
colnames(HL_E135_rep1.velo$spliced)

#test
A<-colnames(HL_E135_rep1.data)
B<-colnames(HL_E135_rep1.velo$ambiguous)
## Test to confirm if the cell names (Colnames) are equal in .data and .velo ##
#if I use all.equal(A,B) without sorting I got this message "3947 string mismatches" 
all.equal(A,B)
#but if I sort B, before running all.equal(A,B) then the answer is TRUE
B <- sort(B)
all.equal(A,B)
#also I've check the similarities by checkin the items in A that are not in B
# it can be done with %in%
A[!(A %in% B)]
# it can be also done with setdiff
setdiff(A,B)
#items in B that are not in A
B[!(B %in% A)]
setdiff(B,A)

## Checking if there are differences in the genes between the file .data (a) and .velo (b) ##
b <- rownames(HL_E135_rep1.velo$spliced)
a <- rownames(HL_E135_rep1.data)
a <- sort(a)
b <- sort(b)
length(a)
length(b)

## Merging both Matrix in one list to generate the Seurat object ##
HL_E135_rep1 <- HL_E135_rep1.velo
HL_E135_rep1$RNA <- HL_E135_rep1.data

## checking the data inside the Matrix ##
HL_E135_rep1$spliced[c("Pitx1", "Shox2", "Hoxa13","mt-Co2"),1:30]
HL_E135_rep1$RNA[c("Pitx1", "Shox2", "Hoxa13","mt-Co2"),1:30]

#create Seurat object
HL_E135_rep1 <- as.Seurat(x = HL_E135_rep1)
head(HL_E135_rep1@meta.data,5)
tail(HL_E135_rep1@meta.data,5)
HL_E135_rep1

DefaultAssay(object = HL_E135_rep1) <- "RNA"
HL_E135_rep1

rm(a,b,A,B,HL_E135_rep1.data,HL_E135_rep1.velo)

# Rename orig.ident for future analysis #
levels(HL_E135_rep1$orig.ident) <- "HL_E135_rep1"
head(HL_E135_rep1@meta.data,5)
#Rename levels also
levels(HL_E135_rep1)
HL_E135_rep1 <- RenameIdents(HL_E135_rep1, "SeuratProject" = "HL_E135_rep1")
levels(HL_E135_rep1)

#calculate % of reads belonging to Mitochondrial genes
HL_E135_rep1[["percent.mt"]] <- PercentageFeatureSet(HL_E135_rep1, pattern = "^mt-",assay = "RNA")
head(HL_E135_rep1@meta.data,5)
summary(HL_E135_rep1@meta.data)

#count cells in the seurat object
t1 <- table(HL_E135_rep1@meta.data[["orig.ident"]])
capture.output(bind_rows(t1), file=paste0("numbers_total.txt"))
t1

### A.1.2 QUALITY CONTROL 1 (QC_1) ###
FeatureScatter(HL_E135_rep1, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(HL_E135_rep1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
VlnPlot(HL_E135_rep1, features = "nCount_RNA", pt.size=0)
VlnPlot(HL_E135_rep1, features = "nFeature_RNA", pt.size=0)

pdf(file=paste0("HL_58_E135_QC_rep1_comparison_before_filtering.pdf"))
VlnPlot(HL_E135_rep1, features = "nCount_RNA", pt.size=0)+xlab("orig.ident")+ylab("# UMI")+labs(title="# UMI")
VlnPlot(HL_E135_rep1, features = "nFeature_RNA", pt.size=0)+xlab("orig.ident")+ylab("# genes")+labs(title="# genes")
VlnPlot(HL_E135_rep1, features = "percent.mt", pt.size=0)+xlab("orig.ident")+ylab("% mito.reads")+labs(title="% mito.reads")
FeatureScatter(HL_E135_rep1, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size=0.2)+xlab("# UMI")+ylab("# percent.mt")
FeatureScatter(HL_E135_rep1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.2)+xlab("# UMI")+ylab("# genes")
dev.off()

# save as seurat object 
saveRDS(HL_E135_rep1, file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E135_rep1_unfiltered.rds")

### A.1. filtering FINAL --> MT > 0.5 < 5 & nfeatures >200 <7500 ###
HL_E135_rep1 <- readRDS(file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E135_rep1_unfiltered.rds")
HL_E135_rep1 <- subset(HL_E135_rep1, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt >0.5 & percent.mt < 5)
FeatureScatter(HL_E135_rep1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(HL_E135_rep1, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size=0.2)+xlab("# UMI")+ylab("# percent.mt")
table(HL_E135_rep1@meta.data[["orig.ident"]])

pdf(file=paste0("HL_58_E135_QC_rep1_comparison_after_filtering_mt0.5-5_feature200-7500.pdf"))
VlnPlot(HL_E135_rep1, features = "nCount_RNA", pt.size=0)+xlab("orig.ident")+ylab("# UMI")+labs(title="# UMI")
VlnPlot(HL_E135_rep1, features = "nFeature_RNA", pt.size=0)+xlab("orig.ident")+ylab("# genes")+labs(title="# genes")
VlnPlot(HL_E135_rep1, features = "percent.mt", pt.size=0)+xlab("orig.ident")+ylab("% mito.reads")+labs(title="% mito.reads")
FeatureScatter(HL_E135_rep1, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size=0.2)+xlab("# UMI")+ylab("# percent.mt")
FeatureScatter(HL_E135_rep1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.2)+xlab("# UMI")+ylab("# genes")
dev.off()

# save as seurat object 
saveRDS(HL_E135_rep1, file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E135_rep1_filtered.rds")

###### A.7 Selected samples to continue analysis #####
#I decided to continue working only with a replicate from each dataset so the number of cells per stage is more or less the same.
#From the E10.5 and the E11.5 I choose the two replicates that in the previous QC were looking the best.
#After filtering MT > 0.5 < 5 & nfeatures >200 <7500
#E10.5_rep1
HL_E105 <- readRDS(file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E105_rep1_filtered.rds")
#E11.5_rep2
HL_E115 <-readRDS(file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E115_rep2_filtered.rds")
#E12.5_rep1
HL_E125 <- readRDS(file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E125_rep1_filtered.rds")
#E13.5_rep1
HL_E135 <- readRDS(file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E135_rep1_filtered.rds")

#### B.Individual samples analysis and doublets removal  #####
#### B.1 HL_E105_rep1 ####
HL_E105 <- readRDS(file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E105_rep1_filtered.rds")
HL_E105
# log normalization
HL_E105 <-NormalizeData(HL_E105)

#Find 2000 most variable features
HL_E105 <-FindVariableFeatures(HL_E105)
top10 <- head(VariableFeatures(HL_E105), 10)
plot1 <- VariableFeaturePlot(HL_E105)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
plot3 <- LabelPoints(plot = plot2, points = c("EYFP", "Shox2", "CRE", "Pitx1", "Hoxd13", "dsmCherry"), repel = TRUE)
plot3

# excluding the reporter CRE from the variable genes to avoid driving the PCA
grep ("CRE", HL_E105@assays[["RNA"]]@var.features)
HL_E105@assays[["RNA"]]@var.features <- grep ("CRE", HL_E105@assays[["RNA"]]@var.features, invert = TRUE, value=TRUE)

# excluding the reporter EYFP from the variable genes to avoid driving the PCA
grep ("EYFP", HL_E105@assays[["RNA"]]@var.features)
HL_E105@assays[["RNA"]]@var.features <- grep ("EYFP", HL_E105@assays[["RNA"]]@var.features, invert = TRUE, value=TRUE)

# excluding the reporter dsmCherry from the variable genes to avoid driving the PCA
grep ("dsmCherry", HL_E105@assays[["RNA"]]@var.features)
HL_E105@assays[["RNA"]]@var.features <- grep ("dsmCherry", HL_E105@assays[["RNA"]]@var.features, invert = TRUE, value=TRUE)

length(VariableFeatures(HL_E105))

#scaling
HL_E105 <-ScaleData(HL_E105)

#pca
HL_E105 <-RunPCA(HL_E105,npcs = 50)

pdf(file=paste0("HL_E105_PCA.pdf"))
DimPlot(HL_E105, reduction = "pca", label = FALSE)
VizDimLoadings(HL_E105, dims = 1:2, reduction = "pca")
ElbowPlot(HL_E105, ndims = 50, reduction = "pca")
dev.off()

#UMAP
HL_E105 <-RunUMAP(HL_E105,dims = 1:50)
DimPlot(HL_E105, reduction = "umap", label = FALSE)

#Doublet Finder
temp1hl <- paramSweep_v3(HL_E105, PCs = 1:50, sct = FALSE)
temp2hl <- summarizeSweep(temp1hl, GT = FALSE)
temp3hl <- find.pK(temp2hl) # BCmetric: look for pK value with highest metric score in temp2
round(0.025*ncol(HL_E105)) #nEXp

HL_E105 <- doubletFinder_v3(HL_E105, PCs = 1:50, pN = 0.25, pK =0.3, nExp =89 , reuse.pANN = FALSE, sct = FALSE)
head(HL_E105@meta.data)
table(HL_E105@meta.data$DF.classifications_0.25_0.3_89)

pdf(file=paste0("HL_E105_UMAP_doubletfinder.pdf"))
DimPlot(HL_E105, reduction = "umap", label = FALSE)
FeatureScatter(HL_E105, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "DF.classifications_0.25_0.3_89")
DimPlot(object = HL_E105, reduction = "umap", group.by = "DF.classifications_0.25_0.3_89", split.by = "DF.classifications_0.25_0.3_89")
dev.off()

HL_E105.singlets <- subset(x=HL_E105, subset =DF.classifications_0.25_0.3_89 == "Singlet")

FeatureScatter(HL_E105.singlets, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "DF.classifications_0.25_0.3_89")
DimPlot(HL_E105.singlets, reduction = "umap", label = FALSE)

summary(HL_E105@meta.data)
summary(HL_E105.singlets@meta.data)

# save as seurat object 
saveRDS(HL_E105.singlets, file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E105_singlets.rds")

rm(HL_E105,plot1,plot2,plot3,temp1hl,temp2hl,temp3hl, top10)

#### B.2 HL_E115_rep2 ####
HL_E115 <- readRDS(file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E115_rep2_filtered.rds")
HL_E115
# log normalization
HL_E115 <-NormalizeData(HL_E115)

#Find 2000 most variable features
HL_E115 <-FindVariableFeatures(HL_E115)
top10 <- head(VariableFeatures(HL_E115), 10)
plot1 <- VariableFeaturePlot(HL_E115)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
plot3 <- LabelPoints(plot = plot2, points = c("EYFP", "Shox2", "CRE", "Pitx1", "Hoxd13", "dsmCherry"), repel = TRUE)
plot3

# excluding the reporter CRE from the variable genes to avoid driving the PCA
grep ("CRE", HL_E115@assays[["RNA"]]@var.features)
HL_E115@assays[["RNA"]]@var.features <- grep ("CRE", HL_E115@assays[["RNA"]]@var.features, invert = TRUE, value=TRUE)

# excluding the reporter EYFP from the variable genes to avoid driving the PCA
grep ("EYFP", HL_E115@assays[["RNA"]]@var.features)
HL_E115@assays[["RNA"]]@var.features <- grep ("EYFP", HL_E115@assays[["RNA"]]@var.features, invert = TRUE, value=TRUE)

# excluding the reporter dsmCherry from the variable genes to avoid driving the PCA
grep ("dsmCherry", HL_E115@assays[["RNA"]]@var.features)
HL_E115@assays[["RNA"]]@var.features <- grep ("dsmCherry", HL_E115@assays[["RNA"]]@var.features, invert = TRUE, value=TRUE)

length(VariableFeatures(HL_E115))

#scaling
HL_E115 <-ScaleData(HL_E115)

#pca
HL_E115 <-RunPCA(HL_E115,npcs = 50)

pdf(file=paste0("HL_E115_PCA.pdf"))
DimPlot(HL_E115, reduction = "pca", label = FALSE)
VizDimLoadings(HL_E115, dims = 1:2, reduction = "pca")
ElbowPlot(HL_E115, ndims = 50, reduction = "pca")
dev.off()

#UMAP
HL_E115 <-RunUMAP(HL_E115,dims = 1:50)
DimPlot(HL_E115, reduction = "umap", label = FALSE)

#Doublet Finder
temp1hl <- paramSweep_v3(HL_E115, PCs = 1:50, sct = FALSE)
temp2hl <- summarizeSweep(temp1hl, GT = FALSE)
temp3hl <- find.pK(temp2hl) # BCmetric: look for pK value with highest metric score in temp2
round(0.025*ncol(HL_E115)) #nEXp

HL_E115 <- doubletFinder_v3(HL_E115, PCs = 1:50, pN = 0.25, pK =0.16, nExp =98 , reuse.pANN = FALSE, sct = FALSE)
head(HL_E115@meta.data)
table(HL_E115@meta.data$DF.classifications_0.25_0.16_98)

pdf(file=paste0("HL_E115_UMAP_doubletfinder.pdf"))
DimPlot(HL_E115, reduction = "umap", label = FALSE)
FeatureScatter(HL_E115, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "DF.classifications_0.25_0.16_98")
DimPlot(object = HL_E115, reduction = "umap", group.by = "DF.classifications_0.25_0.16_98", split.by = "DF.classifications_0.25_0.16_98")
dev.off()

HL_E115.singlets <- subset(x=HL_E115, subset =DF.classifications_0.25_0.16_98 == "Singlet")

FeatureScatter(HL_E115.singlets, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "DF.classifications_0.25_0.16_98")
DimPlot(HL_E115.singlets, reduction = "umap", label = FALSE)

summary(HL_E115@meta.data)
summary(HL_E115.singlets@meta.data)

# save as seurat object 
saveRDS(HL_E115.singlets, file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E115_singlets.rds")

rm(HL_E115,plot1,plot2,plot3,temp1hl,temp2hl,temp3hl, top10)

#### B.3 HL_E125_rep1 ####
HL_E125 <- readRDS(file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E125_rep1_filtered.rds")
HL_E125
# log normalization
HL_E125 <-NormalizeData(HL_E125)

#Find 2000 most variable features
HL_E125 <-FindVariableFeatures(HL_E125)
top10 <- head(VariableFeatures(HL_E125), 10)
plot1 <- VariableFeaturePlot(HL_E125)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
plot3 <- LabelPoints(plot = plot2, points = c("EYFP", "Shox2", "CRE", "Pitx1", "Hoxd13", "dsmCherry"), repel = TRUE)
plot3

# excluding the reporter CRE from the variable genes to avoid driving the PCA
grep ("CRE", HL_E125@assays[["RNA"]]@var.features)
HL_E125@assays[["RNA"]]@var.features <- grep ("CRE", HL_E125@assays[["RNA"]]@var.features, invert = TRUE, value=TRUE)

# excluding the reporter EYFP from the variable genes to avoid driving the PCA
grep ("EYFP", HL_E125@assays[["RNA"]]@var.features)
HL_E125@assays[["RNA"]]@var.features <- grep ("EYFP", HL_E125@assays[["RNA"]]@var.features, invert = TRUE, value=TRUE)

# excluding the reporter dsmCherry from the variable genes to avoid driving the PCA
grep ("dsmCherry", HL_E125@assays[["RNA"]]@var.features)
HL_E125@assays[["RNA"]]@var.features <- grep ("dsmCherry", HL_E125@assays[["RNA"]]@var.features, invert = TRUE, value=TRUE)

length(VariableFeatures(HL_E125))

#scaling
HL_E125 <-ScaleData(HL_E125)

#pca
HL_E125 <-RunPCA(HL_E125,npcs = 50)

pdf(file=paste0("HL_E125_PCA.pdf"))
DimPlot(HL_E125, reduction = "pca", label = FALSE)
VizDimLoadings(HL_E125, dims = 1:2, reduction = "pca")
ElbowPlot(HL_E125, ndims = 50, reduction = "pca")
dev.off()

#UMAP
HL_E125 <-RunUMAP(HL_E125,dims = 1:50)
DimPlot(HL_E125, reduction = "umap", label = FALSE)

#Doublet Finder
temp1hl <- paramSweep_v3(HL_E125, PCs = 1:50, sct = FALSE)
temp2hl <- summarizeSweep(temp1hl, GT = FALSE)
temp3hl <- find.pK(temp2hl) # BCmetric: look for pK value with highest metric score in temp2
round(0.025*ncol(HL_E125)) #nEXp

HL_E125 <- doubletFinder_v3(HL_E125, PCs = 1:50, pN = 0.25, pK =0.25, nExp =71 , reuse.pANN = FALSE, sct = FALSE)
head(HL_E125@meta.data)
table(HL_E125@meta.data$DF.classifications_0.25_0.25_71)

pdf(file=paste0("HL_E125_UMAP_doubletfinder.pdf"))
DimPlot(HL_E125, reduction = "umap", label = FALSE)
FeatureScatter(HL_E125, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "DF.classifications_0.25_0.25_71")
DimPlot(object = HL_E125, reduction = "umap", group.by = "DF.classifications_0.25_0.25_71", split.by = "DF.classifications_0.25_0.25_71")
dev.off()

HL_E125.singlets <- subset(x=HL_E125, subset =DF.classifications_0.25_0.25_71 == "Singlet")

FeatureScatter(HL_E125.singlets, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "DF.classifications_0.25_0.25_71")
DimPlot(HL_E125.singlets, reduction = "umap", label = FALSE)

summary(HL_E125@meta.data)
summary(HL_E125.singlets@meta.data)

# save as seurat object 
saveRDS(HL_E125.singlets, file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E125_singlets.rds")

rm(HL_E125,plot1,plot2,plot3,temp1hl,temp2hl,temp3hl, top10)

#### B.4 HL_E135_rep1 ####
HL_E135 <- readRDS(file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E135_rep1_filtered.rds")
HL_E135
# log normalization
HL_E135 <-NormalizeData(HL_E135)

#Find 2000 most variable features
HL_E135 <-FindVariableFeatures(HL_E135)
top10 <- head(VariableFeatures(HL_E135), 10)
plot1 <- VariableFeaturePlot(HL_E135)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2
plot3 <- LabelPoints(plot = plot2, points = c("EYFP", "Shox2", "CRE", "Pitx1", "Hoxd13", "dsmCherry"), repel = TRUE)
plot3

# excluding the reporter CRE from the variable genes to avoid driving the PCA
grep ("CRE", HL_E135@assays[["RNA"]]@var.features)
HL_E135@assays[["RNA"]]@var.features <- grep ("CRE", HL_E135@assays[["RNA"]]@var.features, invert = TRUE, value=TRUE)

# excluding the reporter EYFP from the variable genes to avoid driving the PCA
grep ("EYFP", HL_E135@assays[["RNA"]]@var.features)
HL_E135@assays[["RNA"]]@var.features <- grep ("EYFP", HL_E135@assays[["RNA"]]@var.features, invert = TRUE, value=TRUE)

# excluding the reporter dsmCherry from the variable genes to avoid driving the PCA
grep ("dsmCherry", HL_E135@assays[["RNA"]]@var.features)
HL_E135@assays[["RNA"]]@var.features <- grep ("dsmCherry", HL_E135@assays[["RNA"]]@var.features, invert = TRUE, value=TRUE)

length(VariableFeatures(HL_E135))

#scaling
HL_E135 <-ScaleData(HL_E135)

#pca
HL_E135 <-RunPCA(HL_E135,npcs = 50)

pdf(file=paste0("HL_E135_PCA.pdf"))
DimPlot(HL_E135, reduction = "pca", label = FALSE)
VizDimLoadings(HL_E135, dims = 1:2, reduction = "pca")
ElbowPlot(HL_E135, ndims = 50, reduction = "pca")
dev.off()

#UMAP
HL_E135 <-RunUMAP(HL_E135,dims = 1:50)
DimPlot(HL_E135, reduction = "umap", label = FALSE)

#Doublet Finder
temp1hl <- paramSweep_v3(HL_E135, PCs = 1:50, sct = FALSE)
temp2hl <- summarizeSweep(temp1hl, GT = FALSE)
temp3hl <- find.pK(temp2hl) # BCmetric: look for pK value with highest metric score in temp2
round(0.025*ncol(HL_E135)) #nEXp

HL_E135 <- doubletFinder_v3(HL_E135, PCs = 1:50, pN = 0.25, pK =0.1, nExp =88 , reuse.pANN = FALSE, sct = FALSE)
head(HL_E135@meta.data)
table(HL_E135@meta.data$DF.classifications_0.25_0.1_88)

pdf(file=paste0("HL_E135_UMAP_doubletfinder.pdf"))
DimPlot(HL_E135, reduction = "umap", label = FALSE)
FeatureScatter(HL_E135, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "DF.classifications_0.25_0.1_88")
DimPlot(object = HL_E135, reduction = "umap", group.by = "DF.classifications_0.25_0.1_88", split.by = "DF.classifications_0.25_0.1_88")
dev.off()

HL_E135.singlets <- subset(x=HL_E135, subset =DF.classifications_0.25_0.1_88 == "Singlet")

FeatureScatter(HL_E135.singlets, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "DF.classifications_0.25_0.1_88")
DimPlot(HL_E135.singlets, reduction = "umap", label = FALSE)

summary(HL_E135@meta.data)
summary(HL_E135.singlets@meta.data)

# save as seurat object 
saveRDS(HL_E135.singlets, file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E135_singlets.rds")

rm(HL_E135,plot1,plot2,plot3,temp1hl,temp2hl,temp3hl, top10)

#### C.QC check for CRE, dsmCherry, EYFP & Shox2 genes & variablefeatures annotation ####
#### C.1 -> HL_E105.singlets ####
HL_E105_rep1.singlets <- readRDS("/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E105_singlets.rds")
HL_E105_rep1.singlets

# flag and count CRE, dsmCherry, EYFP & Shox2 expressing cells
CRE_counts <- FetchData(HL_E105_rep1.singlets, vars="CRE", slot="counts")
dsmCherry_counts <- FetchData(HL_E105_rep1.singlets, vars="dsmCherry", slot="counts")
EYFP_counts <- FetchData(HL_E105_rep1.singlets, vars="EYFP", slot="counts")
Shox2_counts <- FetchData(HL_E105_rep1.singlets, vars="Shox2", slot="counts")

pdf(file=paste0("HL_E105_rep1.singlets_CRE_dsmCherry_Shox2_EYFP_counts_QC.pdf"))
ggplot(data=CRE_counts, aes(CRE_counts$`CRE`)) + geom_histogram(binwidth = 0.5) + labs(x="CRE")
ggplot(data=dsmCherry_counts, aes(dsmCherry_counts$`dsmCherry`)) + geom_histogram(binwidth = 0.5) + labs(x="dsmCherry")
ggplot(data=Shox2_counts, aes(Shox2_counts$`Shox2`)) + geom_histogram(binwidth = 0.5) + labs(x="Shox2")
ggplot(data=EYFP_counts, aes(EYFP_counts$`EYFP`)) + geom_histogram(binwidth = 0.5) + labs(x="EYFP")
ggplot(data=CRE_counts, aes(CRE_counts$`CRE`)) + xlim(0,20) + geom_histogram(binwidth = 0.5) + labs(x="CRE UMIs") 
ggplot(data=dsmCherry_counts, aes(dsmCherry_counts$`dsmCherry`)) + xlim(0,20) + geom_histogram(binwidth = 0.5) + labs(x="dsmCherry UMIs") 
ggplot(data=Shox2_counts, aes(Shox2_counts$`Shox2`)) + xlim(0,20) + geom_histogram(binwidth = 0.5) + labs(x="Shox2 UMIs") 
ggplot(data=EYFP_counts, aes(EYFP_counts$`EYFP`)) + xlim(0,20) + geom_histogram(binwidth = 0.5) + labs(x="EYFP UMIs") 
dev.off()

#CRE_counts
CRE_selected_cells <- WhichCells(HL_E105_rep1.singlets, expression = `CRE`>=1, slot="counts")
length(CRE_selected_cells)
HL_E105_rep1.singlets[["CRE_positive"]] <- colnames(HL_E105_rep1.singlets) %in% CRE_selected_cells

#dsmCherry_counts
dsmCherry_selected_cells <- WhichCells(HL_E105_rep1.singlets, expression = `dsmCherry`>=1, slot="counts")
length(dsmCherry_selected_cells)
HL_E105_rep1.singlets[["dsmCherry_positive"]] <- colnames(HL_E105_rep1.singlets) %in% dsmCherry_selected_cells

#EYFP_counts
EYFP_selected_cells <- WhichCells(HL_E105_rep1.singlets, expression = `EYFP`>=1, slot="counts")
length(EYFP_selected_cells)
HL_E105_rep1.singlets[["EYFP_positive"]] <- colnames(HL_E105_rep1.singlets) %in% EYFP_selected_cells

#Shox2_counts
Shox2_selected_cells <- WhichCells(HL_E105_rep1.singlets, expression = `Shox2`>=1, slot="counts")
length(Shox2_selected_cells)
HL_E105_rep1.singlets[["Shox2_positive"]] <- colnames(HL_E105_rep1.singlets) %in% Shox2_selected_cells

#Positive and negative distribution
Shox2_pos_cells <- WhichCells(HL_E105_rep1.singlets, expression = `Shox2`>=1, slot="counts")
Shox2_pos <- length(Shox2_pos_cells)
Shox2_neg_cells <- WhichCells(HL_E105_rep1.singlets, expression = `Shox2`<1, slot="counts")
Shox2_neg <- length(Shox2_neg_cells)
EYFP_pos_cells <- WhichCells(HL_E105_rep1.singlets, expression = `EYFP`>=1, slot="counts")
EYFP_pos <- length(EYFP_pos_cells)
EYFP_neg_cells <- WhichCells(HL_E105_rep1.singlets, expression = `EYFP`<1, slot="counts")
EYFP_neg <- length(EYFP_neg_cells)
CRE_pos_cells <- WhichCells(HL_E105_rep1.singlets, expression = `CRE`>=1, slot="counts")
CRE_pos <- length(CRE_pos_cells)
CRE_neg_cells <- WhichCells(HL_E105_rep1.singlets, expression = `CRE`<1, slot="counts")
CRE_neg <- length(CRE_neg_cells)

Name <- c("Shox2","EYFP","CRE")
Positive <- c(Shox2_pos,EYFP_pos,CRE_pos)
Negative <- c(Shox2_neg,EYFP_neg,CRE_neg)
df1_HL_E105_rep1.singlets <- data.frame(Name,Positive,Negative)
write.csv(df1_HL_E105_rep1.singlets, file = "df1_HL_E105_rep1.singlets.csv")

EYFPneg_Shox2pos_cells <- WhichCells(HL_E105_rep1.singlets, cells = Shox2_pos_cells , expression = `EYFP` <1, slot="counts")
EYFPneg_Shox2pos <- length(EYFPneg_Shox2pos_cells)

EYFPneg_Shox2neg_cells <- WhichCells(HL_E105_rep1.singlets, cells = Shox2_neg_cells , expression = `EYFP` <1, slot="counts")
EYFPneg_Shox2neg <- length(EYFPneg_Shox2neg_cells)

EYFPpos_Shox2pos_cells <- WhichCells(HL_E105_rep1.singlets, cells = Shox2_pos_cells , expression = `EYFP` >=1, slot="counts")
EYFPpos_Shox2pos <- length(EYFPpos_Shox2pos_cells)

EYFPpos_Shox2neg_cells <- WhichCells(HL_E105_rep1.singlets, cells = Shox2_neg_cells , expression = `EYFP` >=1, slot="counts")
EYFPpos_Shox2neg <- length(EYFPpos_Shox2neg_cells)

EYFPpos_Shox2pos_CREpos_cells <- WhichCells(HL_E105_rep1.singlets, cells = EYFPpos_Shox2pos_cells, expression = `CRE` >=1, slot="counts")
EYFPpos_Shox2pos_CREpos <- length(EYFPpos_Shox2pos_CREpos_cells)

EYFPneg_Shox2pos_CREpos_cells <- WhichCells(HL_E105_rep1.singlets, cells = EYFPneg_Shox2pos_cells, expression = `CRE` >=1, slot="counts")
EYFPneg_Shox2pos_CREpos <- length(EYFPneg_Shox2pos_CREpos_cells)

EYFPpos_Shox2neg_CREpos_cells <- WhichCells(HL_E105_rep1.singlets, cells = EYFPpos_Shox2neg_cells, expression = `CRE` >=1, slot="counts")
EYFPpos_Shox2neg_CREpos <- length(EYFPpos_Shox2neg_CREpos_cells)

EYFPneg_Shox2neg_CREpos_cells <- WhichCells(HL_E105_rep1.singlets, cells = EYFPneg_Shox2neg_cells, expression = `CRE` >=1, slot="counts")
EYFPneg_Shox2neg_CREpos <- length(EYFPneg_Shox2neg_CREpos_cells)

EYFPpos_Shox2pos_CREneg_cells <- WhichCells(HL_E105_rep1.singlets, cells = EYFPpos_Shox2pos_cells, expression = `CRE` <1, slot="counts")
EYFPpos_Shox2pos_CREneg <- length(EYFPpos_Shox2pos_CREneg_cells)

EYFPneg_Shox2pos_CREneg_cells <- WhichCells(HL_E105_rep1.singlets, cells = EYFPneg_Shox2pos_cells, expression = `CRE` <1, slot="counts")
EYFPneg_Shox2pos_CREneg <- length(EYFPneg_Shox2pos_CREneg_cells)

EYFPpos_Shox2neg_CREneg_cells <- WhichCells(HL_E105_rep1.singlets, cells = EYFPpos_Shox2neg_cells, expression = `CRE` <1, slot="counts")
EYFPpos_Shox2neg_CREneg <- length(EYFPpos_Shox2neg_CREneg_cells)

EYFPneg_Shox2neg_CREneg_cells <- WhichCells(HL_E105_rep1.singlets, cells = EYFPneg_Shox2neg_cells, expression = `CRE` <1, slot="counts")
EYFPneg_Shox2neg_CREneg <- length(EYFPneg_Shox2neg_CREneg_cells)

sample1 <- c("EYFPneg_Shox2pos", "EYFPneg_Shox2neg","EYFPpos_Shox2pos","EYFPpos_Shox2neg")
n_cells1 <- c(EYFPneg_Shox2pos, EYFPneg_Shox2neg, EYFPpos_Shox2pos,EYFPpos_Shox2neg)

sample2 <- c("EYFPneg_Shox2pos_CREpos","EYFPneg_Shox2neg_CREpos","EYFPpos_Shox2pos_CREpos", "EYFPpos_Shox2neg_CREpos")
n_cells2 <- c(EYFPneg_Shox2pos_CREpos,EYFPneg_Shox2neg_CREpos,EYFPpos_Shox2pos_CREpos,EYFPpos_Shox2neg_CREpos)

sample3 <- c("EYFPneg_Shox2pos_CREneg","EYFPneg_Shox2neg_CREneg","EYFPpos_Shox2pos_CREneg", "EYFPpos_Shox2neg_CREneg")
n_cells3 <- c(EYFPneg_Shox2pos_CREneg,EYFPneg_Shox2neg_CREneg,EYFPpos_Shox2pos_CREneg,EYFPpos_Shox2neg_CREneg)

df2_HL_E105_rep1.singlets <- data.frame(sample1,n_cells1,sample2,n_cells2,sample3,n_cells3)
write.csv(df2_HL_E105_rep1.singlets, file = "df2_HL_E105_rep1.singlets.csv")

d1 <- data.frame(EYFPpos_Shox2pos_CREpos_cells, "EYFPpos_Shox2pos_CREpos_cells", "A")
colnames(d1) <- c("cells","gene_class","letter")
d2 <- data.frame(EYFPneg_Shox2pos_CREpos_cells, "EYFPneg_Shox2pos_CREpos_cells", "B")
colnames(d2) <- c("cells","gene_class","letter")
d3 <- data.frame(EYFPpos_Shox2neg_CREpos_cells, "EYFPpos_Shox2neg_CREpos_cells", "C")
colnames(d3) <- c("cells","gene_class","letter")
d4 <- data.frame(EYFPneg_Shox2neg_CREpos_cells, "EYFPneg_Shox2neg_CREpos_cells", "D")
colnames(d4) <- c("cells","gene_class","letter")
d5 <- data.frame(EYFPpos_Shox2pos_CREneg_cells, "EYFPpos_Shox2pos_CREneg_cells", "E")
colnames(d5) <- c("cells","gene_class","letter")
d6 <- data.frame(EYFPneg_Shox2pos_CREneg_cells, "EYFPneg_Shox2pos_CREneg_cells", "F")
colnames(d6) <- c("cells","gene_class","letter")
d7 <- data.frame(EYFPpos_Shox2neg_CREneg_cells, "EYFPpos_Shox2neg_CREneg_cells", "G")
colnames(d7) <- c("cells","gene_class","letter")
d8 <- data.frame(EYFPneg_Shox2neg_CREneg_cells, "EYFPneg_Shox2neg_CREneg_cells", "H")
colnames(d8) <- c("cells","gene_class","letter")

cells_class <- rbind(d1,d2,d3,d4,d5,d6,d7,d8)
row.names(x = cells_class) <- cells_class$cells
cells_class2 <- cells_class[,c("gene_class","letter")]
head(cells_class2)

HL_E105_rep1.singlets[[colnames(x = cells_class2)]] <- cells_class2
head(HL_E105_rep1.singlets@meta.data)

pdf(file=paste0("HL_E105_rep1.singlets_CRE_dsmCherry_Shox2_EYFP_proportions.pdf"))

DimPlot(HL_E105_rep1.singlets, reduction = "umap",label = FALSE, group.by = "gene_class")
VlnPlot(object = HL_E105_rep1.singlets, assay = "RNA", features = "Shox2", group.by = "gene_class",
        pt.size = 0.1) 
VlnPlot(object = HL_E105_rep1.singlets, assay = "RNA", features = "CRE", group.by = "gene_class",
        pt.size = 0.1) 
VlnPlot(object = HL_E105_rep1.singlets, assay = "RNA", features = "EYFP", group.by = "gene_class",
        pt.size = 0.1) 

t1 <- table(HL_E105_rep1.singlets@meta.data$gene_class)
total <- sum(t1)
total
t1 <- as.data.frame(t1)
t1
t1$ratio <- t1[,2]/total
sum(t1$ratio)

ggplot(t1, aes(fill= Var1, x=Var1,y=Freq)) + 
  geom_bar(position="dodge", stat = "identity") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1))+
  geom_text(aes(label=Freq), vjust=-0.3, color="black", size=3.5)

ggplot(t1, aes(fill= Var1, x=Var1,y=ratio)) + 
  geom_bar(position="dodge", stat = "identity") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1))+
  geom_text(aes(label = sprintf("%0.2f", round(ratio, 2))), vjust=-0.3, color="black", size=3.5)
dev.off()

# save as seurat object 
saveRDS(HL_E105_rep1.singlets, file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E105.singlets.rds")

#### C.2 -> HL_E115.singlets ####
HL_E115_rep1.singlets <- readRDS("/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E115_singlets.rds")
HL_E115_rep1.singlets

# flag and count CRE, dsmCherry, EYFP & Shox2 expressing cells
CRE_counts <- FetchData(HL_E115_rep1.singlets, vars="CRE", slot="counts")
dsmCherry_counts <- FetchData(HL_E115_rep1.singlets, vars="dsmCherry", slot="counts")
EYFP_counts <- FetchData(HL_E115_rep1.singlets, vars="EYFP", slot="counts")
Shox2_counts <- FetchData(HL_E115_rep1.singlets, vars="Shox2", slot="counts")

pdf(file=paste0("HL_E115_rep1.singlets_CRE_dsmCherry_Shox2_EYFP_counts_QC.pdf"))
ggplot(data=CRE_counts, aes(CRE_counts$`CRE`)) + geom_histogram(binwidth = 0.5) + labs(x="CRE")
ggplot(data=dsmCherry_counts, aes(dsmCherry_counts$`dsmCherry`)) + geom_histogram(binwidth = 0.5) + labs(x="dsmCherry")
ggplot(data=Shox2_counts, aes(Shox2_counts$`Shox2`)) + geom_histogram(binwidth = 0.5) + labs(x="Shox2")
ggplot(data=EYFP_counts, aes(EYFP_counts$`EYFP`)) + geom_histogram(binwidth = 0.5) + labs(x="EYFP")
ggplot(data=CRE_counts, aes(CRE_counts$`CRE`)) + xlim(0,20) + geom_histogram(binwidth = 0.5) + labs(x="CRE UMIs") 
ggplot(data=dsmCherry_counts, aes(dsmCherry_counts$`dsmCherry`)) + xlim(0,20) + geom_histogram(binwidth = 0.5) + labs(x="dsmCherry UMIs") 
ggplot(data=Shox2_counts, aes(Shox2_counts$`Shox2`)) + xlim(0,20) + geom_histogram(binwidth = 0.5) + labs(x="Shox2 UMIs") 
ggplot(data=EYFP_counts, aes(EYFP_counts$`EYFP`)) + xlim(0,20) + geom_histogram(binwidth = 0.5) + labs(x="EYFP UMIs") 
dev.off()

#CRE_counts
CRE_selected_cells <- WhichCells(HL_E115_rep1.singlets, expression = `CRE`>=1, slot="counts")
length(CRE_selected_cells)
HL_E115_rep1.singlets[["CRE_positive"]] <- colnames(HL_E115_rep1.singlets) %in% CRE_selected_cells

#dsmCherry_counts
dsmCherry_selected_cells <- WhichCells(HL_E115_rep1.singlets, expression = `dsmCherry`>=1, slot="counts")
length(dsmCherry_selected_cells)
HL_E115_rep1.singlets[["dsmCherry_positive"]] <- colnames(HL_E115_rep1.singlets) %in% dsmCherry_selected_cells

#EYFP_counts
EYFP_selected_cells <- WhichCells(HL_E115_rep1.singlets, expression = `EYFP`>=1, slot="counts")
length(EYFP_selected_cells)
HL_E115_rep1.singlets[["EYFP_positive"]] <- colnames(HL_E115_rep1.singlets) %in% EYFP_selected_cells

#Shox2_counts
Shox2_selected_cells <- WhichCells(HL_E115_rep1.singlets, expression = `Shox2`>=1, slot="counts")
length(Shox2_selected_cells)
HL_E115_rep1.singlets[["Shox2_positive"]] <- colnames(HL_E115_rep1.singlets) %in% Shox2_selected_cells

#Positive and negative distribution
Shox2_pos_cells <- WhichCells(HL_E115_rep1.singlets, expression = `Shox2`>=1, slot="counts")
Shox2_pos <- length(Shox2_pos_cells)
Shox2_neg_cells <- WhichCells(HL_E115_rep1.singlets, expression = `Shox2`<1, slot="counts")
Shox2_neg <- length(Shox2_neg_cells)
EYFP_pos_cells <- WhichCells(HL_E115_rep1.singlets, expression = `EYFP`>=1, slot="counts")
EYFP_pos <- length(EYFP_pos_cells)
EYFP_neg_cells <- WhichCells(HL_E115_rep1.singlets, expression = `EYFP`<1, slot="counts")
EYFP_neg <- length(EYFP_neg_cells)
CRE_pos_cells <- WhichCells(HL_E115_rep1.singlets, expression = `CRE`>=1, slot="counts")
CRE_pos <- length(CRE_pos_cells)
CRE_neg_cells <- WhichCells(HL_E115_rep1.singlets, expression = `CRE`<1, slot="counts")
CRE_neg <- length(CRE_neg_cells)

Name <- c("Shox2","EYFP","CRE")
Positive <- c(Shox2_pos,EYFP_pos,CRE_pos)
Negative <- c(Shox2_neg,EYFP_neg,CRE_neg)
df1_HL_E115_rep1.singlets <- data.frame(Name,Positive,Negative)
write.csv(df1_HL_E115_rep1.singlets, file = "df1_HL_E115_rep2.singlets.csv")

EYFPneg_Shox2pos_cells <- WhichCells(HL_E115_rep1.singlets, cells = Shox2_pos_cells , expression = `EYFP` <1, slot="counts")
EYFPneg_Shox2pos <- length(EYFPneg_Shox2pos_cells)

EYFPneg_Shox2neg_cells <- WhichCells(HL_E115_rep1.singlets, cells = Shox2_neg_cells , expression = `EYFP` <1, slot="counts")
EYFPneg_Shox2neg <- length(EYFPneg_Shox2neg_cells)

EYFPpos_Shox2pos_cells <- WhichCells(HL_E115_rep1.singlets, cells = Shox2_pos_cells , expression = `EYFP` >=1, slot="counts")
EYFPpos_Shox2pos <- length(EYFPpos_Shox2pos_cells)

EYFPpos_Shox2neg_cells <- WhichCells(HL_E115_rep1.singlets, cells = Shox2_neg_cells , expression = `EYFP` >=1, slot="counts")
EYFPpos_Shox2neg <- length(EYFPpos_Shox2neg_cells)

EYFPpos_Shox2pos_CREpos_cells <- WhichCells(HL_E115_rep1.singlets, cells = EYFPpos_Shox2pos_cells, expression = `CRE` >=1, slot="counts")
EYFPpos_Shox2pos_CREpos <- length(EYFPpos_Shox2pos_CREpos_cells)

EYFPneg_Shox2pos_CREpos_cells <- WhichCells(HL_E115_rep1.singlets, cells = EYFPneg_Shox2pos_cells, expression = `CRE` >=1, slot="counts")
EYFPneg_Shox2pos_CREpos <- length(EYFPneg_Shox2pos_CREpos_cells)

EYFPpos_Shox2neg_CREpos_cells <- WhichCells(HL_E115_rep1.singlets, cells = EYFPpos_Shox2neg_cells, expression = `CRE` >=1, slot="counts")
EYFPpos_Shox2neg_CREpos <- length(EYFPpos_Shox2neg_CREpos_cells)

EYFPneg_Shox2neg_CREpos_cells <- WhichCells(HL_E115_rep1.singlets, cells = EYFPneg_Shox2neg_cells, expression = `CRE` >=1, slot="counts")
EYFPneg_Shox2neg_CREpos <- length(EYFPneg_Shox2neg_CREpos_cells)

EYFPpos_Shox2pos_CREneg_cells <- WhichCells(HL_E115_rep1.singlets, cells = EYFPpos_Shox2pos_cells, expression = `CRE` <1, slot="counts")
EYFPpos_Shox2pos_CREneg <- length(EYFPpos_Shox2pos_CREneg_cells)

EYFPneg_Shox2pos_CREneg_cells <- WhichCells(HL_E115_rep1.singlets, cells = EYFPneg_Shox2pos_cells, expression = `CRE` <1, slot="counts")
EYFPneg_Shox2pos_CREneg <- length(EYFPneg_Shox2pos_CREneg_cells)

EYFPpos_Shox2neg_CREneg_cells <- WhichCells(HL_E115_rep1.singlets, cells = EYFPpos_Shox2neg_cells, expression = `CRE` <1, slot="counts")
EYFPpos_Shox2neg_CREneg <- length(EYFPpos_Shox2neg_CREneg_cells)

EYFPneg_Shox2neg_CREneg_cells <- WhichCells(HL_E115_rep1.singlets, cells = EYFPneg_Shox2neg_cells, expression = `CRE` <1, slot="counts")
EYFPneg_Shox2neg_CREneg <- length(EYFPneg_Shox2neg_CREneg_cells)

sample1 <- c("EYFPneg_Shox2pos", "EYFPneg_Shox2neg","EYFPpos_Shox2pos","EYFPpos_Shox2neg")
n_cells1 <- c(EYFPneg_Shox2pos, EYFPneg_Shox2neg, EYFPpos_Shox2pos,EYFPpos_Shox2neg)

sample2 <- c("EYFPneg_Shox2pos_CREpos","EYFPneg_Shox2neg_CREpos","EYFPpos_Shox2pos_CREpos", "EYFPpos_Shox2neg_CREpos")
n_cells2 <- c(EYFPneg_Shox2pos_CREpos,EYFPneg_Shox2neg_CREpos,EYFPpos_Shox2pos_CREpos,EYFPpos_Shox2neg_CREpos)

sample3 <- c("EYFPneg_Shox2pos_CREneg","EYFPneg_Shox2neg_CREneg","EYFPpos_Shox2pos_CREneg", "EYFPpos_Shox2neg_CREneg")
n_cells3 <- c(EYFPneg_Shox2pos_CREneg,EYFPneg_Shox2neg_CREneg,EYFPpos_Shox2pos_CREneg,EYFPpos_Shox2neg_CREneg)

df2_HL_E115_rep1.singlets <- data.frame(sample1,n_cells1,sample2,n_cells2,sample3,n_cells3)
write.csv(df2_HL_E115_rep1.singlets, file = "df2_HL_E115_rep2.singlets.csv")

d1 <- data.frame(EYFPpos_Shox2pos_CREpos_cells, "EYFPpos_Shox2pos_CREpos_cells", "A")
colnames(d1) <- c("cells","gene_class","letter")
d2 <- data.frame(EYFPneg_Shox2pos_CREpos_cells, "EYFPneg_Shox2pos_CREpos_cells", "B")
colnames(d2) <- c("cells","gene_class","letter")
d3 <- data.frame(EYFPpos_Shox2neg_CREpos_cells, "EYFPpos_Shox2neg_CREpos_cells", "C")
colnames(d3) <- c("cells","gene_class","letter")
d4 <- data.frame(EYFPneg_Shox2neg_CREpos_cells, "EYFPneg_Shox2neg_CREpos_cells", "D")
colnames(d4) <- c("cells","gene_class","letter")
d5 <- data.frame(EYFPpos_Shox2pos_CREneg_cells, "EYFPpos_Shox2pos_CREneg_cells", "E")
colnames(d5) <- c("cells","gene_class","letter")
d6 <- data.frame(EYFPneg_Shox2pos_CREneg_cells, "EYFPneg_Shox2pos_CREneg_cells", "F")
colnames(d6) <- c("cells","gene_class","letter")
d7 <- data.frame(EYFPpos_Shox2neg_CREneg_cells, "EYFPpos_Shox2neg_CREneg_cells", "G")
colnames(d7) <- c("cells","gene_class","letter")
d8 <- data.frame(EYFPneg_Shox2neg_CREneg_cells, "EYFPneg_Shox2neg_CREneg_cells", "H")
colnames(d8) <- c("cells","gene_class","letter")

cells_class <- rbind(d1,d2,d3,d4,d5,d6,d7,d8)
row.names(x = cells_class) <- cells_class$cells
cells_class2 <- cells_class[,c("gene_class","letter")]
head(cells_class2)

HL_E115_rep1.singlets[[colnames(x = cells_class2)]] <- cells_class2
head(HL_E115_rep1.singlets@meta.data)

pdf(file=paste0("HL_E115_rep1.singlets_CRE_dsmCherry_Shox2_EYFP_proportions.pdf"))

DimPlot(HL_E115_rep1.singlets, reduction = "umap",label = FALSE, group.by = "gene_class")
VlnPlot(object = HL_E115_rep1.singlets, assay = "RNA", features = "Shox2", group.by = "gene_class",
        pt.size = 0.1) 
VlnPlot(object = HL_E115_rep1.singlets, assay = "RNA", features = "CRE", group.by = "gene_class",
        pt.size = 0.1) 
VlnPlot(object = HL_E115_rep1.singlets, assay = "RNA", features = "EYFP", group.by = "gene_class",
        pt.size = 0.1) 

t1 <- table(HL_E115_rep1.singlets@meta.data$gene_class)
total <- sum(t1)
total
t1 <- as.data.frame(t1)
t1
t1$ratio <- t1[,2]/total
sum(t1$ratio)

ggplot(t1, aes(fill= Var1, x=Var1,y=Freq)) + 
  geom_bar(position="dodge", stat = "identity") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1))+
  geom_text(aes(label=Freq), vjust=-0.3, color="black", size=3.5)

ggplot(t1, aes(fill= Var1, x=Var1,y=ratio)) + 
  geom_bar(position="dodge", stat = "identity") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1))+
  geom_text(aes(label = sprintf("%0.2f", round(ratio, 2))), vjust=-0.3, color="black", size=3.5)
dev.off()

# save as seurat object 
saveRDS(HL_E115_rep1.singlets, file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E115.singlets.rds")

#### C.3 -> HL_E125.singlets ####
HL_E125_rep1.singlets <- readRDS("/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E125_singlets.rds")
HL_E125_rep1.singlets

# flag and count CRE, dsmCherry, EYFP & Shox2 expressing cells
CRE_counts <- FetchData(HL_E125_rep1.singlets, vars="CRE", slot="counts")
dsmCherry_counts <- FetchData(HL_E125_rep1.singlets, vars="dsmCherry", slot="counts")
EYFP_counts <- FetchData(HL_E125_rep1.singlets, vars="EYFP", slot="counts")
Shox2_counts <- FetchData(HL_E125_rep1.singlets, vars="Shox2", slot="counts")

pdf(file=paste0("HL_E125_rep1.singlets_CRE_dsmCherry_Shox2_EYFP_counts_QC.pdf"))
ggplot(data=CRE_counts, aes(CRE_counts$`CRE`)) + geom_histogram(binwidth = 0.5) + labs(x="CRE")
ggplot(data=dsmCherry_counts, aes(dsmCherry_counts$`dsmCherry`)) + geom_histogram(binwidth = 0.5) + labs(x="dsmCherry")
ggplot(data=Shox2_counts, aes(Shox2_counts$`Shox2`)) + geom_histogram(binwidth = 0.5) + labs(x="Shox2")
ggplot(data=EYFP_counts, aes(EYFP_counts$`EYFP`)) + geom_histogram(binwidth = 0.5) + labs(x="EYFP")
ggplot(data=CRE_counts, aes(CRE_counts$`CRE`)) + xlim(0,20) + geom_histogram(binwidth = 0.5) + labs(x="CRE UMIs") 
ggplot(data=dsmCherry_counts, aes(dsmCherry_counts$`dsmCherry`)) + xlim(0,20) + geom_histogram(binwidth = 0.5) + labs(x="dsmCherry UMIs") 
ggplot(data=Shox2_counts, aes(Shox2_counts$`Shox2`)) + xlim(0,20) + geom_histogram(binwidth = 0.5) + labs(x="Shox2 UMIs") 
ggplot(data=EYFP_counts, aes(EYFP_counts$`EYFP`)) + xlim(0,20) + geom_histogram(binwidth = 0.5) + labs(x="EYFP UMIs") 
dev.off()

#CRE_counts
CRE_selected_cells <- WhichCells(HL_E125_rep1.singlets, expression = `CRE`>=1, slot="counts")
length(CRE_selected_cells)
HL_E125_rep1.singlets[["CRE_positive"]] <- colnames(HL_E125_rep1.singlets) %in% CRE_selected_cells

#dsmCherry_counts
dsmCherry_selected_cells <- WhichCells(HL_E125_rep1.singlets, expression = `dsmCherry`>=1, slot="counts")
length(dsmCherry_selected_cells)
HL_E125_rep1.singlets[["dsmCherry_positive"]] <- colnames(HL_E125_rep1.singlets) %in% dsmCherry_selected_cells

#EYFP_counts
EYFP_selected_cells <- WhichCells(HL_E125_rep1.singlets, expression = `EYFP`>=1, slot="counts")
length(EYFP_selected_cells)
HL_E125_rep1.singlets[["EYFP_positive"]] <- colnames(HL_E125_rep1.singlets) %in% EYFP_selected_cells

#Shox2_counts
Shox2_selected_cells <- WhichCells(HL_E125_rep1.singlets, expression = `Shox2`>=1, slot="counts")
length(Shox2_selected_cells)
HL_E125_rep1.singlets[["Shox2_positive"]] <- colnames(HL_E125_rep1.singlets) %in% Shox2_selected_cells

#Positive and negative distribution
Shox2_pos_cells <- WhichCells(HL_E125_rep1.singlets, expression = `Shox2`>=1, slot="counts")
Shox2_pos <- length(Shox2_pos_cells)
Shox2_neg_cells <- WhichCells(HL_E125_rep1.singlets, expression = `Shox2`<1, slot="counts")
Shox2_neg <- length(Shox2_neg_cells)
EYFP_pos_cells <- WhichCells(HL_E125_rep1.singlets, expression = `EYFP`>=1, slot="counts")
EYFP_pos <- length(EYFP_pos_cells)
EYFP_neg_cells <- WhichCells(HL_E125_rep1.singlets, expression = `EYFP`<1, slot="counts")
EYFP_neg <- length(EYFP_neg_cells)
CRE_pos_cells <- WhichCells(HL_E125_rep1.singlets, expression = `CRE`>=1, slot="counts")
CRE_pos <- length(CRE_pos_cells)
CRE_neg_cells <- WhichCells(HL_E125_rep1.singlets, expression = `CRE`<1, slot="counts")
CRE_neg <- length(CRE_neg_cells)

Name <- c("Shox2","EYFP","CRE")
Positive <- c(Shox2_pos,EYFP_pos,CRE_pos)
Negative <- c(Shox2_neg,EYFP_neg,CRE_neg)
df1_HL_E125_rep1.singlets <- data.frame(Name,Positive,Negative)
write.csv(df1_HL_E125_rep1.singlets, file = "df1_HL_E125_rep1.singlets.csv")

EYFPneg_Shox2pos_cells <- WhichCells(HL_E125_rep1.singlets, cells = Shox2_pos_cells , expression = `EYFP` <1, slot="counts")
EYFPneg_Shox2pos <- length(EYFPneg_Shox2pos_cells)

EYFPneg_Shox2neg_cells <- WhichCells(HL_E125_rep1.singlets, cells = Shox2_neg_cells , expression = `EYFP` <1, slot="counts")
EYFPneg_Shox2neg <- length(EYFPneg_Shox2neg_cells)

EYFPpos_Shox2pos_cells <- WhichCells(HL_E125_rep1.singlets, cells = Shox2_pos_cells , expression = `EYFP` >=1, slot="counts")
EYFPpos_Shox2pos <- length(EYFPpos_Shox2pos_cells)

EYFPpos_Shox2neg_cells <- WhichCells(HL_E125_rep1.singlets, cells = Shox2_neg_cells , expression = `EYFP` >=1, slot="counts")
EYFPpos_Shox2neg <- length(EYFPpos_Shox2neg_cells)

EYFPpos_Shox2pos_CREpos_cells <- WhichCells(HL_E125_rep1.singlets, cells = EYFPpos_Shox2pos_cells, expression = `CRE` >=1, slot="counts")
EYFPpos_Shox2pos_CREpos <- length(EYFPpos_Shox2pos_CREpos_cells)

EYFPneg_Shox2pos_CREpos_cells <- WhichCells(HL_E125_rep1.singlets, cells = EYFPneg_Shox2pos_cells, expression = `CRE` >=1, slot="counts")
EYFPneg_Shox2pos_CREpos <- length(EYFPneg_Shox2pos_CREpos_cells)

EYFPpos_Shox2neg_CREpos_cells <- WhichCells(HL_E125_rep1.singlets, cells = EYFPpos_Shox2neg_cells, expression = `CRE` >=1, slot="counts")
EYFPpos_Shox2neg_CREpos <- length(EYFPpos_Shox2neg_CREpos_cells)

EYFPneg_Shox2neg_CREpos_cells <- WhichCells(HL_E125_rep1.singlets, cells = EYFPneg_Shox2neg_cells, expression = `CRE` >=1, slot="counts")
EYFPneg_Shox2neg_CREpos <- length(EYFPneg_Shox2neg_CREpos_cells)

EYFPpos_Shox2pos_CREneg_cells <- WhichCells(HL_E125_rep1.singlets, cells = EYFPpos_Shox2pos_cells, expression = `CRE` <1, slot="counts")
EYFPpos_Shox2pos_CREneg <- length(EYFPpos_Shox2pos_CREneg_cells)

EYFPneg_Shox2pos_CREneg_cells <- WhichCells(HL_E125_rep1.singlets, cells = EYFPneg_Shox2pos_cells, expression = `CRE` <1, slot="counts")
EYFPneg_Shox2pos_CREneg <- length(EYFPneg_Shox2pos_CREneg_cells)

EYFPpos_Shox2neg_CREneg_cells <- WhichCells(HL_E125_rep1.singlets, cells = EYFPpos_Shox2neg_cells, expression = `CRE` <1, slot="counts")
EYFPpos_Shox2neg_CREneg <- length(EYFPpos_Shox2neg_CREneg_cells)

EYFPneg_Shox2neg_CREneg_cells <- WhichCells(HL_E125_rep1.singlets, cells = EYFPneg_Shox2neg_cells, expression = `CRE` <1, slot="counts")
EYFPneg_Shox2neg_CREneg <- length(EYFPneg_Shox2neg_CREneg_cells)

sample1 <- c("EYFPneg_Shox2pos", "EYFPneg_Shox2neg","EYFPpos_Shox2pos","EYFPpos_Shox2neg")
n_cells1 <- c(EYFPneg_Shox2pos, EYFPneg_Shox2neg, EYFPpos_Shox2pos,EYFPpos_Shox2neg)

sample2 <- c("EYFPneg_Shox2pos_CREpos","EYFPneg_Shox2neg_CREpos","EYFPpos_Shox2pos_CREpos", "EYFPpos_Shox2neg_CREpos")
n_cells2 <- c(EYFPneg_Shox2pos_CREpos,EYFPneg_Shox2neg_CREpos,EYFPpos_Shox2pos_CREpos,EYFPpos_Shox2neg_CREpos)

sample3 <- c("EYFPneg_Shox2pos_CREneg","EYFPneg_Shox2neg_CREneg","EYFPpos_Shox2pos_CREneg", "EYFPpos_Shox2neg_CREneg")
n_cells3 <- c(EYFPneg_Shox2pos_CREneg,EYFPneg_Shox2neg_CREneg,EYFPpos_Shox2pos_CREneg,EYFPpos_Shox2neg_CREneg)

df2_HL_E125_rep1.singlets <- data.frame(sample1,n_cells1,sample2,n_cells2,sample3,n_cells3)
write.csv(df2_HL_E125_rep1.singlets, file = "df2_HL_E125_rep1.singlets.csv")

d1 <- data.frame(EYFPpos_Shox2pos_CREpos_cells, "EYFPpos_Shox2pos_CREpos_cells", "A")
colnames(d1) <- c("cells","gene_class","letter")
d2 <- data.frame(EYFPneg_Shox2pos_CREpos_cells, "EYFPneg_Shox2pos_CREpos_cells", "B")
colnames(d2) <- c("cells","gene_class","letter")
d3 <- data.frame(EYFPpos_Shox2neg_CREpos_cells, "EYFPpos_Shox2neg_CREpos_cells", "C")
colnames(d3) <- c("cells","gene_class","letter")
d4 <- data.frame(EYFPneg_Shox2neg_CREpos_cells, "EYFPneg_Shox2neg_CREpos_cells", "D")
colnames(d4) <- c("cells","gene_class","letter")
d5 <- data.frame(EYFPpos_Shox2pos_CREneg_cells, "EYFPpos_Shox2pos_CREneg_cells", "E")
colnames(d5) <- c("cells","gene_class","letter")
d6 <- data.frame(EYFPneg_Shox2pos_CREneg_cells, "EYFPneg_Shox2pos_CREneg_cells", "F")
colnames(d6) <- c("cells","gene_class","letter")
d7 <- data.frame(EYFPpos_Shox2neg_CREneg_cells, "EYFPpos_Shox2neg_CREneg_cells", "G")
colnames(d7) <- c("cells","gene_class","letter")
d8 <- data.frame(EYFPneg_Shox2neg_CREneg_cells, "EYFPneg_Shox2neg_CREneg_cells", "H")
colnames(d8) <- c("cells","gene_class","letter")

cells_class <- rbind(d1,d2,d3,d4,d5,d6,d7,d8)
row.names(x = cells_class) <- cells_class$cells
cells_class2 <- cells_class[,c("gene_class","letter")]
head(cells_class2)

HL_E125_rep1.singlets[[colnames(x = cells_class2)]] <- cells_class2
head(HL_E125_rep1.singlets@meta.data)

pdf(file=paste0("HL_E125_rep1.singlets_CRE_dsmCherry_Shox2_EYFP_proportions.pdf"))

DimPlot(HL_E125_rep1.singlets, reduction = "umap",label = FALSE, group.by = "gene_class")
VlnPlot(object = HL_E125_rep1.singlets, assay = "RNA", features = "Shox2", group.by = "gene_class",
        pt.size = 0.1) 
VlnPlot(object = HL_E125_rep1.singlets, assay = "RNA", features = "CRE", group.by = "gene_class",
        pt.size = 0.1) 
VlnPlot(object = HL_E125_rep1.singlets, assay = "RNA", features = "EYFP", group.by = "gene_class",
        pt.size = 0.1) 

t1 <- table(HL_E125_rep1.singlets@meta.data$gene_class)
total <- sum(t1)
total
t1 <- as.data.frame(t1)
t1
t1$ratio <- t1[,2]/total
sum(t1$ratio)

ggplot(t1, aes(fill= Var1, x=Var1,y=Freq)) + 
  geom_bar(position="dodge", stat = "identity") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1))+
  geom_text(aes(label=Freq), vjust=-0.3, color="black", size=3.5)

ggplot(t1, aes(fill= Var1, x=Var1,y=ratio)) + 
  geom_bar(position="dodge", stat = "identity") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1))+
  geom_text(aes(label = sprintf("%0.2f", round(ratio, 2))), vjust=-0.3, color="black", size=3.5)
dev.off()

# save as seurat object 
saveRDS(HL_E125_rep1.singlets, file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E125.singlets.rds")

#### C.4 -> HL_E135.singlets ####
HL_E135_rep1.singlets <- readRDS("/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E135_singlets.rds")
HL_E135_rep1.singlets

# flag and count CRE, dsmCherry, EYFP & Shox2 expressing cells
CRE_counts <- FetchData(HL_E135_rep1.singlets, vars="CRE", slot="counts")
dsmCherry_counts <- FetchData(HL_E135_rep1.singlets, vars="dsmCherry", slot="counts")
EYFP_counts <- FetchData(HL_E135_rep1.singlets, vars="EYFP", slot="counts")
Shox2_counts <- FetchData(HL_E135_rep1.singlets, vars="Shox2", slot="counts")

pdf(file=paste0("HL_E135_rep1.singlets_CRE_dsmCherry_Shox2_EYFP_counts_QC.pdf"))
ggplot(data=CRE_counts, aes(CRE_counts$`CRE`)) + geom_histogram(binwidth = 0.5) + labs(x="CRE")
ggplot(data=dsmCherry_counts, aes(dsmCherry_counts$`dsmCherry`)) + geom_histogram(binwidth = 0.5) + labs(x="dsmCherry")
ggplot(data=Shox2_counts, aes(Shox2_counts$`Shox2`)) + geom_histogram(binwidth = 0.5) + labs(x="Shox2")
ggplot(data=EYFP_counts, aes(EYFP_counts$`EYFP`)) + geom_histogram(binwidth = 0.5) + labs(x="EYFP")
ggplot(data=CRE_counts, aes(CRE_counts$`CRE`)) + xlim(0,20) + geom_histogram(binwidth = 0.5) + labs(x="CRE UMIs") 
ggplot(data=dsmCherry_counts, aes(dsmCherry_counts$`dsmCherry`)) + xlim(0,20) + geom_histogram(binwidth = 0.5) + labs(x="dsmCherry UMIs") 
ggplot(data=Shox2_counts, aes(Shox2_counts$`Shox2`)) + xlim(0,20) + geom_histogram(binwidth = 0.5) + labs(x="Shox2 UMIs") 
ggplot(data=EYFP_counts, aes(EYFP_counts$`EYFP`)) + xlim(0,20) + geom_histogram(binwidth = 0.5) + labs(x="EYFP UMIs") 
dev.off()

#CRE_counts
CRE_selected_cells <- WhichCells(HL_E135_rep1.singlets, expression = `CRE`>=1, slot="counts")
length(CRE_selected_cells)
HL_E135_rep1.singlets[["CRE_positive"]] <- colnames(HL_E135_rep1.singlets) %in% CRE_selected_cells

#dsmCherry_counts
dsmCherry_selected_cells <- WhichCells(HL_E135_rep1.singlets, expression = `dsmCherry`>=1, slot="counts")
length(dsmCherry_selected_cells)
HL_E135_rep1.singlets[["dsmCherry_positive"]] <- colnames(HL_E135_rep1.singlets) %in% dsmCherry_selected_cells

#EYFP_counts
EYFP_selected_cells <- WhichCells(HL_E135_rep1.singlets, expression = `EYFP`>=1, slot="counts")
length(EYFP_selected_cells)
HL_E135_rep1.singlets[["EYFP_positive"]] <- colnames(HL_E135_rep1.singlets) %in% EYFP_selected_cells

#Shox2_counts
Shox2_selected_cells <- WhichCells(HL_E135_rep1.singlets, expression = `Shox2`>=1, slot="counts")
length(Shox2_selected_cells)
HL_E135_rep1.singlets[["Shox2_positive"]] <- colnames(HL_E135_rep1.singlets) %in% Shox2_selected_cells

#Positive and negative distribution
Shox2_pos_cells <- WhichCells(HL_E135_rep1.singlets, expression = `Shox2`>=1, slot="counts")
Shox2_pos <- length(Shox2_pos_cells)
Shox2_neg_cells <- WhichCells(HL_E135_rep1.singlets, expression = `Shox2`<1, slot="counts")
Shox2_neg <- length(Shox2_neg_cells)
EYFP_pos_cells <- WhichCells(HL_E135_rep1.singlets, expression = `EYFP`>=1, slot="counts")
EYFP_pos <- length(EYFP_pos_cells)
EYFP_neg_cells <- WhichCells(HL_E135_rep1.singlets, expression = `EYFP`<1, slot="counts")
EYFP_neg <- length(EYFP_neg_cells)
CRE_pos_cells <- WhichCells(HL_E135_rep1.singlets, expression = `CRE`>=1, slot="counts")
CRE_pos <- length(CRE_pos_cells)
CRE_neg_cells <- WhichCells(HL_E135_rep1.singlets, expression = `CRE`<1, slot="counts")
CRE_neg <- length(CRE_neg_cells)

Name <- c("Shox2","EYFP","CRE")
Positive <- c(Shox2_pos,EYFP_pos,CRE_pos)
Negative <- c(Shox2_neg,EYFP_neg,CRE_neg)
df1_HL_E135_rep1.singlets <- data.frame(Name,Positive,Negative)
write.csv(df1_HL_E135_rep1.singlets, file = "df1_HL_E135_rep1.singlets.csv")

EYFPneg_Shox2pos_cells <- WhichCells(HL_E135_rep1.singlets, cells = Shox2_pos_cells , expression = `EYFP` <1, slot="counts")
EYFPneg_Shox2pos <- length(EYFPneg_Shox2pos_cells)

EYFPneg_Shox2neg_cells <- WhichCells(HL_E135_rep1.singlets, cells = Shox2_neg_cells , expression = `EYFP` <1, slot="counts")
EYFPneg_Shox2neg <- length(EYFPneg_Shox2neg_cells)

EYFPpos_Shox2pos_cells <- WhichCells(HL_E135_rep1.singlets, cells = Shox2_pos_cells , expression = `EYFP` >=1, slot="counts")
EYFPpos_Shox2pos <- length(EYFPpos_Shox2pos_cells)

EYFPpos_Shox2neg_cells <- WhichCells(HL_E135_rep1.singlets, cells = Shox2_neg_cells , expression = `EYFP` >=1, slot="counts")
EYFPpos_Shox2neg <- length(EYFPpos_Shox2neg_cells)

EYFPpos_Shox2pos_CREpos_cells <- WhichCells(HL_E135_rep1.singlets, cells = EYFPpos_Shox2pos_cells, expression = `CRE` >=1, slot="counts")
EYFPpos_Shox2pos_CREpos <- length(EYFPpos_Shox2pos_CREpos_cells)

EYFPneg_Shox2pos_CREpos_cells <- WhichCells(HL_E135_rep1.singlets, cells = EYFPneg_Shox2pos_cells, expression = `CRE` >=1, slot="counts")
EYFPneg_Shox2pos_CREpos <- length(EYFPneg_Shox2pos_CREpos_cells)

EYFPpos_Shox2neg_CREpos_cells <- WhichCells(HL_E135_rep1.singlets, cells = EYFPpos_Shox2neg_cells, expression = `CRE` >=1, slot="counts")
EYFPpos_Shox2neg_CREpos <- length(EYFPpos_Shox2neg_CREpos_cells)

EYFPneg_Shox2neg_CREpos_cells <- WhichCells(HL_E135_rep1.singlets, cells = EYFPneg_Shox2neg_cells, expression = `CRE` >=1, slot="counts")
EYFPneg_Shox2neg_CREpos <- length(EYFPneg_Shox2neg_CREpos_cells)

EYFPpos_Shox2pos_CREneg_cells <- WhichCells(HL_E135_rep1.singlets, cells = EYFPpos_Shox2pos_cells, expression = `CRE` <1, slot="counts")
EYFPpos_Shox2pos_CREneg <- length(EYFPpos_Shox2pos_CREneg_cells)

EYFPneg_Shox2pos_CREneg_cells <- WhichCells(HL_E135_rep1.singlets, cells = EYFPneg_Shox2pos_cells, expression = `CRE` <1, slot="counts")
EYFPneg_Shox2pos_CREneg <- length(EYFPneg_Shox2pos_CREneg_cells)

EYFPpos_Shox2neg_CREneg_cells <- WhichCells(HL_E135_rep1.singlets, cells = EYFPpos_Shox2neg_cells, expression = `CRE` <1, slot="counts")
EYFPpos_Shox2neg_CREneg <- length(EYFPpos_Shox2neg_CREneg_cells)

EYFPneg_Shox2neg_CREneg_cells <- WhichCells(HL_E135_rep1.singlets, cells = EYFPneg_Shox2neg_cells, expression = `CRE` <1, slot="counts")
EYFPneg_Shox2neg_CREneg <- length(EYFPneg_Shox2neg_CREneg_cells)

sample1 <- c("EYFPneg_Shox2pos", "EYFPneg_Shox2neg","EYFPpos_Shox2pos","EYFPpos_Shox2neg")
n_cells1 <- c(EYFPneg_Shox2pos, EYFPneg_Shox2neg, EYFPpos_Shox2pos,EYFPpos_Shox2neg)

sample2 <- c("EYFPneg_Shox2pos_CREpos","EYFPneg_Shox2neg_CREpos","EYFPpos_Shox2pos_CREpos", "EYFPpos_Shox2neg_CREpos")
n_cells2 <- c(EYFPneg_Shox2pos_CREpos,EYFPneg_Shox2neg_CREpos,EYFPpos_Shox2pos_CREpos,EYFPpos_Shox2neg_CREpos)

sample3 <- c("EYFPneg_Shox2pos_CREneg","EYFPneg_Shox2neg_CREneg","EYFPpos_Shox2pos_CREneg", "EYFPpos_Shox2neg_CREneg")
n_cells3 <- c(EYFPneg_Shox2pos_CREneg,EYFPneg_Shox2neg_CREneg,EYFPpos_Shox2pos_CREneg,EYFPpos_Shox2neg_CREneg)

df2_HL_E135_rep1.singlets <- data.frame(sample1,n_cells1,sample2,n_cells2,sample3,n_cells3)
write.csv(df2_HL_E135_rep1.singlets, file = "df2_HL_E135_rep1.singlets.csv")

d1 <- data.frame(EYFPpos_Shox2pos_CREpos_cells, "EYFPpos_Shox2pos_CREpos_cells", "A")
colnames(d1) <- c("cells","gene_class","letter")
d2 <- data.frame(EYFPneg_Shox2pos_CREpos_cells, "EYFPneg_Shox2pos_CREpos_cells", "B")
colnames(d2) <- c("cells","gene_class","letter")
d3 <- data.frame(EYFPpos_Shox2neg_CREpos_cells, "EYFPpos_Shox2neg_CREpos_cells", "C")
colnames(d3) <- c("cells","gene_class","letter")
d4 <- data.frame(EYFPneg_Shox2neg_CREpos_cells, "EYFPneg_Shox2neg_CREpos_cells", "D")
colnames(d4) <- c("cells","gene_class","letter")
d5 <- data.frame(EYFPpos_Shox2pos_CREneg_cells, "EYFPpos_Shox2pos_CREneg_cells", "E")
colnames(d5) <- c("cells","gene_class","letter")
d6 <- data.frame(EYFPneg_Shox2pos_CREneg_cells, "EYFPneg_Shox2pos_CREneg_cells", "F")
colnames(d6) <- c("cells","gene_class","letter")
d7 <- data.frame(EYFPpos_Shox2neg_CREneg_cells, "EYFPpos_Shox2neg_CREneg_cells", "G")
colnames(d7) <- c("cells","gene_class","letter")
d8 <- data.frame(EYFPneg_Shox2neg_CREneg_cells, "EYFPneg_Shox2neg_CREneg_cells", "H")
colnames(d8) <- c("cells","gene_class","letter")

cells_class <- rbind(d1,d2,d3,d4,d5,d6,d7,d8)
row.names(x = cells_class) <- cells_class$cells
cells_class2 <- cells_class[,c("gene_class","letter")]
head(cells_class2)

HL_E135_rep1.singlets[[colnames(x = cells_class2)]] <- cells_class2
head(HL_E135_rep1.singlets@meta.data)

pdf(file=paste0("HL_E135_rep1.singlets_CRE_dsmCherry_Shox2_EYFP_proportions.pdf"))

DimPlot(HL_E135_rep1.singlets, reduction = "umap",label = FALSE, group.by = "gene_class")
VlnPlot(object = HL_E135_rep1.singlets, assay = "RNA", features = "Shox2", group.by = "gene_class",
        pt.size = 0.1) 
VlnPlot(object = HL_E135_rep1.singlets, assay = "RNA", features = "CRE", group.by = "gene_class",
        pt.size = 0.1) 
VlnPlot(object = HL_E135_rep1.singlets, assay = "RNA", features = "EYFP", group.by = "gene_class",
        pt.size = 0.1) 

t1 <- table(HL_E135_rep1.singlets@meta.data$gene_class)
total <- sum(t1)
total
t1 <- as.data.frame(t1)
t1
t1$ratio <- t1[,2]/total
sum(t1$ratio)

ggplot(t1, aes(fill= Var1, x=Var1,y=Freq)) + 
  geom_bar(position="dodge", stat = "identity") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1))+
  geom_text(aes(label=Freq), vjust=-0.3, color="black", size=3.5)

ggplot(t1, aes(fill= Var1, x=Var1,y=ratio)) + 
  geom_bar(position="dodge", stat = "identity") + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust=1))+
  geom_text(aes(label = sprintf("%0.2f", round(ratio, 2))), vjust=-0.3, color="black", size=3.5)
dev.off()

# save as seurat object 
saveRDS(HL_E135_rep1.singlets, file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E135.singlets.rds")

#### D Create a single MERGED Seurat Object #####
HL_E105.singlets <- readRDS("/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E105.singlets.rds")
HL_E105.singlets
summary(HL_E105.singlets@meta.data)
table(HL_E105.singlets@meta.data$gene_class)

HL_E115.singlets <- readRDS("/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E115.singlets.rds")
HL_E115.singlets
summary(HL_E115.singlets@meta.data)
table(HL_E115.singlets@meta.data$gene_class)

HL_E125.singlets <- readRDS("/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E125.singlets.rds")
HL_E125.singlets
summary(HL_E125.singlets@meta.data)
table(HL_E125.singlets@meta.data$gene_class)

HL_E135.singlets <- readRDS("/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/HL_E135.singlets.rds")
HL_E135.singlets
summary(HL_E135.singlets@meta.data)
table(HL_E135.singlets@meta.data$gene_class)


all.merged <- merge(HL_E105.singlets, y= c(HL_E115.singlets,HL_E125.singlets,HL_E135.singlets),
                    add.cell.ids = c("HL_E105","HL_E115","HL_E125","HL_E135"))

all.merged
summary(all.merged@meta.data)
table(all.merged$orig.ident)
table(all.merged$gene_class)

## add $stage @metadata ##
all.merged@meta.data$stage <- gsub("[0-9]$","", all.merged$orig.ident)
all.merged@meta.data$stage <- gsub("_rep","", all.merged$stage)
all.merged@meta.data$stage <- gsub("HL_","", all.merged$stage)
table(all.merged$stage)

saveRDS(all.merged, file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/all.merged.rds")

#### E. SCT Transform prior to CC regression ####
all.merged <- readRDS("/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/all.merged.rds")
all.merged
summary(all.merged@meta.data)
table(all.merged$orig.ident)
table(all.merged$gene_class)
table(all.merged$stage)

DefaultAssay(object = all.merged) <- "spliced"
all.merged

#sct
all.merged <- SCTransform(all.merged, assay = "spliced")
all.merged
length(VariableFeatures(all.merged))

# excluding the reporter CRE from the variable genes to avoid driving the PCA
grep ("CRE", all.merged@assays[["SCT"]]@var.features)
all.merged@assays[["SCT"]]@var.features <- grep ("CRE", all.merged@assays[["SCT"]]@var.features, invert = TRUE, value=TRUE)

# excluding the reporter EYFP from the variable genes to avoid driving the PCA
grep ("EYFP", all.merged@assays[["SCT"]]@var.features)
all.merged@assays[["SCT"]]@var.features <- grep ("EYFP", all.merged@assays[["SCT"]]@var.features, invert = TRUE, value=TRUE)

# excluding the reporter dsmCherry from the variable genes to avoid driving the PCA
grep ("dsmCherry", all.merged@assays[["SCT"]]@var.features)
all.merged@assays[["SCT"]]@var.features <- grep ("dsmCherry", all.merged@assays[["SCT"]]@var.features, invert = TRUE, value=TRUE)

length(VariableFeatures(all.merged))

# Check cell-cycle
s.genes <- str_to_title(cc.genes$s.genes)
g2m.genes <- str_to_title(cc.genes$g2m.genes)
all.merged <- CellCycleScoring(all.merged, s.features = s.genes, g2m.features = g2m.genes, assay="SCT", set.ident = FALSE)

pdf(file=paste0("all.merged_CellCycle.pdf"))
Idents(all.merged) <- "Phase"
RidgePlot(all.merged, features = c("Pcna", "Top2a", "Mcm6", "Mki67"), ncol = 2)
dev.off()

Idents(all.merged) <- "orig.ident"
Idents(all.merged)
summary(all.merged@meta.data)

saveRDS(all.merged, file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/all.merged.beforeSCTCCregression.rds")

#### F. SCTransform, regressing cell-cycle and stage and excluding the EYFP, CRE, dsmCherry from variable genes ####
all.merged <- readRDS("/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/all.merged.beforeSCTCCregression.rds")
all.merged
summary(all.merged@meta.data)

DefaultAssay(object = all.merged) <- "spliced"
all.merged

all.merged.SCTCC <- SCTransform(all.merged,vars.to.regress = c("S.Score", "G2M.Score","stage"),assay = "spliced")
all.merged.SCTCC
summary(all.merged.SCTCC@meta.data)

# excluding the reporter CRE from the variable genes to avoid driving the PCA
grep ("CRE", all.merged.SCTCC@assays[["SCT"]]@var.features)
all.merged.SCTCC@assays[["SCT"]]@var.features <- grep ("CRE", all.merged.SCTCC@assays[["SCT"]]@var.features, invert = TRUE, value=TRUE)

# excluding the reporter EYFP from the variable genes to avoid driving the PCA
grep ("EYFP", all.merged.SCTCC@assays[["SCT"]]@var.features)
all.merged.SCTCC@assays[["SCT"]]@var.features <- grep ("EYFP", all.merged.SCTCC@assays[["SCT"]]@var.features, invert = TRUE, value=TRUE)

# excluding the reporter dsmCherry from the variable genes to avoid driving the PCA
grep ("dsmCherry", all.merged.SCTCC@assays[["SCT"]]@var.features)
all.merged.SCTCC@assays[["SCT"]]@var.features <- grep ("dsmCherry", all.merged.SCTCC@assays[["SCT"]]@var.features, invert = TRUE, value=TRUE)

length(VariableFeatures(all.merged.SCTCC))

saveRDS(all.merged.SCTCC, file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.merged.SCTCCstage.rds")

#### G.1 Final PCA ####
all.merged.SCTCCSt <- readRDS("/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.merged.SCTCCstage.rds")
all.merged.SCTCCSt
summary(all.merged.SCTCCSt@meta.data)
table(all.merged.SCTCCSt$gene_class)
table(all.merged.SCTCCSt$stage)

#PCA on the SCT with CC regression and stage
all.merged.SCTCCSt <- RunPCA(object = all.merged.SCTCCSt, npcs = 50)

pdf(file=paste0("all.merged.SCTCCSt_PCA.pdf"))
DimPlot(all.merged.SCTCCSt, reduction = "pca")
Idents(all.merged.SCTCCSt) <- "Phase"
DimPlot(all.merged.SCTCCSt, reduction = "pca")
Idents(all.merged.SCTCCSt) <- "stage"
Idents(all.merged.SCTCCSt)
dev.off()

saveRDS(all.merged.SCTCCSt, "/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.merged.SCTCCSt.PCA.rds")

##PCA on the SCT prior CC regression to compare
all.merged <- readRDS("/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.merged.SCTCCstage.rds")
all.merged
all.merged <- RunPCA(object = all.merged, npcs = 50)

pdf(file=paste0("all.merged_PCA.pdf"))
DimPlot(all.merged, reduction = "pca")
Idents(all.merged) <- "Phase"
DimPlot(all.merged, reduction = "pca")
Idents(all.merged) <- "orig.ident"
Idents(all.merged)
dev.off()

##### G.2 Final UMAP ####
all.merged.SCTCCSt <- readRDS("/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.merged.SCTCCSt.PCA.rds")
all.merged.SCTCCSt

#UMAP on the SCT
all.merged.SCTCCSt <- RunUMAP(object = all.merged.SCTCCSt, dims = 1:50, n.neighbors = 50L, metric = "cosine",min.dist = 0.2)
Idents(all.merged.SCTCCSt) <- "stage"
Idents(all.merged.SCTCCSt)
DimPlot(all.merged.SCTCCSt, reduction = "umap")

pdf(file=paste0("all.merged.SCTCCSt_UMAP.pdf"))
DimPlot(all.merged.SCTCCSt, reduction = "umap")
Idents(all.merged.SCTCCSt) <- "Phase"
Idents(all.merged.SCTCCSt)
DimPlot(all.merged.SCTCCSt, reduction = "umap")
Idents(all.merged.SCTCCSt) <- "gene_class"
Idents(all.merged.SCTCCSt)
DimPlot(all.merged.SCTCCSt, reduction = "umap")
DimPlot(all.merged.SCTCCSt, reduction = "umap", split.by = "stage")
DimPlot(all.merged.SCTCCSt, reduction = "umap", split.by = "gene_class")
DimPlot(all.merged.SCTCCSt, reduction = "umap", split.by = "Phase")
Idents(all.merged.SCTCCSt) <- "stage"
Idents(all.merged.SCTCCSt)
DefaultAssay(object = all.merged.SCTCCSt) <- "RNA"
all.merged.SCTCCSt
genes <- c("Pitx1","Shox2","Hoxd13","CRE","dsmCherry","EYFP","Wnt6","Myod1","Krt14","Cdh5","Hba-a1")
FeaturePlot(object = all.merged.SCTCCSt, features = genes ,cols = c("grey95","darkred"), combine = FALSE)
genes <- c("Dcn","Gdf5","Shh","Ihh","Sox9","Irx5","Msx1","Osr1","Ccr1")
FeaturePlot(object = all.merged.SCTCCSt, features = genes ,cols = c("grey95","darkred"), combine = FALSE)
genes2 <- c("Shox2","Hoxd13")
FeaturePlot(object = all.merged.SCTCCSt, features = genes2 ,cols = c("grey95","darkred"), combine = FALSE, split.by ="gene_class" )
DefaultAssay(object = all.merged.SCTCCSt) <- "SCT"
all.merged.SCTCCSt
dev.off()

#Find Neighbors on the SCT
all.merged.SCTCCSt <- FindNeighbors(object = all.merged.SCTCCSt)

# save as seurat object 
saveRDS(all.merged.SCTCCSt, file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.merged.SCTCCSt.UMAP.rds")

#### H.1 Trial Clustering ####
all.merged.SCTCCSt <- readRDS("/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.merged.SCTCCSt.UMAP.rds")
all.merged.SCTCCSt
summary(all.merged.SCTCCSt@meta.data)

all.merged.SCTCCSt <- FindClusters(object = all.merged.SCTCCSt, resolution = seq(from = 0.1, to = 2, by = 0.1))
clustree(x = all.merged.SCTCCSt, prefix = "SCT_snn_res.")
summary(all.merged.SCTCCSt@meta.data)
pdf(file=paste0("all.merged.STCC_UMAP_clustering_options.pdf"))
DimPlot(object = all.merged.SCTCCSt, reduction = "umap", group.by = "SCT_snn_res.0.1", label = TRUE, label.size = 6)
DimPlot(object = all.merged.SCTCCSt, reduction = "umap", group.by = "SCT_snn_res.0.2", label = TRUE, label.size = 6) 
DimPlot(object = all.merged.SCTCCSt, reduction = "umap", group.by = "SCT_snn_res.0.3", label = TRUE, label.size = 6)
DimPlot(object = all.merged.SCTCCSt, reduction = "umap", group.by = "SCT_snn_res.0.4", label = TRUE, label.size = 6)
DimPlot(object = all.merged.SCTCCSt, reduction = "umap", group.by = "SCT_snn_res.0.5", label = TRUE, label.size = 6)
DimPlot(object = all.merged.SCTCCSt, reduction = "umap", group.by = "SCT_snn_res.0.6", label = TRUE, label.size = 6)
DimPlot(object = all.merged.SCTCCSt, reduction = "umap", group.by = "SCT_snn_res.0.7", label = TRUE, label.size = 6)
Idents(all.merged.SCTCCSt) <- "stage"
Idents(all.merged.SCTCCSt)
DimPlot(all.merged.SCTCCSt, reduction = "umap")
Idents(all.merged.SCTCCSt) <- "Phase"
Idents(all.merged.SCTCCSt)
DimPlot(all.merged.SCTCCSt, reduction = "umap")
FeaturePlot(object = all.merged.SCTCCSt, features = "Cdh5" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.merged.SCTCCSt, features = "Cldn5" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.merged.SCTCCSt, features = "Wnt6" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.merged.SCTCCSt, features = "Krt14" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.merged.SCTCCSt, features = "Cd52" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.merged.SCTCCSt, features = "Ttn" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.merged.SCTCCSt, features = "C1qa" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.merged.SCTCCSt, features = "Ccr1" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.merged.SCTCCSt, features = "Hba-a1" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.merged.SCTCCSt, features = "Pitx1" ,cols = c("grey95","darkred"), combine = FALSE)
dev.off()

#### H.2 FINAL CLUSTERING ####
all.merged.SCTCCSt <- readRDS("/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.merged.SCTCCSt.UMAP.rds")
all.merged.SCTCCSt
summary(all.merged.SCTCCSt@meta.data)

all.merged.SCTCCSt <- FindClusters(object = all.merged.SCTCCSt, resolution = 0.7)

cell.num <- table(Idents(all.merged.SCTCCSt))
clusterLabels <- paste(names(cell.num), paste0("(", cell.num, ")"))
clusterBreaks <- names(cell.num)

pdf(file=paste0("all.merged.SCTCCSt_UMAP_allclusters.pdf"))
DimPlot(all.merged.SCTCCSt, reduction = "umap",label = TRUE, label.size = 5)
DimPlot(all.merged.SCTCCSt, reduction = "umap", label=TRUE) + scale_colour_discrete(breaks = clusterBreaks,labels = clusterLabels) + labs(title = "UMAP")
ptsize <- 0.5
DimPlot(all.merged.SCTCCSt, group.by = "CRE_positive", cols = c("grey95","darkred"), reduction = "umap", label=F, order=TRUE, pt.size=ptsize) + labs(title = "CRE_positive_cells") + NoLegend()
DimPlot(all.merged.SCTCCSt, group.by = "EYFP_positive", cols = c("grey95","darkred"), reduction = "umap", label=F, order=TRUE, pt.size=ptsize) + labs(title = "EYFP_positive_cells") + NoLegend()
DimPlot(all.merged.SCTCCSt, group.by = "Shox2_positive", cols = c("grey95","darkred"), reduction = "umap", label=F, order=TRUE, pt.size=ptsize) + labs(title = "Shox2_positive_cells") + NoLegend()
DimPlot(all.merged.SCTCCSt, group.by = "Phase", reduction = "umap", label=F)
DimPlot(all.merged.SCTCCSt, group.by = "gene_class", reduction = "umap", label=F)
dev.off()

# save as seurat object 
saveRDS(all.merged.SCTCCSt, file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.merged.SCTCCSt.clus.rds")

#### H.3 Merging clusters ####
#Since by checking the features markers expression I can already know which clusters corresponds to the sattelite and which ones are mesenchyme I'm gonna merge them before labelling
all.merged.SCTCCSt <- readRDS("/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.merged.SCTCCSt.clus.rds")
all.merged.SCTCCSt
Idents(all.merged.SCTCCSt) <- "seurat_clusters"
Idents(all.merged.SCTCCSt)
summary(all.merged.SCTCCSt@meta.data)
DimPlot(all.merged.SCTCCSt, reduction = "umap",label = TRUE, label.size = 6)
summary(all.merged.SCTCCSt@meta.data)

all.merged.SCTCCSt@meta.data$merged.clusters <- all.merged.SCTCCSt@meta.data$seurat_clusters
levels(all.merged.SCTCCSt$merged.clusters) <- c("Mesenchyme","Mesenchyme","Mesenchyme","Mesenchyme","Mesenchyme","Mesenchyme","Mesenchyme","Mesenchyme","Mesenchyme",
                                                "Epithelium","Mesenchyme","Mesenchyme","Mesenchyme",
                                                "Muscle","Endothelium","Mesenchyme","Muscle","Immune Cells",
                                                "Epithelium","Epithelium","Blood Cells")

DimPlot(all.merged.SCTCCSt, reduction = "umap",label = TRUE, label.size = 3, group.by = "merged.clusters", repel = TRUE)
DimPlot(all.merged.SCTCCSt, reduction = "umap",label = FALSE, group.by = "stage")
DimPlot(all.merged.SCTCCSt, reduction = "umap",label = FALSE, group.by = "merged.clusters", split.by = "stage")

levels(all.merged.SCTCCSt)

Idents(all.merged.SCTCCSt) <- "merged.clusters"

saveRDS(all.merged.SCTCCSt, "/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.merged.SCTCCSt.labelled.rds")

#### H.4 Marker genes per cluster ####
all.merged.SCTCCSt.labelled <- readRDS("/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.merged.SCTCCSt.labelled.rds")
all.merged.SCTCCSt.labelled
summary(all.merged.SCTCCSt.labelled@meta.data)
DimPlot(all.merged.SCTCCSt.labelled, reduction = "umap",label = TRUE, label.size = 3,repel=TRUE)
DimPlot(all.merged.SCTCCSt.labelled, reduction = "umap",label = TRUE, label.size = 3,repel=TRUE, split.by = "stage")

DefaultAssay(all.merged.SCTCCSt.labelled) <- "RNA"
all.merged.SCTCCSt.labelled

Mesenchyme_cluster.markers = FindMarkers(all.merged.SCTCCSt.labelled, ident.1 = "Mesenchyme", grouping.var = "stage", only.pos =TRUE)
head (Mesenchyme_cluster.markers,15)
write.csv(Mesenchyme_cluster.markers,file=("Mesenchyme_cluster.markers.csv"))

Epithelium_cluster.markers = FindMarkers(all.merged.SCTCCSt.labelled, ident.1 = "Epithelium", grouping.var = "stage", only.pos =TRUE)
head (Epithelium_cluster.markers,15)
write.csv(Epithelium_cluster.markers,file=("Epithelium_cluster.markers.csv"))

Endothelium_cluster.markers = FindMarkers(all.merged.SCTCCSt.labelled, ident.1 = "Endothelium", grouping.var = "stage", only.pos =TRUE)
head (Endothelium_cluster.markers,15)
write.csv(Endothelium_cluster.markers,file=("Endothelium_cluster.markers.csv"))

Muscle_cluster.markers = FindMarkers(all.merged.SCTCCSt.labelled, ident.1 ="Muscle", grouping.var = "stage", only.pos =TRUE)
head (Muscle_cluster.markers,15)
write.csv(Muscle_cluster.markers,file=("Muscle_cluster.markers.csv"))

Immune_cells.markers = FindMarkers(all.merged.SCTCCSt.labelled, ident.1 = "Immune Cells", grouping.var = "stage", only.pos =TRUE)
head (Immune_cells.markers,15)
write.csv(Immune_cells.markers,file=("Immune_cells_cluster.markers.csv"))

Blood_cells.markers = FindMarkers(all.merged.SCTCCSt.labelled, ident.1 = "Blood Cells", grouping.var = "stage", only.pos =TRUE)
head (Blood_cells.markers,15)
write.csv(Blood_cells.markers,file=("Blood_cells_cluster.markers.csv"))


#### I. Proportion analysis merged Mesenchyme and Epithelium clusters ####
all.merged.SCTCCSt.labelled <- readRDS("/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.merged.SCTCCSt.labelled.rds")
all.merged.SCTCCSt.labelled
summary(all.merged.SCTCCSt.labelled@meta.data)
DimPlot(all.merged.SCTCCSt.labelled, reduction = "umap",label = TRUE, label.size = 3)

#number of cells per cluster
table(all.merged.SCTCCSt.labelled$merged.clusters)

#proportion of cells per cluster
prop.table(table(Idents(all.merged.SCTCCSt.labelled)))

#number of cells per stage 
total_cells <- table(all.merged.SCTCCSt.labelled$stage)
total_cells

#number of cells per cluster and sample stage
number_cells_percluster <- table(all.merged.SCTCCSt.labelled$merged.clusters, all.merged.SCTCCSt.labelled$stage)
number_cells_percluster

## DPA = Differential Proportion Analysis ##
#Install functions needed for the analysis
source("~/Dropbox/DPA/diffprop_functions.R");

## Read in file of counts of cells in each population (Cluster) across conditions (datasets)
obs.counts <- t(number_cells_percluster)

## Run an example using error (p) of 0.1 and with 100,000 iterations
tip.exp <- generateNull(obs.counts, n=100000, p=0.1)    # Generate the null distribution based on sampling

### P-value tests for Sham vs MI-Day 3
pp <- two.class.test(obs.counts, tip.exp, cond.control="E105", cond.treatment="E115",to.plot=T)
pp
write.csv(pp,"UMAP_MERGED_E105vsE115-pvalues_all_clusters.csv")

pp1 <- two.class.test(obs.counts, tip.exp, cond.control="E115", cond.treatment="E125",to.plot=T)
pp1
write.csv(pp1,"UMAP_MERGED_E115vsE125-pvalues_all_clusters.csv")

pp2 <- two.class.test(obs.counts, tip.exp, cond.control="E125", cond.treatment="E135",to.plot=T)
pp2
write.csv(pp2,"UMAP_MERGED_E125vsE135-pvalues_all_clusters.csv")

pp3 <- two.class.test(obs.counts, tip.exp, cond.control="E105", cond.treatment="E125",to.plot=T)
pp3
write.csv(pp3,"UMAP_MERGED_E105vsE125-pvalues_all_clusters.csv")

pp4 <- two.class.test(obs.counts, tip.exp, cond.control="E105", cond.treatment="E135",to.plot=T)
pp4
write.csv(pp4,"UMAP_MERGED_E105vsE135-pvalues_all_clusters.csv")

pp5 <- two.class.test(obs.counts, tip.exp, cond.control="E115", cond.treatment="E135",to.plot=T)
pp5
write.csv(pp5,"UMAP_MERGED_E115vsE135-pvalues_all_clusters.csv")

#proportion of cells per cluster and sample
propotions_cells_percluster <- prop.table(table(all.merged.SCTCCSt.labelled$merged.clusters, all.merged.SCTCCSt.labelled$stage), margin = 2)
propotions_cells_percluster.df <- as.data.frame(propotions_cells_percluster)
t(propotions_cells_percluster)

write.csv(propotions_cells_percluster.df,"UMAP_MERGED_cell_proportions_all_clusters.csv")

ggplot(propotions_cells_percluster.df, aes(fill= Var2, color = Var2, x=Var1,y=Freq)) + 
  geom_bar(position="dodge", stat = "identity") + theme_bw() 

ggplot(propotions_cells_percluster.df, aes(fill= Var1, x=Var2,y=Freq)) + 
  geom_bar(position="fill", stat = "identity") + 
  geom_text(aes(label = sprintf("%0.2f", round(Freq, 2))), color="black", size=2,position = position_stack(vjust = 0.5))

#some plots  & basic proportion of cells per cluster and sample

pdf(file=paste0("all.merged.SCTCCSt.labelled_cell_proportion.pdf"))

DimPlot(all.merged.SCTCCSt.labelled, reduction = "umap",label = TRUE, label.size = 3)
cell.num <- table(Idents(all.merged.SCTCCSt.labelled))
clusterLabels <- paste(names(cell.num), paste0("(", cell.num, ")"))
clusterBreaks <- names(cell.num)
DimPlot(all.merged.SCTCCSt.labelled, reduction = "umap", label=TRUE) + scale_colour_discrete(breaks = clusterBreaks,labels = clusterLabels) 

propotions_cells_percluster <- prop.table(table(all.merged.SCTCCSt.labelled$merged.clusters, all.merged.SCTCCSt.labelled$stage), margin = 2)
propotions_cells_percluster.df <- as.data.frame(propotions_cells_percluster)
t(propotions_cells_percluster)

ggplot(propotions_cells_percluster.df, aes(fill= Var2, x=Var1,y=Freq)) + 
  geom_bar(position="dodge", stat = "identity") + theme_bw() + 
  geom_text(aes(label = sprintf("%0.2f", round(Freq, 2))), color="black", size=2, position = position_dodge(0.9), vjust = 0)

ggplot(propotions_cells_percluster.df, aes(fill= Var1, x=Var2,y=Freq)) + 
  geom_bar(position="fill", stat = "identity") +
  geom_text(aes(label = sprintf("%0.2f", round(Freq, 2))), color="black", size=2,position = position_stack(vjust = 0.5))

propotions_cells_percluster <- prop.table(table(all.merged.SCTCCSt.labelled$merged.clusters, all.merged.SCTCCSt.labelled$gene_class), margin = 2)
propotions_cells_percluster.df <- as.data.frame(propotions_cells_percluster)
t(propotions_cells_percluster)

ggplot(propotions_cells_percluster.df, aes(fill= Var2, x=Var1,y=Freq)) + 
  geom_bar(position="dodge", stat = "identity") + theme_bw() +
  geom_text(aes(label = sprintf("%0.2f", round(Freq, 2))), color="black", size=2.5, position = position_dodge(0.9), vjust = 0)

ggplot(propotions_cells_percluster.df, aes(fill= Var1, x=Var2,y=Freq)) + 
  geom_bar(position="fill", stat = "identity") + 
  geom_text(aes(label = sprintf("%0.2f", round(Freq, 2))), color="black", size=2,position = position_stack(vjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1))
dev.off()

## Cell proportion based on EYFP, Shox2, CRE expression classification

propotions_cells_percluster3 <- prop.table(table(all.merged.SCTCCSt.labelled$gene_class, all.merged.SCTCCSt.labelled$merged.clusters), margin = 2)
propotions_cells_percluster3.df <- as.data.frame(propotions_cells_percluster3)
t(propotions_cells_percluster3)

ggplot(propotions_cells_percluster3.df, aes(fill= Var1, x=Var2,y=Freq)) + 
  geom_bar(position="fill", stat = "identity") + 
  geom_text(aes(label = sprintf("%0.2f", round(Freq, 2))), color="black", size=2,position = position_stack(vjust = 0.5))

propotions_cells_percluster4 <- prop.table(table(all.merged.SCTCCSt.labelled$gene_class, all.merged.SCTCCSt.labelled$stage), margin = 2)
propotions_cells_percluster4.df <- as.data.frame(propotions_cells_percluster4)
t(propotions_cells_percluster4)

ggplot(propotions_cells_percluster4.df, aes(fill= Var1, x=Var2,y=Freq)) + 
  geom_bar(position="fill", stat = "identity") + 
  geom_text(aes(label = sprintf("%0.2f", round(Freq, 2))), color="black", size=2,position = position_stack(vjust = 0.5))

write.csv(propotions_cells_percluster4.df,"UMAP_MERGED_all_clusters_cell_proportions_perstage_of_EYFP_Shox2_CRE_expressing_cells.csv")

all.merged.SCTCCSt.labelled$geneclass.stage <- paste(all.merged.SCTCCSt.labelled$gene_class,all.merged.SCTCCSt.labelled$stage, sep = "_")
summary(all.merged.SCTCCSt.labelled@meta.data)

propotions_cells_percluster5 <- prop.table(table(all.merged.SCTCCSt.labelled$geneclass.stage, all.merged.SCTCCSt.labelled$merged.clusters), margin = 2)
propotions_cells_percluster5.df <- as.data.frame(propotions_cells_percluster5)
t(propotions_cells_percluster5)

propotions_cells_percluster5.df$stage <- propotions_cells_percluster5.df$Var1
propotions_cells_percluster5.df$stage <- as.character(gsub(".*(E105|E115|E125|E135).*", "\\1", propotions_cells_percluster5.df$stage))
head(propotions_cells_percluster5.df)

propotions_cells_percluster5.df$gene_class <- propotions_cells_percluster5.df$Var1
propotions_cells_percluster5.df$gene_class <- gsub("_E[0-9]+$","", propotions_cells_percluster5.df$gene_class)
head(propotions_cells_percluster5.df)

ggplot(propotions_cells_percluster5.df, aes(fill= gene_class, x=Var2,y=Freq)) + 
  geom_bar(position="fill", stat = "identity") + 
  facet_grid(~stage)


#### J. Some plots #####
all.merged.SCTCCSt.labelled <- readRDS("/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.merged.SCTCCSt.labelled.rds")
all.merged.SCTCCSt.labelled
summary(all.merged.SCTCCSt.labelled@meta.data)
DimPlot(all.merged.SCTCCSt.labelled, reduction = "umap",label = TRUE, label.size = 3)

DefaultAssay(all.merged.SCTCCSt.labelled) <- "RNA"
Idents(all.merged.SCTCCSt.labelled) <- "gene_class"
levels(all.merged.SCTCCSt.labelled)
levels(all.merged.SCTCCSt.labelled) <-c("EYFPneg_Shox2neg_CREneg_cells","EYFPneg_Shox2neg_CREpos_cells","EYFPneg_Shox2pos_CREneg_cells",
                                        "EYFPneg_Shox2pos_CREpos_cells","EYFPpos_Shox2neg_CREneg_cells","EYFPpos_Shox2neg_CREpos_cells",
                                        "EYFPpos_Shox2pos_CREneg_cells","EYFPpos_Shox2pos_CREpos_cells")

FeatureScatter(all.merged.SCTCCSt.labelled, feature1 = "EYFP", feature2 = "Shox2")
FeatureScatter(all.merged.SCTCCSt.labelled, feature1 = "CRE", feature2 = "Shox2")

FeatureScatter(all.merged.SCTCCSt.labelled,slot = "counts", feature1 = "EYFP", feature2 = "CRE")
FeatureScatter(all.merged.SCTCCSt.labelled,slot = "data", feature1 = "EYFP", feature2 = "CRE")

Idents(all.merged.SCTCCSt.labelled) <- "merged.clusters"

FeaturePlot(object = all.merged.SCTCCSt.labelled, features = "Shox2" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.merged.SCTCCSt.labelled, features = "Col1a1" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.merged.SCTCCSt.labelled, features = "EYFP" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.merged.SCTCCSt.labelled, features = "CRE" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.merged.SCTCCSt.labelled, features = "Rsrc1" ,cols = c("grey95","darkred"), combine = FALSE)

# Violin Plots
## We finally choose to use the Noise option = TRUE
DefaultAssay(all.merged.SCTCCSt.labelled) <- "RNA"
all.merged.SCTCCSt.labelled

VlnPlot(object = all.merged.SCTCCSt.labelled, assay = "RNA", features = "Shox2", group.by = "stage", 
        pt.size = 0, y.max = 4, add.noise = TRUE) + ggtitle("Shox2") 

ggsave("Vlnplot_Shox2_all_clusters_noiseT_split_by_stage.png", height = 5, width = 5)
ggsave("Vlnplot_Shox2_all_clusters_noiseT_split_by_stage.eps", height = 5, width = 5)

VlnPlot(object = all.merged.SCTCCSt.labelled, assay = "RNA", features = "CRE", group.by = "stage", 
        pt.size = 0, y.max = 4, add.noise = TRUE) + ggtitle("CRE") 

ggsave("Vlnplot_CRE_all_clusters_noiseT_split_by_stage.png", height = 5, width = 5)
ggsave("Vlnplot_CRE_all_clusters_noiseT_split_by_stage.eps", height = 5, width = 5)


VlnPlot(object = all.merged.SCTCCSt.labelled, assay = "RNA", features = "Shox2", group.by = "merged.clusters", 
        pt.size = 0, y.max = 4, add.noise = TRUE) + ggtitle("Shox2") 

ggsave("Vlnplot_Shox2_all_clusters_noiseT_split_by_allclusters.png", height = 5, width = 5)
ggsave("Vlnplot_Shox2_all_clusters_noiseT_split_by_allclusters.eps", height = 5, width = 5)

VlnPlot(object = all.merged.SCTCCSt.labelled, assay = "RNA", features = "CRE", group.by = "merged.clusters", 
        pt.size = 0, y.max = 4, add.noise = TRUE) + ggtitle("CRE") 

ggsave("Vlnplot_CRE_all_clusters_noiseT_split_by_allclusters.png", height = 5, width = 5)
ggsave("Vlnplot_CRE_all_clusters_noiseT_split_by_allclusters.eps", height = 5, width = 5)

#### J.1 Shox2 density plot per stage considering all the clusters ####
all.merged.SCTCCSt.labelled <- readRDS("/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.merged.SCTCCSt.labelled.rds")
all.merged.SCTCCSt.labelled
summary(all.merged.SCTCCSt.labelled@meta.data)
DimPlot(all.merged.SCTCCSt.labelled, reduction = "umap",label = TRUE, label.size = 3)

##per sample
Idents(all.merged.SCTCCSt.labelled) <- "stage"
levels(all.merged.SCTCCSt.labelled)

Shox2.expression <- all.merged.SCTCCSt.labelled@assays$RNA@data["Shox2",]
df1 <- data.frame(Shox2.expression)
df1$stage <- row.names(df1)
df1$stage <- gsub("_[A-Z]+$","", df1$stage)
df1$stage <- gsub("HL_","", df1$stage)
table(df1$stage)

#Using y=..density.. scales the histograms so the area under each is 1
ggplot(df1, aes(Shox2.expression,fill=stage, color=stage)) + 
  geom_histogram(aes(y=..density..),position = "dodge", alpha=0.5) + 
  geom_density(alpha=0.05,trim=TRUE) +
  scale_y_continuous(trans = 'sqrt') +
  ggtitle("Shox2 global")

ggplot(df1, aes(Shox2.expression,fill=stage, color=stage)) + 
  geom_density(alpha=0.05,trim=TRUE) +
  scale_y_continuous(trans = 'sqrt') +
  ggtitle("Shox2 global")

#### J.2 CRE density plot per stage considering all the clusters ####
all.merged.SCTCCSt.labelled <- readRDS("/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.merged.SCTCCSt.labelled.rds")
all.merged.SCTCCSt.labelled
summary(all.merged.SCTCCSt.labelled@meta.data)
DimPlot(all.merged.SCTCCSt.labelled, reduction = "umap",label = TRUE, label.size = 3)

##per sample
Idents(all.merged.SCTCCSt.labelled) <- "stage"
levels(all.merged.SCTCCSt.labelled)

CRE.expression <- all.merged.SCTCCSt.labelled@assays$RNA@data["CRE",]
df1 <- data.frame(CRE.expression)
df1$stage <- row.names(df1)
df1$stage <- gsub("_[A-Z]+$","", df1$stage)
df1$stage <- gsub("HL_","", df1$stage)
table(df1$stage)

#Using y=..density.. scales the histograms so the area under each is 1
ggplot(df1, aes(CRE.expression,fill=stage, color=stage)) + 
  geom_histogram(aes(y=..density..),position = "dodge", alpha=0.5) + 
  geom_density(alpha=0.05,trim=TRUE) +
  scale_y_continuous(trans = 'sqrt') +
  ggtitle("CRE global")

ggplot(df1, aes(CRE.expression,fill=stage, color=stage)) + 
  geom_density(alpha=0.05,trim=TRUE) +
  scale_y_continuous(trans = 'sqrt') +
  ggtitle("CRE global")

#### J.3 EYFP density plot per stage considering all the clusters ####
all.merged.SCTCCSt.labelled <- readRDS("/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.merged.SCTCCSt.labelled.rds")
all.merged.SCTCCSt.labelled
summary(all.merged.SCTCCSt.labelled@meta.data)
DimPlot(all.merged.SCTCCSt.labelled, reduction = "umap",label = TRUE, label.size = 3)

##per sample
Idents(all.merged.SCTCCSt.labelled) <- "stage"
levels(all.merged.SCTCCSt.labelled)

EYFP.expression <- all.merged.SCTCCSt.labelled@assays$RNA@data["EYFP",]
df1 <- data.frame(EYFP.expression)
df1$stage <- row.names(df1)
df1$stage <- gsub("_[A-Z]+$","", df1$stage)
df1$stage <- gsub("HL_","", df1$stage)
table(df1$stage)

#Using y=..density.. scales the histograms so the area under each is 1
ggplot(df1, aes(EYFP.expression,fill=stage, color=stage)) + 
  geom_histogram(aes(y=..density..),position = "dodge", alpha=0.5) + 
  geom_density(alpha=0.05,trim=TRUE) +
  scale_y_continuous(trans = 'sqrt') +
  ggtitle("EYFP global")

ggplot(df1, aes(EYFP.expression,fill=stage, color=stage)) + 
  geom_density(alpha=0.05,trim=TRUE) +
  scale_y_continuous(trans = 'sqrt') +
  ggtitle("EYFP global")

#### J.4 Shox2 density plot per gene class considering all the clusters ####
all.merged.SCTCCSt.labelled <- readRDS("/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.merged.SCTCCSt.labelled.rds")
all.merged.SCTCCSt.labelled
summary(all.merged.SCTCCSt.labelled@meta.data)
table(all.merged.SCTCCSt.labelled$gene_class)
DimPlot(all.merged.SCTCCSt.labelled, reduction = "umap",label = TRUE, label.size = 3)

##per sample
Idents(all.merged.SCTCCSt.labelled) <- "gene_class"
levels(all.merged.SCTCCSt.labelled)

#EYFPpos_Shox2pos_CREpos_cells
EYFPpos_Shox2pos_CREpos_cells <- subset(all.merged.SCTCCSt.labelled, idents = "EYFPpos_Shox2pos_CREpos_cells")

Shox2.EYFPpos_Shox2pos_CREpos_cells <- EYFPpos_Shox2pos_CREpos_cells@assays$RNA@data["Shox2",]
dfa <- data.frame(Shox2.EYFPpos_Shox2pos_CREpos_cells)
dfa$stage <- row.names(dfa)
dfa$stage <- gsub("_[A-Z]+$","", dfa$stage)
dfa$stage <- gsub("HL_","", dfa$stage)
table(dfa$stage)

#Using y=..density.. scales the histograms so the area under each is 1
ggplot(dfa, aes(Shox2.EYFPpos_Shox2pos_CREpos_cells,fill=stage, color=stage)) + 
  geom_density(alpha=0.05,trim=TRUE) +
  ggtitle("Shox2 EYFPpos_Shox2pos_CREpos_cells")

#EYFPpos_Shox2pos_CREneg_cells
EYFPpos_Shox2pos_CREneg_cells <- subset(all.merged.SCTCCSt.labelled, idents = "EYFPpos_Shox2pos_CREneg_cells")

Shox2.EYFPpos_Shox2pos_CREneg_cells <- EYFPpos_Shox2pos_CREneg_cells@assays$RNA@data["Shox2",]
dfa <- data.frame(Shox2.EYFPpos_Shox2pos_CREneg_cells)
dfa$stage <- row.names(dfa)
dfa$stage <- gsub("_[A-Z]+$","", dfa$stage)
dfa$stage <- gsub("HL_","", dfa$stage)
table(dfa$stage)

#Using y=..density.. scales the histograms so the area under each is 1
ggplot(dfa, aes(Shox2.EYFPpos_Shox2pos_CREneg_cells,fill=stage, color=stage)) + 
  geom_density(alpha=0.05,trim=TRUE) +
  ggtitle("Shox2 EYFPpos_Shox2pos_CREneg_cells")

#EYFPpos_Shox2neg_CREpos_cells
EYFPpos_Shox2neg_CREpos_cells <- subset(all.merged.SCTCCSt.labelled, idents = "EYFPpos_Shox2neg_CREpos_cells")

Shox2.EYFPpos_Shox2neg_CREpos_cells <- EYFPpos_Shox2neg_CREpos_cells@assays$RNA@data["Shox2",]
dfa <- data.frame(Shox2.EYFPpos_Shox2neg_CREpos_cells)
dfa$stage <- row.names(dfa)
dfa$stage <- gsub("_[A-Z]+$","", dfa$stage)
dfa$stage <- gsub("HL_","", dfa$stage)
table(dfa$stage)

#Using y=..density.. scales the histograms so the area under each is 1
ggplot(dfa, aes(Shox2.EYFPpos_Shox2neg_CREpos_cells,fill=stage, color=stage)) + 
  geom_density(alpha=0.05,trim=TRUE) +
  ggtitle("Shox2 EYFPpos_Shox2neg_CREpos_cells")

#EYFPpos_Shox2neg_CREneg_cells
EYFPpos_Shox2neg_CREneg_cells <- subset(all.merged.SCTCCSt.labelled, idents = "EYFPpos_Shox2neg_CREneg_cells")

Shox2.EYFPpos_Shox2neg_CREneg_cells <- EYFPpos_Shox2neg_CREneg_cells@assays$RNA@data["Shox2",]
dfa <- data.frame(Shox2.EYFPpos_Shox2neg_CREneg_cells)
dfa$stage <- row.names(dfa)
dfa$stage <- gsub("_[A-Z]+$","", dfa$stage)
dfa$stage <- gsub("HL_","", dfa$stage)
table(dfa$stage)

#Using y=..density.. scales the histograms so the area under each is 1
ggplot(dfa, aes(Shox2.EYFPpos_Shox2neg_CREneg_cells,fill=stage, color=stage)) + 
  geom_density(alpha=0.05,trim=TRUE) +
  ggtitle("Shox2 EYFPpos_Shox2neg_CREneg_cells")

#EYFPneg_Shox2pos_CREpos_cells
EYFPneg_Shox2pos_CREpos_cells <- subset(all.merged.SCTCCSt.labelled, idents = "EYFPneg_Shox2pos_CREpos_cells")

Shox2.EYFPneg_Shox2pos_CREpos_cells <- EYFPneg_Shox2pos_CREpos_cells@assays$RNA@data["Shox2",]
dfa <- data.frame(Shox2.EYFPneg_Shox2pos_CREpos_cells)
dfa$stage <- row.names(dfa)
dfa$stage <- gsub("_[A-Z]+$","", dfa$stage)
dfa$stage <- gsub("HL_","", dfa$stage)
table(dfa$stage)

#Using y=..density.. scales the histograms so the area under each is 1
ggplot(dfa, aes(Shox2.EYFPneg_Shox2pos_CREpos_cells,fill=stage, color=stage)) + 
  geom_density(alpha=0.05,trim=TRUE) +
  ggtitle("Shox2 EYFPneg_Shox2pos_CREpos_cells")

#EYFPneg_Shox2neg_CREpos_cells
EYFPneg_Shox2neg_CREpos_cells <- subset(all.merged.SCTCCSt.labelled, idents = "EYFPneg_Shox2neg_CREpos_cells")

Shox2.EYFPneg_Shox2neg_CREpos_cells <- EYFPneg_Shox2neg_CREpos_cells@assays$RNA@data["Shox2",]
dfa <- data.frame(Shox2.EYFPneg_Shox2neg_CREpos_cells)
dfa$stage <- row.names(dfa)
dfa$stage <- gsub("_[A-Z]+$","", dfa$stage)
dfa$stage <- gsub("HL_","", dfa$stage)
table(dfa$stage)

#Using y=..density.. scales the histograms so the area under each is 1
ggplot(dfa, aes(Shox2.EYFPneg_Shox2neg_CREpos_cells,fill=stage, color=stage)) + 
  geom_density(alpha=0.05,trim=TRUE) +
  ggtitle("Shox2 EYFPneg_Shox2neg_CREpos_cells")

#EYFPneg_Shox2neg_CREneg_cells
EYFPneg_Shox2neg_CREneg_cells <- subset(all.merged.SCTCCSt.labelled, idents = "EYFPneg_Shox2neg_CREneg_cells")

Shox2.EYFPneg_Shox2neg_CREneg_cells <- EYFPneg_Shox2neg_CREneg_cells@assays$RNA@data["Shox2",]
dfa <- data.frame(Shox2.EYFPneg_Shox2neg_CREneg_cells)
dfa$stage <- row.names(dfa)
dfa$stage <- gsub("_[A-Z]+$","", dfa$stage)
dfa$stage <- gsub("HL_","", dfa$stage)
table(dfa$stage)

#Using y=..density.. scales the histograms so the area under each is 1
ggplot(dfa, aes(Shox2.EYFPneg_Shox2neg_CREneg_cells,fill=stage, color=stage)) + 
  geom_density(alpha=0.05,trim=TRUE) +
  ggtitle("Shox2 EYFPneg_Shox2neg_CREneg_cells")

#EYFPneg_Shox2pos_CREneg_cells
EYFPneg_Shox2pos_CREneg_cells <- subset(all.merged.SCTCCSt.labelled, idents = "EYFPneg_Shox2pos_CREneg_cells")

Shox2.EYFPneg_Shox2pos_CREneg_cells <- EYFPneg_Shox2pos_CREneg_cells@assays$RNA@data["Shox2",]
dfa <- data.frame(Shox2.EYFPneg_Shox2pos_CREneg_cells)
dfa$stage <- row.names(dfa)
dfa$stage <- gsub("_[A-Z]+$","", dfa$stage)
dfa$stage <- gsub("HL_","", dfa$stage)
table(dfa$stage)

#Using y=..density.. scales the histograms so the area under each is 1
ggplot(dfa, aes(Shox2.EYFPneg_Shox2pos_CREneg_cells,fill=stage, color=stage)) + 
  geom_density(alpha=0.05,trim=TRUE) +
  ggtitle("Shox2 EYFPneg_Shox2pos_CREneg_cells")

#### J.5 EYFP density plot per gene class considering all the clusters ####
all.merged.SCTCCSt.labelled <- readRDS("/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.merged.SCTCCSt.labelled.rds")
all.merged.SCTCCSt.labelled
summary(all.merged.SCTCCSt.labelled@meta.data)
table(all.merged.SCTCCSt.labelled$gene_class)
DimPlot(all.merged.SCTCCSt.labelled, reduction = "umap",label = TRUE, label.size = 3)

##per sample
Idents(all.merged.SCTCCSt.labelled) <- "gene_class"
levels(all.merged.SCTCCSt.labelled)

#EYFPpos_Shox2pos_CREpos_cells
EYFPpos_Shox2pos_CREpos_cells <- subset(all.merged.SCTCCSt.labelled, idents = "EYFPpos_Shox2pos_CREpos_cells")

Shox2.EYFPpos_Shox2pos_CREpos_cells <- EYFPpos_Shox2pos_CREpos_cells@assays$RNA@data["EYFP",]
dfa <- data.frame(Shox2.EYFPpos_Shox2pos_CREpos_cells)
dfa$stage <- row.names(dfa)
dfa$stage <- gsub("_[A-Z]+$","", dfa$stage)
dfa$stage <- gsub("HL_","", dfa$stage)
table(dfa$stage)

#Using y=..density.. scales the histograms so the area under each is 1
ggplot(dfa, aes(Shox2.EYFPpos_Shox2pos_CREpos_cells,fill=stage, color=stage)) + 
  geom_density(alpha=0.05,trim=TRUE) +
  ggtitle("EYFP EYFPpos_Shox2pos_CREpos_cells")

#EYFPpos_Shox2pos_CREneg_cells
EYFPpos_Shox2pos_CREneg_cells <- subset(all.merged.SCTCCSt.labelled, idents = "EYFPpos_Shox2pos_CREneg_cells")

Shox2.EYFPpos_Shox2pos_CREneg_cells <- EYFPpos_Shox2pos_CREneg_cells@assays$RNA@data["EYFP",]
dfa <- data.frame(Shox2.EYFPpos_Shox2pos_CREneg_cells)
dfa$stage <- row.names(dfa)
dfa$stage <- gsub("_[A-Z]+$","", dfa$stage)
dfa$stage <- gsub("HL_","", dfa$stage)
table(dfa$stage)

#Using y=..density.. scales the histograms so the area under each is 1
ggplot(dfa, aes(Shox2.EYFPpos_Shox2pos_CREneg_cells,fill=stage, color=stage)) + 
  geom_density(alpha=0.05,trim=TRUE) +
  ggtitle("EYFP EYFPpos_Shox2pos_CREneg_cells")

#EYFPpos_Shox2neg_CREpos_cells
EYFPpos_Shox2neg_CREpos_cells <- subset(all.merged.SCTCCSt.labelled, idents = "EYFPpos_Shox2neg_CREpos_cells")

Shox2.EYFPpos_Shox2neg_CREpos_cells <- EYFPpos_Shox2neg_CREpos_cells@assays$RNA@data["EYFP",]
dfa <- data.frame(Shox2.EYFPpos_Shox2neg_CREpos_cells)
dfa$stage <- row.names(dfa)
dfa$stage <- gsub("_[A-Z]+$","", dfa$stage)
dfa$stage <- gsub("HL_","", dfa$stage)
table(dfa$stage)

#Using y=..density.. scales the histograms so the area under each is 1
ggplot(dfa, aes(Shox2.EYFPpos_Shox2neg_CREpos_cells,fill=stage, color=stage)) + 
  geom_density(alpha=0.05,trim=TRUE) +
  ggtitle("EYFP EYFPpos_Shox2neg_CREpos_cells")

#EYFPpos_Shox2neg_CREneg_cells
EYFPpos_Shox2neg_CREneg_cells <- subset(all.merged.SCTCCSt.labelled, idents = "EYFPpos_Shox2neg_CREneg_cells")

Shox2.EYFPpos_Shox2neg_CREneg_cells <- EYFPpos_Shox2neg_CREneg_cells@assays$RNA@data["EYFP",]
dfa <- data.frame(Shox2.EYFPpos_Shox2neg_CREneg_cells)
dfa$stage <- row.names(dfa)
dfa$stage <- gsub("_[A-Z]+$","", dfa$stage)
dfa$stage <- gsub("HL_","", dfa$stage)
table(dfa$stage)

#Using y=..density.. scales the histograms so the area under each is 1
ggplot(dfa, aes(Shox2.EYFPpos_Shox2neg_CREneg_cells,fill=stage, color=stage)) + 
  geom_density(alpha=0.05,trim=TRUE) +
  ggtitle("EYFP EYFPpos_Shox2neg_CREneg_cells")

#EYFPneg_Shox2pos_CREpos_cells
EYFPneg_Shox2pos_CREpos_cells <- subset(all.merged.SCTCCSt.labelled, idents = "EYFPneg_Shox2pos_CREpos_cells")

Shox2.EYFPneg_Shox2pos_CREpos_cells <- EYFPneg_Shox2pos_CREpos_cells@assays$RNA@data["EYFP",]
dfa <- data.frame(Shox2.EYFPneg_Shox2pos_CREpos_cells)
dfa$stage <- row.names(dfa)
dfa$stage <- gsub("_[A-Z]+$","", dfa$stage)
dfa$stage <- gsub("HL_","", dfa$stage)
table(dfa$stage)

#Using y=..density.. scales the histograms so the area under each is 1
ggplot(dfa, aes(Shox2.EYFPneg_Shox2pos_CREpos_cells,fill=stage, color=stage)) + 
  geom_density(alpha=0.05,trim=TRUE) +
  ggtitle("EYFP EYFPneg_Shox2pos_CREpos_cells")

#EYFPneg_Shox2neg_CREpos_cells
EYFPneg_Shox2neg_CREpos_cells <- subset(all.merged.SCTCCSt.labelled, idents = "EYFPneg_Shox2neg_CREpos_cells")

Shox2.EYFPneg_Shox2neg_CREpos_cells <- EYFPneg_Shox2neg_CREpos_cells@assays$RNA@data["EYFP",]
dfa <- data.frame(Shox2.EYFPneg_Shox2neg_CREpos_cells)
dfa$stage <- row.names(dfa)
dfa$stage <- gsub("_[A-Z]+$","", dfa$stage)
dfa$stage <- gsub("HL_","", dfa$stage)
table(dfa$stage)

#Using y=..density.. scales the histograms so the area under each is 1
ggplot(dfa, aes(Shox2.EYFPneg_Shox2neg_CREpos_cells,fill=stage, color=stage)) + 
  geom_density(alpha=0.05,trim=TRUE) +
  ggtitle("EYFP EYFPneg_Shox2neg_CREpos_cells")

#EYFPneg_Shox2neg_CREneg_cells
EYFPneg_Shox2neg_CREneg_cells <- subset(all.merged.SCTCCSt.labelled, idents = "EYFPneg_Shox2neg_CREneg_cells")

Shox2.EYFPneg_Shox2neg_CREneg_cells <- EYFPneg_Shox2neg_CREneg_cells@assays$RNA@data["EYFP",]
dfa <- data.frame(Shox2.EYFPneg_Shox2neg_CREneg_cells)
dfa$stage <- row.names(dfa)
dfa$stage <- gsub("_[A-Z]+$","", dfa$stage)
dfa$stage <- gsub("HL_","", dfa$stage)
table(dfa$stage)

#Using y=..density.. scales the histograms so the area under each is 1
ggplot(dfa, aes(Shox2.EYFPneg_Shox2neg_CREneg_cells,fill=stage, color=stage)) + 
  geom_density(alpha=0.05,trim=TRUE) +
  ggtitle("EYFP EYFPneg_Shox2neg_CREneg_cells")

#EYFPneg_Shox2pos_CREneg_cells
EYFPneg_Shox2pos_CREneg_cells <- subset(all.merged.SCTCCSt.labelled, idents = "EYFPneg_Shox2pos_CREneg_cells")

Shox2.EYFPneg_Shox2pos_CREneg_cells <- EYFPneg_Shox2pos_CREneg_cells@assays$RNA@data["EYFP",]
dfa <- data.frame(Shox2.EYFPneg_Shox2pos_CREneg_cells)
dfa$stage <- row.names(dfa)
dfa$stage <- gsub("_[A-Z]+$","", dfa$stage)
dfa$stage <- gsub("HL_","", dfa$stage)
table(dfa$stage)

#Using y=..density.. scales the histograms so the area under each is 1
ggplot(dfa, aes(Shox2.EYFPneg_Shox2pos_CREneg_cells,fill=stage, color=stage)) + 
  geom_density(alpha=0.05,trim=TRUE) +
  ggtitle("EYFP EYFPneg_Shox2pos_CREneg_cells")

#### K.1 Subsetting mesenchyme cluster ####
all.merged.SCTCCSt.labelled <- readRDS("/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.merged.SCTCCSt.labelled.rds")
all.merged.SCTCCSt.labelled
summary(all.merged.SCTCCSt.labelled@meta.data)
DimPlot(all.merged.SCTCCSt.labelled, reduction = "umap",label = TRUE, label.size = 3)

all.mesenchyme <- subset(all.merged.SCTCCSt.labelled, idents = c("Mesenchyme"))
all.mesenchyme
table(all.mesenchyme$stage)
table(all.mesenchyme$merged.clusters)
DimPlot(all.mesenchyme, reduction = "umap", label = TRUE)

# save as seurat object 
saveRDS(all.mesenchyme, file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.mesenchyme.rds")

#### K.2 UMAP Mesenchyme cells ####
all.mesenchyme <- readRDS("/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.mesenchyme.rds")
all.mesenchyme
summary(all.mesenchyme@meta.data)

DimPlot(all.mesenchyme, reduction = "umap", label = TRUE)
Idents(all.mesenchyme) <- "stage"
Idents(all.mesenchyme)
DimPlot(all.mesenchyme, reduction = "umap")

#all.mesenchyme <- RunUMAP(object = all.mesenchyme, dims = 1:50)

#all.mesenchyme <- RunUMAP(object = all.mesenchyme, dims = 1:50,n.neighbors = 50L, metric = "euclidean",min.dist = 0.5)

#all.mesenchyme <- RunUMAP(object = all.mesenchyme, dims = 1:50,n.neighbors = 50L, metric = "cosine",min.dist = 0.5)

#all.mesenchyme <- RunUMAP(object = all.mesenchyme, dims = 1:50,n.neighbors = 50L, metric = "cosine",min.dist = 1)

#all.mesenchyme <- RunUMAP(object = all.mesenchyme, dims = 1:10,n.neighbors = 30L, metric = "cosine",min.dist = 0.1,spread = 2)

#all.mesenchyme <- RunUMAP(all.mesenchyme,dims = c(1:10), n.neighbors = 40L, min.dist = 0.5, metric = "euclidean", spread = 1)

#all.mesenchyme <- RunUMAP(all.mesenchyme,dims = c(1:10), n.neighbors = 40L, min.dist = 0.2, metric = "euclidean", spread = 1)

#all.mesenchyme <- RunUMAP(all.mesenchyme,dims = c(1:10), n.neighbors = 30L, min.dist = 0.2, metric = "euclidean", spread = 1)

all.mesenchyme <- RunUMAP(all.mesenchyme,dims = c(1:10), n.neighbors = 30L, min.dist = 0.5, metric = "euclidean", spread = 1)

DimPlot(all.mesenchyme, reduction = "umap", label = FALSE)

DefaultAssay(all.mesenchyme) <- "RNA"

FeaturePlot(object = all.mesenchyme, features = "Hoxd13" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Shox2" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "EYFP" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "CRE" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Sox9" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Ihh" ,cols = c("grey95","darkred"), combine = FALSE)

DefaultAssay(all.mesenchyme) <- "SCT"
all.mesenchyme <- FindNeighbors(object = all.mesenchyme)

# save as seurat object 
saveRDS(all.mesenchyme, file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.UMAP_mesenchyme")

#### K.3 Trial Re-clustering Mesenchyme ####
all.mesenchyme <- readRDS(file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.UMAP_mesenchyme")
all.mesenchyme
summary(all.mesenchyme@meta.data)
DimPlot(all.mesenchyme, reduction = "umap", label = TRUE)

all.mesenchyme <- FindClusters(object = all.mesenchyme, resolution = seq(from = 0.1, to = 2, by = 0.1))
clustree(x = all.mesenchyme, prefix = "SCT_snn_res.")
summary(all.mesenchyme@meta.data)
pdf(file=paste0("all.mesenchyme_UMAP_clustering_options.pdf"))
DimPlot(object = all.mesenchyme, reduction = "umap", group.by = "SCT_snn_res.0.1", label = TRUE, label.size = 6)
DimPlot(object = all.mesenchyme, reduction = "umap", group.by = "SCT_snn_res.0.2", label = TRUE, label.size = 6) 
DimPlot(object = all.mesenchyme, reduction = "umap", group.by = "SCT_snn_res.0.3", label = TRUE, label.size = 6)
DimPlot(object = all.mesenchyme, reduction = "umap", group.by = "SCT_snn_res.0.4", label = TRUE, label.size = 6)
DimPlot(object = all.mesenchyme, reduction = "umap", group.by = "SCT_snn_res.0.5", label = TRUE, label.size = 6)
DimPlot(object = all.mesenchyme, reduction = "umap", group.by = "SCT_snn_res.0.6", label = TRUE, label.size = 6)
DimPlot(object = all.mesenchyme, reduction = "umap", group.by = "SCT_snn_res.0.7", label = TRUE, label.size = 6)
DimPlot(object = all.mesenchyme, reduction = "umap", group.by = "SCT_snn_res.0.8", label = TRUE, label.size = 6)
DimPlot(object = all.mesenchyme, reduction = "umap", group.by = "SCT_snn_res.0.9", label = TRUE, label.size = 6)
DimPlot(object = all.mesenchyme, reduction = "umap", group.by = "SCT_snn_res.1", label = TRUE, label.size = 6)
DimPlot(object = all.mesenchyme, reduction = "umap", group.by = "SCT_snn_res.1.1", label = TRUE, label.size = 6)
dev.off()

pdf(file=paste0("all.mesenchyme_Features.pdf"))
DefaultAssay(all.mesenchyme) <- "RNA"
FeaturePlot(object = all.mesenchyme, features = "Shox2" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Hoxd13" ,cols = c("grey95","darkred"), combine = FALSE)

FeaturePlot(object = all.mesenchyme, features = "Irx1" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Irx2" ,cols = c("grey95","darkred"), combine = FALSE)

FeaturePlot(object = all.mesenchyme, features = "Irx3" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Irx5" ,cols = c("grey95","darkred"), combine = FALSE)

FeaturePlot(object = all.mesenchyme, features = "Osr1" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Lum" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Kera" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Scx" ,cols = c("grey95","darkred"), combine = FALSE)

FeaturePlot(object = all.mesenchyme, features = "Shh" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Msx1" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Top2a" ,cols = c("grey95","darkred"), combine = FALSE)

FeaturePlot(object = all.mesenchyme, features = "Tcf4" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Tbx4" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Osr2" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Cxcl12" ,cols = c("grey95","darkred"), combine = FALSE)

FeaturePlot(object = all.mesenchyme, features = "Gdf5" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Ihh" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Runx2" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Sox9" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Acan" ,cols = c("grey95","darkred"), combine = FALSE)

FeaturePlot(object = all.mesenchyme, features = "Grem1" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Hand2" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Meis2" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Meis1" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Isl1" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Dpt" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Hoxa11" ,cols = c("grey95","darkred"), combine = FALSE)

FeaturePlot(object = all.mesenchyme, features = "EYFP" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "CRE" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "dsmCherry" ,cols = c("grey95","darkred"), order = TRUE, combine = FALSE)
dev.off()

#### K.4 Final Re-clustering Mesenchyme ####
all.mesenchyme <- readRDS(file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.UMAP_mesenchyme.rds")
all.mesenchyme
summary(all.mesenchyme@meta.data)
DimPlot(all.mesenchyme, reduction = "umap", label = TRUE)

all.mesenchyme <- FindClusters(object = all.mesenchyme, resolution = 1.1)
DimPlot(all.mesenchyme, reduction = "umap", label = TRUE)

DefaultAssay(all.mesenchyme) <- "RNA"
all.mesenchyme

subcluster0.markers = FindMarkers(all.mesenchyme, assay = "RNA", ident.1 = 0, grouping.var = "stage", 
                                  only.pos =TRUE, min.diff.pct=0.1, logfc.threshold = 0.3)
head (subcluster0.markers,15)
write.csv(subcluster0.markers,file=("subcluster0.markers.csv"))

subcluster1.markers = FindMarkers(all.mesenchyme, assay = "RNA", ident.1 = 1, grouping.var = "stage",
                                  only.pos =TRUE,min.diff.pct=0.1, logfc.threshold = 0.3)
head (subcluster1.markers,15)
write.csv(subcluster1.markers,file=("subcluster1.markers.csv"))

subcluster2.markers = FindMarkers(all.mesenchyme, assay = "RNA", ident.1 = 2, grouping.var = "stage",
                                  only.pos =TRUE, min.diff.pct=0.1, logfc.threshold = 0.3)
head (subcluster2.markers,15)
write.csv(subcluster2.markers,file=("subcluster2.markers.csv"))

subcluster3.markers = FindMarkers(all.mesenchyme, assay = "RNA", ident.1 = 3, grouping.var = "stage",
                                  only.pos =TRUE, min.diff.pct=0.1, logfc.threshold = 0.3)
head (subcluster3.markers,15)
write.csv(subcluster3.markers,file=("subcluster3.markers.csv"))

subcluster4.markers = FindMarkers(all.mesenchyme, assay = "RNA", ident.1 = 4, grouping.var = "stage",
                                  only.pos =TRUE, min.diff.pct=0.1, logfc.threshold = 0.3)
head (subcluster4.markers,15)
write.csv(subcluster4.markers,file=("subcluster4.markers.csv"))

subcluster5.markers = FindMarkers(all.mesenchyme, assay = "RNA", ident.1 = 5, grouping.var = "stage",
                                  only.pos =TRUE, min.diff.pct=0.1, logfc.threshold = 0.3)
head (subcluster5.markers,15)
write.csv(subcluster5.markers,file=("subcluster5.markers.csv"))

subcluster6.markers = FindMarkers(all.mesenchyme, assay = "RNA", ident.1 = 6, grouping.var = "stage",
                                  only.pos =TRUE, min.diff.pct=0.1, logfc.threshold = 0.3)
head (subcluster6.markers,15)
write.csv(subcluster6.markers,file=("subcluster6.markers.csv"))

subcluster7.markers = FindMarkers(all.mesenchyme, assay = "RNA", ident.1 = 7, grouping.var = "stage",
                                  only.pos =TRUE, min.diff.pct=0.1, logfc.threshold = 0.3)
head (subcluster7.markers,15)
write.csv(subcluster7.markers,file=("subcluster7.markers.csv"))

subcluster8.markers = FindMarkers(all.mesenchyme, assay = "RNA", ident.1 = 8, grouping.var = "stage",
                                  only.pos =TRUE, min.diff.pct=0.1, logfc.threshold = 0.3)
head (subcluster8.markers,15)
write.csv(subcluster8.markers,file=("subcluster8.markers.csv"))

subcluster9.markers = FindMarkers(all.mesenchyme, assay = "RNA", ident.1 = 9, grouping.var = "stage",
                                  only.pos =TRUE, min.diff.pct=0.1, logfc.threshold = 0.3)
head (subcluster9.markers,15)
write.csv(subcluster9.markers,file=("subcluster9.markers.csv"))

subcluster10.markers = FindMarkers(all.mesenchyme, assay = "RNA", ident.1 = 10, grouping.var = "stage",
                                   only.pos =TRUE, min.diff.pct=0.1, logfc.threshold = 0.3)
head (subcluster10.markers,15)
write.csv(subcluster10.markers,file=("subcluster10.markers.csv"))

subcluster11.markers = FindMarkers(all.mesenchyme, assay = "RNA", ident.1 = 11, grouping.var = "stage",
                                   only.pos =TRUE, min.diff.pct=0.1, logfc.threshold = 0.3)
head (subcluster11.markers,15)
write.csv(subcluster11.markers,file=("subcluster11.markers.csv"))

subcluster12.markers = FindMarkers(all.mesenchyme, assay = "RNA", ident.1 = 12, grouping.var = "stage",
                                   only.pos =TRUE, min.diff.pct=0.1, logfc.threshold = 0.3)
head (subcluster12.markers,15)
write.csv(subcluster12.markers,file=("subcluster12.markers.csv"))

subcluster13.markers = FindMarkers(all.mesenchyme, assay = "RNA", ident.1 = 13, grouping.var = "stage",
                                   only.pos =TRUE, min.diff.pct=0.1, logfc.threshold = 0.3)
head (subcluster13.markers,15)
write.csv(subcluster13.markers,file=("subcluster13.markers.csv"))

subcluster14.markers = FindMarkers(all.mesenchyme, assay = "RNA", ident.1 = 14, grouping.var = "stage",
                                   only.pos =TRUE, min.diff.pct=0.1, logfc.threshold = 0.3)
head (subcluster14.markers,15)
write.csv(subcluster14.markers,file=("subcluster14.markers.csv"))

subcluster15.markers = FindMarkers(all.mesenchyme, assay = "RNA", ident.1 = 15, grouping.var = "stage",
                                   only.pos =TRUE, min.diff.pct=0.1, logfc.threshold = 0.3)
head (subcluster15.markers,15)
write.csv(subcluster15.markers,file=("subcluster15.markers.csv"))

subcluster16.markers = FindMarkers(all.mesenchyme, assay = "RNA", ident.1 = 16, grouping.var = "stage",
                                   only.pos =TRUE, min.diff.pct=0.1, logfc.threshold = 0.3)
head (subcluster16.markers,15)
write.csv(subcluster16.markers,file=("subcluster16.markers.csv"))

subcluster17.markers = FindMarkers(all.mesenchyme, assay = "RNA", ident.1 = 17, grouping.var = "stage",
                                   only.pos =TRUE, min.diff.pct=0.1, logfc.threshold = 0.3)
head (subcluster17.markers,15)
write.csv(subcluster17.markers,file=("subcluster17.markers.csv"))


decluster.markers11_06 = FindMarkers(all.mesenchyme, ident.1 =11, ident.2 =c(0,6), only.pos = TRUE,
                                    min.diff.pct=0.2, logfc.threshold = 0.3)
head (decluster.markers11_06,15)

# save as seurat object 
saveRDS(all.mesenchyme, file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.mesenchyme_UMAP1-1.rds")

#### L.1 Some plots UMAP1.1 ####
all.mesenchyme <- readRDS(file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.mesenchyme_UMAP1-1.rds")
all.mesenchyme
summary(all.mesenchyme@meta.data)
DimPlot(all.mesenchyme, reduction = "umap", label = TRUE)

##Feature plots
DefaultAssay(all.mesenchyme) <- "RNA"

FeaturePlot(object = all.mesenchyme, features = "Fgf2" ,cols = c("grey95","darkred"), combine = FALSE)

## Vlnplots ##
DefaultAssay(all.mesenchyme) <- "RNA"
summary(all.mesenchyme@meta.data)
Idents(all.mesenchyme) <- "SCT_snn_res.1.1"

fun1 = function (x) {log1p(mean(x = expm1(x = x)))}

pdf(file=paste0("all.mesenchyme_vlnplots.pdf"))
VlnPlot(object = all.mesenchyme, assay = "RNA", features = "Shox2", group.by = "SCT_snn_res.1.1", pt.size = 0, y.max = 4) + ggtitle("Shox2") +
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4,colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(object = all.mesenchyme, assay = "RNA", features = "EYFP", group.by = "SCT_snn_res.1.1", pt.size = 0, y.max = 4) + ggtitle("EYFP") +
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4,colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(object = all.mesenchyme, assay = "RNA", features = "CRE", group.by = "SCT_snn_res.1.1", pt.size = 0, y.max = 4) + ggtitle("CRE") +
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4,colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(object = all.mesenchyme, assay = "RNA", features = "Hoxd13", group.by = "SCT_snn_res.1.1", pt.size = 0, y.max = 4) + ggtitle("Hoxd13") +
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4,colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 


Idents(all.mesenchyme) <- "stage"
VlnPlot(all.mesenchyme,assay = "RNA", features = "Shox2", group.by = "SCT_snn_res.1.1", idents = "E105", pt.size = 0, y.max = 4) + ggtitle("Shox2 E10.5")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "Shox2", group.by = "SCT_snn_res.1.1", idents = "E115", pt.size = 0, y.max = 4) + ggtitle("Shox2 E11.5")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "Shox2", group.by = "SCT_snn_res.1.1", idents = "E125", pt.size = 0, y.max = 4) + ggtitle("Shox2 E12.5")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "Shox2", group.by = "SCT_snn_res.1.1", idents = "E135", pt.size = 0, y.max = 4) + ggtitle("Shox2 E13.5")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "EYFP", group.by = "SCT_snn_res.1.1", idents = "E105", pt.size = 0, y.max = 4) + ggtitle("EYFP E10.5")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "EYFP", group.by = "SCT_snn_res.1.1", idents = "E115", pt.size = 0, y.max = 4) + ggtitle("EYFP E11.5")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "EYFP", group.by = "SCT_snn_res.1.1", idents = "E125", pt.size = 0, y.max = 4) + ggtitle("EYFP E12.5")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "EYFP", group.by = "SCT_snn_res.1.1", idents = "E135", pt.size = 0, y.max = 4) + ggtitle("EYFP E13.5")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "CRE", group.by = "SCT_snn_res.1.1", idents = "E105", pt.size = 0, y.max = 4) + ggtitle("CRE E10.5")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "CRE", group.by = "SCT_snn_res.1.1", idents = "E115", pt.size = 0, y.max = 4) + ggtitle("CRE E11.5")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "CRE", group.by = "SCT_snn_res.1.1", idents = "E125", pt.size = 0, y.max = 4) + ggtitle("CRE E12.5")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "CRE", group.by = "SCT_snn_res.1.1", idents = "E135", pt.size = 0, y.max = 4) + ggtitle("CRE E13.5")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "Hoxd13", group.by = "SCT_snn_res.1.1", idents = "E105", pt.size = 0, y.max = 4) + ggtitle("Hoxd13 E10.5")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "Hoxd13", group.by = "SCT_snn_res.1.1", idents = "E115", pt.size = 0, y.max = 4) + ggtitle("Hoxd13 E11.5")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "Hoxd13", group.by = "SCT_snn_res.1.1", idents = "E125", pt.size = 0, y.max = 4) + ggtitle("Hoxd13 E12.5")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "Hoxd13", group.by = "SCT_snn_res.1.1", idents = "E135", pt.size = 0, y.max = 4) + ggtitle("Hoxd13 E13.5")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 


Idents(all.mesenchyme) <- "gene_class"
VlnPlot(all.mesenchyme,assay = "RNA", features = "Shox2", group.by = "SCT_snn_res.1.1", idents = "EYFPneg_Shox2neg_CREneg_cells", pt.size = 0, y.max = 4) + ggtitle("Shox2 EYFPneg_Shox2neg_CREneg_cells")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "Shox2", group.by = "SCT_snn_res.1.1", idents = "EYFPneg_Shox2neg_CREpos_cells", pt.size = 0, y.max = 4) + ggtitle("Shox2 EYFPneg_Shox2neg_CREpos_cells")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "Shox2", group.by = "SCT_snn_res.1.1", idents = "EYFPneg_Shox2pos_CREneg_cells", pt.size = 0, y.max = 4) + ggtitle("Shox2 EYFPneg_Shox2pos_CREneg_cells")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "Shox2", group.by = "SCT_snn_res.1.1", idents = "EYFPneg_Shox2pos_CREpos_cells", pt.size = 0, y.max = 4) + ggtitle("Shox2 EYFPneg_Shox2pos_CREpos_cells")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "Shox2", group.by = "SCT_snn_res.1.1", idents = "EYFPpos_Shox2neg_CREneg_cells", pt.size = 0, y.max = 4) + ggtitle("Shox2 EYFPpos_Shox2neg_CREneg_cells")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "Shox2", group.by = "SCT_snn_res.1.1", idents = "EYFPpos_Shox2neg_CREpos_cells", pt.size = 0, y.max = 4) + ggtitle("Shox2 EYFPpos_Shox2neg_CREpos_cells")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "Shox2", group.by = "SCT_snn_res.1.1", idents = "EYFPpos_Shox2pos_CREneg_cells", pt.size = 0, y.max = 4) + ggtitle("Shox2 EYFPpos_Shox2pos_CREneg_cells")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "Shox2", group.by = "SCT_snn_res.1.1", idents = "EYFPpos_Shox2pos_CREpos_cells", pt.size = 0, y.max = 4) + ggtitle("Shox2 EYFPpos_Shox2pos_CREpos_cells")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "EYFP", group.by = "SCT_snn_res.1.1", idents = "EYFPneg_EYFPneg_CREneg_cells", pt.size = 0, y.max = 4) + ggtitle("EYFP EYFPneg_EYFPneg_CREneg_cells")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "EYFP", group.by = "SCT_snn_res.1.1", idents = "EYFPneg_EYFPneg_CREpos_cells", pt.size = 0, y.max = 4) + ggtitle("EYFP EYFPneg_EYFPneg_CREpos_cells")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "EYFP", group.by = "SCT_snn_res.1.1", idents = "EYFPneg_EYFPpos_CREneg_cells", pt.size = 0, y.max = 4) + ggtitle("EYFP EYFPneg_EYFPpos_CREneg_cells")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "EYFP", group.by = "SCT_snn_res.1.1", idents = "EYFPneg_EYFPpos_CREpos_cells", pt.size = 0, y.max = 4) + ggtitle("EYFP EYFPneg_EYFPpos_CREpos_cells")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "EYFP", group.by = "SCT_snn_res.1.1", idents = "EYFPpos_EYFPneg_CREneg_cells", pt.size = 0, y.max = 4) + ggtitle("EYFP EYFPpos_EYFPneg_CREneg_cells")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "EYFP", group.by = "SCT_snn_res.1.1", idents = "EYFPpos_EYFPneg_CREpos_cells", pt.size = 0, y.max = 4) + ggtitle("EYFP EYFPpos_EYFPneg_CREpos_cells")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "EYFP", group.by = "SCT_snn_res.1.1", idents = "EYFPpos_EYFPpos_CREneg_cells", pt.size = 0, y.max = 4) + ggtitle("EYFP EYFPpos_EYFPpos_CREneg_cells")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "EYFP", group.by = "SCT_snn_res.1.1", idents = "EYFPpos_EYFPpos_CREpos_cells", pt.size = 0, y.max = 4) + ggtitle("EYFP EYFPpos_EYFPpos_CREpos_cells")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "CRE", group.by = "SCT_snn_res.1.1", idents = "CREneg_CREneg_CREneg_cells", pt.size = 0, y.max = 4) + ggtitle("CRE CREneg_CREneg_CREneg_cells")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "CRE", group.by = "SCT_snn_res.1.1", idents = "CREneg_CREneg_CREpos_cells", pt.size = 0, y.max = 4) + ggtitle("CRE CREneg_CREneg_CREpos_cells")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "CRE", group.by = "SCT_snn_res.1.1", idents = "CREneg_CREpos_CREneg_cells", pt.size = 0, y.max = 4) + ggtitle("CRE CREneg_CREpos_CREneg_cells")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "CRE", group.by = "SCT_snn_res.1.1", idents = "CREneg_CREpos_CREpos_cells", pt.size = 0, y.max = 4) + ggtitle("CRE CREneg_CREpos_CREpos_cells")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "CRE", group.by = "SCT_snn_res.1.1", idents = "CREpos_CREneg_CREneg_cells", pt.size = 0, y.max = 4) + ggtitle("CRE CREpos_CREneg_CREneg_cells")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "CRE", group.by = "SCT_snn_res.1.1", idents = "CREpos_CREneg_CREpos_cells", pt.size = 0, y.max = 4) + ggtitle("CRE CREpos_CREneg_CREpos_cells")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "CRE", group.by = "SCT_snn_res.1.1", idents = "CREpos_CREpos_CREneg_cells", pt.size = 0, y.max = 4) + ggtitle("CRE CREpos_CREpos_CREneg_cells")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(all.mesenchyme,assay = "RNA", features = "CRE", group.by = "SCT_snn_res.1.1", idents = "CREpos_CREpos_CREpos_cells", pt.size = 0, y.max = 4) + ggtitle("CRE CREpos_CREpos_CREpos_cells")+ 
  stat_summary(fun.y = fun1, geom = "point", shape=95, size=4, colour="red" ) + stat_summary(fun.y = fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 
dev.off()


##Dotplot
#ggplot colors
library(scales)

#extract hex color codes for a plot with three elements in ggplot2 
hex <- hue_pal()(18)

#overlay hex color codes on actual colors
show_col(hex)

#levels(mensenchyme.clus.labelled) <-c("TP","EDC","ICT","DP","DPP","PPP","PC","Ms","LDC")
levels(all.mesenchyme) <-c("5","2","11","6","0","1","9","13","10","3","8","12","4","16","7","14","15","17")
colors <- c("#F8766D","#E88526","#D39200","#B79F00","#93AA00","#5EB300","#00BA38","#00BF74","#00C19F","#00BFC4",
            "#00B9E3","#00ADFA","#619CFF","#AE87FF","#DB72FB","#F564E3","#FF61C3","#FF699C")

markers.to.plot <- c("Top2a","Isl1","Irx3","Shh","Meis1","Msx1","Osr1","Hoxa11","Lum","Kera","Sall1","Scx","Gdf5","Ebf2",
                     "Foxp1","Sox9","Runx2","Col9a3","Col2a1","Foxc1","Irx1","Shox2","Hoxd13", "Pitx1")
DotPlot(object = all.mesenchyme, features = markers.to.plot, cols = colors,split.by = "SCT_snn_res.1.1")

markers.to.plot2 <- c("Top2a","Gsc","Isl1","Irx3","Meis1","Msx1","Grem1","Sall3","Tbx4",
                      "Ebf3","Foxc1","Meg3","Fibin","Lhx2","Etv4",
                      "Scx","Hoxa11","Hoxd9", "Hoxd10","Hoxd13","Shox2")
DotPlot(object = all.mesenchyme, features = markers.to.plot2, cols = colors,split.by = "SCT_snn_res.1.1")

markers.to.plot3 <- c("Top2a","Lin28b","Col2a1","Col1a1","Foxp1","Msx1","Meis2","Irx3","Isl1",
                     "Nfib", "Fibin",
                      "Scx","Hoxa11","Hoxd9", "Hoxd10","Hoxd13","Shox2")
DotPlot(object = all.mesenchyme, features = markers.to.plot3, cols = colors,split.by = "SCT_snn_res.1.1")

markers.to.plot4 <- c("Hist1h2a0","Igdcc3","Asb4","Hmga2","Shox2","Hoxd13","Msx1","Lhx9","Aldh1a2",
                      "Dlk1","Zfhx3","Scx","Lum","Zfhx4","Dcn","Col5a2","Lgals1","Tpm1","Col1a1","Col3a1",
                      "Osr1","Sparc","Mgp","Ebf1","Sfrp2","Wwp2","Sox9","Col2a1","Col9a1","Col9a2","Col11a1",
                      "Cnmd","Matn1")
DotPlot(object = all.mesenchyme, features = markers.to.plot4, cols = colors,split.by = "SCT_snn_res.1.1")

markers.to.plot3 <- c("Top2a","Lin28b","Msx1","Grem1","Sall3","Lhx2","Lhx9","Etv4",
                      "Hoxa11","Hoxd9", "Hoxd10","Hoxd13","Shox2","Lhx2","Lhx9","Etv4",
                      "Gsc","Isl1","Irx3","Meis1",
                      "Foxp1", "Fibin","Scx","Col1a1","Osr1","Lum","Dcn","Nfib",
                      "Gdf5","Sox9","Col2a1","Col9a1","Matn1","Irx1")

##Heatmap
all.mesenchyme.markers <- FindAllMarkers(all.mesenchyme, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

DefaultAssay(all.mesenchyme) <- "SCT"
all.mesenchyme.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5
DoHeatmap(all.mesenchyme, features = top5$gene) + NoLegend()

DoHeatmap(all.mesenchyme, features = markers.to.plot3) + NoLegend()

                     

#### M. Nebulosa package Density Plots UMAP1.1 ####
all.mesenchyme <- readRDS(file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.mesenchyme_UMAP1-1.rds")
all.mesenchyme
summary(all.mesenchyme@meta.data)
DimPlot(all.mesenchyme, reduction = "umap", label = TRUE)

Plot_Density_Custom(seurat_object = all.mesenchyme, features = "Shox2")
Plot_Density_Custom(seurat_object = all.mesenchyme, features = "CRE")
Plot_Density_Custom(seurat_object = all.mesenchyme, features = "Hoxd13")
Plot_Density_Custom(seurat_object = all.mesenchyme, features = "EYFP")

#### N. Scvelo UMAP1.1 ####
all.mesenchyme <- readRDS(file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.mesenchyme_UMAP1-1.rds")
all.mesenchyme
summary(all.mesenchyme@meta.data)
DimPlot(all.mesenchyme, reduction = "umap", label = TRUE)


##ScVelo ALL SAMPLES##
DefaultAssay(all.mesenchyme) <- "RNA"
all.mesenchyme

SaveH5Seurat(all.mesenchyme,filename = "all_clusters.UMAP1_1.h5Seurat")
Convert("all_clusters.UMAP1_1.h5Seurat", dest="h5ad")


#Create an environment in conda and install scvelo
conda create -n en_scvelo
conda activate en_scvelo
pip install -U scvelo


#How to run scvelo
python
import scvelo as scv
adata = scv.read("all_clusters.UMAP1_1.h5ad")
adata

ident_colours = ['#F8766D','#E88526','#D39200','#B79F00','#93AA00','#5EB300','#00BA38','#00BF74','#00C19F','#00BFC4','#00B9E3','#00ADFA','#619CFF','#AE87FF','#DB72FB','#F564E3','#FF61C3','#FF699C']
adata.uns['SCT_snn_res.1.1_colors']=ident_colours

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=10, n_neighbors=30)

scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)

scv.pl.velocity_embedding_stream(adata, basis="umap", color="SCT_snn_res.1.1", palette = ident_colours)

#clusters = ['TP','EDC','ICT','DP','DPP','PPP','PC','Ms','LDC']
#adata.obs['seurat_clusters'].cat.categories = clusters
#scv.pl.velocity_embedding_stream(adata, basis="umap", color="seurat_clusters", palette = ident_colours, legend_loc = 'right margin')

#scv.pl.velocity_embedding_stream(adata, basis="umap", color="seurat_clusters", palette = ident_colours)
#scv.pl.velocity_embedding_stream(adata, basis="umap", color="seurat_clusters", palette = ident_colours,size = 50,legend_fontsize = 15)

#### O. Cell Proportion UMAP 1.1 ####
all.mesenchyme <- readRDS(file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.mesenchyme_UMAP1-1.rds")
all.mesenchyme
summary(all.mesenchyme@meta.data)
DimPlot(all.mesenchyme, reduction = "umap", label = TRUE)


#number of cells per cluster and sample stage
number_cells_percluster <- table(all.mesenchyme$SCT_snn_res.1.1, all.mesenchyme$stage)
number_cells_percluster

#proportion of cells per cluster and sample
propotions_cells_percluster <- prop.table(table(all.mesenchyme$SCT_snn_res.1.1, all.mesenchyme$stage), margin = 2)
propotions_cells_percluster.df <- as.data.frame(propotions_cells_percluster)
t(propotions_cells_percluster)

write.csv(propotions_cells_percluster.df,"UMAP_MERGED_cell_proportions_all_clusters.csv")

ggplot(propotions_cells_percluster.df, aes(fill= Var2, color = Var2, x=Var1,y=Freq)) + 
  geom_bar(position="dodge", stat = "identity") + theme_bw() 

ggplot(propotions_cells_percluster.df, aes(fill= Var1, x=Var2,y=Freq)) + 
  geom_bar(position="fill", stat = "identity") + 
  geom_text(aes(label = sprintf("%0.2f", round(Freq, 2))), color="black", size=2,position = position_stack(vjust = 0.5))

#proportion of cells per cluster and sample
propotions_cells_percluster <- prop.table(table(all.mesenchyme$stage, all.mesenchyme$SCT_snn_res.1.1), margin = 2)
propotions_cells_percluster.df <- as.data.frame(propotions_cells_percluster)
t(propotions_cells_percluster)

write.csv(propotions_cells_percluster.df,"UMAP_MERGED_cell_proportions_all_clusters.csv")

ggplot(propotions_cells_percluster.df, aes(fill= Var2, color = Var2, x=Var1,y=Freq)) + 
  geom_bar(position="dodge", stat = "identity") + theme_bw() 

ggplot(propotions_cells_percluster.df, aes(fill= Var1, x=Var2,y=Freq)) + 
  geom_bar(position="fill", stat = "identity") + 
  geom_text(aes(label = sprintf("%0.2f", round(Freq, 2))), color="black", size=2,position = position_stack(vjust = 0.5))

#### P. Merging similar clusters and cluster naming ####
all.mesenchyme <- readRDS(file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.mesenchyme_UMAP1-1.rds")
all.mesenchyme
summary(all.mesenchyme@meta.data)
DimPlot(all.mesenchyme, reduction = "umap", label = TRUE)

#Potential cluster identification
#These are the names of the clusters based on what we discuss the other day:
#0+6+11: Late Proximal Progenitors (LPP)
#1: Proximal Conetive Tissue (PCT)
#2: Early Proximal Progenitors (EPP)
#3: Distal Proliferative Progenitors (DPP)
#4: Interdigit Mesenchyme (IM)
#5: Limb Progenitors (LP)
#7: Early Proximal Condensation (EPC)
#8: Mesopodium (Ms)
#9: Tendon Progenitors (TP)
#10 + 13: Irregular Connective Tissue (ICT)
#12: Distal Progenitors (DP)
#14: Late Digit Condensation (LDC)
#15: Proximal Growth Plate (PGP)
#16: Early Digit Condensatio (EDC)
#17: Proximal Condensation (PC)

all.mesenchyme@meta.data$merged.MES.clusters <- all.mesenchyme@meta.data$SCT_snn_res.1.1
levels(all.mesenchyme$merged.MES.clusters) <- c("LPP","PCT","EPP","DPP","IM","LP","LPP","EPC",
                                                "Ms","TP","ICT","LPP","DP",
                                                "ICT","LDC","PGP","EDC","PC")

DimPlot(all.mesenchyme, reduction = "umap",label = TRUE, label.size = 3, group.by = "merged.MES.clusters", repel = TRUE)
DimPlot(all.mesenchyme, reduction = "umap",label = FALSE, group.by = "merged.MES.clusters", split.by = "stage")

levels(all.mesenchyme)
Idents(all.mesenchyme) <- "merged.MES.clusters"
levels(all.mesenchyme)

saveRDS(all.mesenchyme, "/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.mesenchyme_UMAP1.1MERGED.rds")

#### Q. Marker genes MERGED UMAP1.1 ####
all.mesenchyme <- readRDS(file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.mesenchyme_UMAP1.1MERGED.rds")
all.mesenchyme
summary(all.mesenchyme@meta.data)
DimPlot(all.mesenchyme, reduction = "umap", label = TRUE)

DefaultAssay(all.mesenchyme) <- "RNA"
all.mesenchyme
levels(all.mesenchyme)

LPP.markers <-  FindMarkers(all.mesenchyme, ident.1 = "LPP",
                           ident.2 = c("PCT","EPP","DPP","IM","LP","EPC","Ms","TP","ICT","DP","LDC","PGP","EDC","PC"),
                           only.pos =TRUE,
                           logfc.threshold = 0.5,
                           pseudocount.use = 0,
                           min.diff.pct=0.1)
head (LPP.markers,15)
write.csv(LPP.markers,file=("Late.Proximal.Progenitors.markers.csv"))


PCT.markers <-  FindMarkers(all.mesenchyme, ident.1 = "PCT",
                            ident.2 = c("LPP","EPP","DPP","IM","LP","EPC","Ms","TP","ICT","DP","LDC","PGP","EDC","PC"),
                            only.pos =TRUE,
                            logfc.threshold = 0.5,
                            pseudocount.use = 0,
                            min.diff.pct=0.1)
head (PCT.markers,15)
write.csv(PCT.markers,file=("Proximal.Connective.Tissue.markers.csv"))


EPP.markers <-  FindMarkers(all.mesenchyme, ident.1 = "EPP",
                            ident.2 = c("LPP","PCT","DPP","IM","LP","EPC","Ms","TP","ICT","DP","LDC","PGP","EDC","PC"),
                            only.pos =TRUE,
                            logfc.threshold = 0.5,
                            pseudocount.use = 0,
                            min.diff.pct=0.1)
head (EPP.markers,15)
write.csv(EPP.markers,file=("Early.Proximal.Progenitors.markers.csv"))


DPP.markers <-  FindMarkers(all.mesenchyme, ident.1 = "DPP",
                            ident.2 = c("LPP","PCT","EPP","IM","LP","EPC","Ms","TP","ICT","DP","LDC","PGP","EDC","PC"),
                            only.pos =TRUE,
                            logfc.threshold = 0.5,
                            pseudocount.use = 0,
                            min.diff.pct=0.1)
head (DPP.markers,15)
write.csv(DPP.markers,file=("Distal.Proliferative.Progenitors.markers.csv"))


IM.markers <-  FindMarkers(all.mesenchyme, ident.1 = "IM",
                            ident.2 = c("LPP","PCT","DPP","EPP","LP","EPC","Ms","TP","ICT","DP","LDC","PGP","EDC","PC"),
                            only.pos =TRUE,
                            logfc.threshold = 0.5,
                            pseudocount.use = 0,
                            min.diff.pct=0.1)
head (IM.markers,15)
write.csv(IM.markers,file=("Interdigit.Mesenchyme.markers.csv"))


LP.markers <-  FindMarkers(all.mesenchyme, ident.1 = "LP",
                           ident.2 = c("LPP","PCT","DPP","EPP","IM","EPC","Ms","TP","ICT","DP","LDC","PGP","EDC","PC"),
                           only.pos =TRUE,
                           logfc.threshold = 0.5,
                           pseudocount.use = 0,
                           min.diff.pct=0.1)
head (LP.markers,15)
write.csv(LP.markers,file=("Limb.Progenitors.markers.csv"))


EPC.markers <-  FindMarkers(all.mesenchyme, ident.1 = "EPC",
                           ident.2 = c("LPP","PCT","DPP","EPP","IM","LP","Ms","TP","ICT","DP","LDC","PGP","EDC","PC"),
                           only.pos =TRUE,
                           logfc.threshold = 0.5,
                           pseudocount.use = 0,
                           min.diff.pct=0.1)
head (EPC.markers,15)
write.csv(EPC.markers,file=("Early.Proximal.Condensation.markers.csv"))


Ms.markers <-  FindMarkers(all.mesenchyme, ident.1 = "Ms",
                            ident.2 = c("LPP","PCT","DPP","EPP","IM","LP","EPC","TP","ICT","DP","LDC","PGP","EDC","PC"),
                            only.pos =TRUE,
                            logfc.threshold = 0.5,
                            pseudocount.use = 0,
                            min.diff.pct=0.1)
head (Ms.markers,15)
write.csv(Ms.markers,file=("Mesopodium.markers.csv"))


TP.markers <-  FindMarkers(all.mesenchyme, ident.1 = "TP",
                           ident.2 = c("LPP","PCT","DPP","EPP","IM","LP","EPC","Ms","ICT","DP","LDC","PGP","EDC","PC"),
                           only.pos =TRUE,
                           logfc.threshold = 0.5,
                           pseudocount.use = 0,
                           min.diff.pct=0.1)
head (TP.markers,15)
write.csv(TP.markers,file=("Tendon.progenitors.markers.csv"))


ICT.markers <-  FindMarkers(all.mesenchyme, ident.1 = "ICT",
                           ident.2 = c("LPP","PCT","DPP","EPP","IM","LP","EPC","Ms","TP","DP","LDC","PGP","EDC","PC"),
                           only.pos =TRUE,
                           logfc.threshold = 0.5,
                           pseudocount.use = 0,
                           min.diff.pct=0.1)
head (ICT.markers,15)
write.csv(ICT.markers,file=("Irregular.Connective.Tissue.markers.csv"))


DP.markers <-  FindMarkers(all.mesenchyme, ident.1 = "DP",
                            ident.2 = c("LPP","PCT","DPP","EPP","IM","LP","EPC","Ms","TP","ICT","LDC","PGP","EDC","PC"),
                            only.pos =TRUE,
                            logfc.threshold = 0.5,
                            pseudocount.use = 0,
                            min.diff.pct=0.1)
head (DP.markers,15)
write.csv(DP.markers,file=("Distal.Progenitors.markers.csv"))


LDC.markers <-  FindMarkers(all.mesenchyme, ident.1 = "LDC",
                           ident.2 = c("LPP","PCT","DPP","EPP","IM","LP","EPC","Ms","TP","ICT","DP","PGP","EDC","PC"),
                           only.pos =TRUE,
                           logfc.threshold = 0.5,
                           pseudocount.use = 0,
                           min.diff.pct=0.1)
head (LDC.markers,15)
write.csv(LDC.markers,file=("Late.Digit.Condensation.markers.csv"))


PGP.markers <-  FindMarkers(all.mesenchyme, ident.1 = "PGP",
                            ident.2 = c("LPP","PCT","DPP","EPP","IM","LP","EPC","Ms","TP","ICT","DP","LDC","EDC","PC"),
                            only.pos =TRUE,
                            logfc.threshold = 0.5,
                            pseudocount.use = 0,
                            min.diff.pct=0.1)
head (PGP.markers,15)
write.csv(PGP.markers,file=("Proximal.Growth.Plate.markers.csv"))


EDC.markers <-  FindMarkers(all.mesenchyme, ident.1 = "EDC",
                            ident.2 = c("LPP","PCT","DPP","EPP","IM","LP","EPC","Ms","TP","ICT","DP","LDC","PGP","PC"),
                            only.pos =TRUE,
                            logfc.threshold = 0.5,
                            pseudocount.use = 0,
                            min.diff.pct=0.1)
head (EDC.markers,15)
write.csv(EDC.markers,file=("Early.Digit.Condensation.markers.csv"))



PC.markers <-  FindMarkers(all.mesenchyme, ident.1 = "PC",
                            ident.2 = c("LPP","PCT","DPP","EPP","IM","LP","EPC","Ms","TP","ICT","DP","LDC","PGP","EDC"),
                            only.pos =TRUE,
                            logfc.threshold = 0.5,
                            pseudocount.use = 0,
                            min.diff.pct=0.1)
head (PC.markers,15)
write.csv(PC.markers,file=("Proximal.Condensation.markers.csv"))


#### R. Cell Proportion UMAP 1.1 MERGED clusters ####
all.mesenchyme <- readRDS(file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.mesenchyme_UMAP1.1MERGED.rds")
all.mesenchyme
summary(all.mesenchyme@meta.data)
DimPlot(all.mesenchyme, reduction = "umap", label = TRUE)

cell.num <- table(Idents(all.mesenchyme))
clusterLabels <- paste(names(cell.num), paste0("(", cell.num, ")"))
clusterBreaks <- names(cell.num)

DimPlot(all.mesenchyme, reduction = "umap", label=TRUE) + scale_colour_discrete(breaks = clusterBreaks,labels = clusterLabels) + labs(title = "UMAP")

#number of cells per cluster
table(all.mesenchyme$merged.MES.clusters)

#proportion of cells per cluster
prop.table(table(Idents(all.mesenchyme)))

#number of cells per stage 
total_cells <- table(all.mesenchyme$stage)
total_cells

#number of cells per cluster and sample stage
number_cells_percluster <- table(all.mesenchyme$merged.MES.clusters, all.mesenchyme$stage)
number_cells_percluster

## DPA = Differential Proportion Analysis ##
#Install functions needed for the analysis
source("~/Dropbox/DPA/diffprop_functions.R");

## Read in file of counts of cells in each population (Cluster) across conditions (datasets)
obs.counts <- t(number_cells_percluster)

## Run an example using error (p) of 0.1 and with 100,000 iterations
tip.exp <- generateNull(obs.counts, n=100000, p=0.1)    # Generate the null distribution based on sampling

### P-value tests for Sham vs MI-Day 3
pp <- two.class.test(obs.counts, tip.exp, cond.control="E105", cond.treatment="E115",to.plot=T)
pp
write.csv(pp,"UMAP_Mesenchyme_MERGED_E105vsE115-pvalues.csv")

pp1 <- two.class.test(obs.counts, tip.exp, cond.control="E115", cond.treatment="E125",to.plot=T)
pp1
write.csv(pp1,"UMAP_Mesenchyme_MERGED_E115vsE125-pvalues.csv")

pp2 <- two.class.test(obs.counts, tip.exp, cond.control="E125", cond.treatment="E135",to.plot=T)
pp2
write.csv(pp2,"UMAP_Mesenchyme_MERGED_E125vsE135-pvalues.csv")

pp3 <- two.class.test(obs.counts, tip.exp, cond.control="E105", cond.treatment="E125",to.plot=T)
pp3
write.csv(pp3,"UMAP_Mesenchyme_MERGED_E105vsE125-pvalues.csv")

pp4 <- two.class.test(obs.counts, tip.exp, cond.control="E105", cond.treatment="E135",to.plot=T)
pp4
write.csv(pp4,"UMAP_Mesenchyme_MERGED_E105vsE135-pvalues.csv")

pp5 <- two.class.test(obs.counts, tip.exp, cond.control="E115", cond.treatment="E135",to.plot=T)
pp5
write.csv(pp5,"UMAP_Mesenchyme_MERGED_E115vsE135-pvalues.csv")

#number of cells per cluster and sample stage
number_cells_percluster <- table(all.mesenchyme$merged.MES.clusters, all.mesenchyme$stage)
number_cells_percluster

#proportion of cells per cluster and sample
propotions_cells_percluster <- prop.table(table(all.mesenchyme$merged.MES.clusters, all.mesenchyme$stage), margin = 2)
propotions_cells_percluster.df <- as.data.frame(propotions_cells_percluster)
t(propotions_cells_percluster)

write.csv(propotions_cells_percluster.df,"UMAP_Mesenchyme_MERGED_cell_proportions.csv")

ggplot(propotions_cells_percluster.df, aes(fill= Var2, color = Var2, x=Var1,y=Freq)) + 
  geom_bar(position="dodge", stat = "identity") + theme_bw() 


propotions_cells_percluster2 <- prop.table(table(all.mesenchyme$stage, all.mesenchyme$merged.MES.clusters), margin = 2)
propotions_cells_percluster2.df <- as.data.frame(propotions_cells_percluster2)
t(propotions_cells_percluster)

ggplot(propotions_cells_percluster2.df, aes(fill= Var1, x=Var2,y=Freq)) + 
  geom_bar(position="fill", stat = "identity") + 
  geom_text(aes(label = sprintf("%0.2f", round(Freq, 2))), color="black", size=2,position = position_stack(vjust = 0.5))


propotions_cells_percluster2B <- prop.table(table(all.mesenchyme$stage, all.mesenchyme$merged.MES.clusters), margin = 2)
propotions_cells_percluster2B.df <- as.data.frame(propotions_cells_percluster2B)
t(propotions_cells_percluster2B)

positions <- c("LP", "EPP", "LPP", "DPP", "Ms", "EPC", "PCT", "DP", "EDC", "PGP", "TP", "LDC", "ICT", "PC","IM")
propotions_cells_percluster2B.df$Var2 <- factor(propotions_cells_percluster2B.df$Var2, levels = positions)

ggplot(propotions_cells_percluster2B.df, aes(fill= Var1, x=Var2,y=Freq)) + 
  geom_bar(position="fill", stat = "identity") + 
  geom_text(aes(label = sprintf("%0.2f", round(Freq, 2))), color="black", size=2,position = position_stack(vjust = 0.5))


propotions_cells_percluster3 <- prop.table(table(all.mesenchyme$gene_class, all.mesenchyme$merged.MES.clusters), margin = 2)
propotions_cells_percluster3.df <- as.data.frame(propotions_cells_percluster3)
t(propotions_cells_percluster3)

ggplot(propotions_cells_percluster3.df, aes(fill= Var1, x=Var2,y=Freq)) + 
  geom_bar(position="fill", stat = "identity") + 
  geom_text(aes(label = sprintf("%0.2f", round(Freq, 2))), color="black", size=2,position = position_stack(vjust = 0.5))


propotions_cells_percluster4 <- prop.table(table(all.mesenchyme$gene_class, all.mesenchyme$stage), margin = 2)
propotions_cells_percluster4.df <- as.data.frame(propotions_cells_percluster4)
t(propotions_cells_percluster4)

ggplot(propotions_cells_percluster4.df, aes(fill= Var1, x=Var2,y=Freq)) + 
  geom_bar(position="fill", stat = "identity") + 
  geom_text(aes(label = sprintf("%0.2f", round(Freq, 2))), color="black", size=2,position = position_stack(vjust = 0.5))

all.mesenchyme$geneclass.stage <- paste(all.mesenchyme$gene_class,all.mesenchyme$stage, sep = "_")
summary(all.mesenchyme@meta.data)

propotions_cells_percluster5 <- prop.table(table(all.mesenchyme$geneclass.stage, all.mesenchyme$merged.MES.clusters), margin = 2)
propotions_cells_percluster5.df <- as.data.frame(propotions_cells_percluster5)
t(propotions_cells_percluster5)

ggplot(propotions_cells_percluster5.df, aes(fill= Var1, x=Var2,y=Freq)) + 
  geom_bar(position="fill", stat = "identity") + 
  geom_text(aes(label = sprintf("%0.2f", round(Freq, 2))), color="black", size=2,position = position_stack(vjust = 0.5))


write.csv(propotions_cells_percluster5.df,"UMAP_Mesenchyme_MERGED_cell_proportions_geneclass_stage_percluster.csv")


#### S. Nebulosa package Density Plots MERGED UMAP1.1 ####
all.mesenchyme <- readRDS(file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.mesenchyme_UMAP1.1MERGED.rds")
all.mesenchyme
summary(all.mesenchyme@meta.data)
DimPlot(all.mesenchyme, reduction = "umap", label = TRUE)

DefaultAssay(all.mesenchyme) <- "RNA"
all.mesenchyme

Plot_Density_Custom(seurat_object = all.mesenchyme, features = "Shox2")
Plot_Density_Custom(seurat_object = all.mesenchyme, features = "CRE")
Plot_Density_Custom(seurat_object = all.mesenchyme, features = "Hoxd13")
Plot_Density_Custom(seurat_object = all.mesenchyme, features = "EYFP")

#### T.1 Some plots MERGED UMAP1.1: Dimplot +Scatter Plot ####
all.mesenchyme <- readRDS(file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.mesenchyme_UMAP1.1MERGED.rds")
all.mesenchyme
summary(all.mesenchyme@meta.data)
DimPlot(all.mesenchyme, reduction = "umap", label = TRUE)

DefaultAssay(all.mesenchyme) <- "RNA"
all.mesenchyme
levels(all.mesenchyme)

DimPlot(all.mesenchyme, reduction = "umap", label = TRUE)

#Dimplot
Idents(all.mesenchyme) <- "gene_class"
DimPlot(all.mesenchyme, reduction = "umap", label = FALSE)

#Dimplot by stage
Idents(all.mesenchyme) <- "stage"
levels(all.mesenchyme)

DimPlot(all.mesenchyme, reduction = "umap", label = FALSE)


## Scatter plot ##
FeatureScatter(all.mesenchyme, feature1 = "EYFP", feature2 = "Shox2")
FeatureScatter(all.mesenchyme, feature1 = "CRE", feature2 = "Shox2")

Idents(all.mesenchyme) <- "gene_class"
levels(all.mesenchyme)
levels(all.mesenchyme) <-c("EYFPneg_Shox2neg_CREneg_cells","EYFPneg_Shox2neg_CREpos_cells","EYFPneg_Shox2pos_CREneg_cells",
                                        "EYFPneg_Shox2pos_CREpos_cells","EYFPpos_Shox2neg_CREneg_cells","EYFPpos_Shox2neg_CREpos_cells",
                                        "EYFPpos_Shox2pos_CREneg_cells","EYFPpos_Shox2pos_CREpos_cells")

FeatureScatter(all.mesenchyme, feature1 = "EYFP", feature2 = "Shox2")
FeatureScatter(all.mesenchyme, feature1 = "CRE", feature2 = "Shox2")

FeatureScatter(all.mesenchyme,slot = "counts", feature1 = "EYFP", feature2 = "CRE")
FeatureScatter(all.mesenchyme,slot = "data", feature1 = "EYFP", feature2 = "CRE")

#### T.2 Some plots MERGED UMAP1.1: Feature Plot + DotPlot ####
all.mesenchyme <- readRDS(file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.mesenchyme_UMAP1.1MERGED.rds")
all.mesenchyme
summary(all.mesenchyme@meta.data)
DimPlot(all.mesenchyme, reduction = "umap", label = TRUE)

## Feature Plot ##
DefaultAssay(all.mesenchyme) <- "RNA"
all.mesenchyme

#pdf(file=paste0("all.mesenchyme_other_Features.pdf"))
FeaturePlot(object = all.mesenchyme, features = "Acan" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Shox2" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Col1a1" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Col2a1" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Col10a1" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Ihh" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Bmp4" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Bmp2" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Prrx1" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Shh" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "Rsrc1" ,cols = c("grey95","darkred"), combine = FALSE)
FeaturePlot(object = all.mesenchyme, features = "CRE" ,cols = c("grey95","darkred"), combine = FALSE)

FeaturePlot(object = all.mesenchyme, features = "Smarca2" ,cols = c("grey95","darkred"), split.by = "stage", combine = TRUE)
FeaturePlot(object = all.mesenchyme, features = "Smarca4" ,cols = c("grey95","darkred"), split.by = "stage", combine = TRUE)
FeaturePlot(object = all.mesenchyme, features = "Col1a1" ,cols = c("grey95","darkred"), split.by = "stage", combine = TRUE)
dev.off()

FeaturePlot(all.mesenchyme, features = c("Hoxd13", "Shox2"), blend = TRUE, 
            blend.threshold = 0.5, order = TRUE)

## Dotplot
DefaultAssay(all.mesenchyme) <- "RNA"
all.mesenchyme
markers.to.plot2 <- c("Hoxd13","Shox2")
DotPlot(object = all.mesenchyme, features = markers.to.plot2)


levels(all.mesenchyme)
##Dotplot
#ggplot colors
library(scales)

#extract hex color codes for a plot with three elements in ggplot2 
hex <- hue_pal()(15)
#overlay hex color codes on actual colors
show_col(hex)

colors <- c("#F8766D","#E58700","#C99800","#A3A500","#6BB100","#00BA38","#00BF7D","#00C0AF","#00BCD8","#00B0F6",
            "#619CFF","#B983FF","#E76BF3","#FD61D1","#FF67A4")
markers.to.plot <- c("Irx3","Osr1","Sall4","Evx1","Dbx2","Tfap2c","Sox9","Gdf5","Meox1","Kera","Msx1","Acan",
                     "Col9a3","Irx1","Ihh","CRE","Shox2","Hoxd13","EYFP")
DotPlot(object = all.mesenchyme, features = markers.to.plot, cols = colors,split.by = "merged.MES.clusters")


#not used at the moment, I should reorganize the markers
levels(all.mesenchyme) <- c("LDC","EDC","IM","DP","DPP","Ms","LP", "EPP", "LPP", "PCT","ICT","TP", "EPC", "PGP", "PC")
colors <- c("#B983FF","#FD61D1","#6BB100","#619CFF","#A3A500","#00C0AF","#00BA38","#C99800","#F8766D",
            "#E58700","#00B0F6","#00BCD8","#00BF7D","#E76BF3","#FF67A4")

markers.to.plot <- c("Irx3","Osr1","Sall4","Evx1","Dbx2","Tfap2c","Sox9","Gdf5","Meox1","Kera","Msx1","Acan",
                     "Col9a3","Irx1","Ihh","CRE","Shox2","Hoxd13","EYFP")
DotPlot(object = all.mesenchyme, features = markers.to.plot, cols = colors,split.by = "merged.MES.clusters")

#### T.3 Some plots MERGED UMAP1.1: gene expression distribution ####
all.mesenchyme <- readRDS(file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.mesenchyme_UMAP1.1MERGED.rds")
all.mesenchyme
summary(all.mesenchyme@meta.data)
DimPlot(all.mesenchyme, reduction = "umap", label = TRUE)

df1 <- all.mesenchyme@assays$RNA@data[,]
df2 <- data.frame(df1)

genes <- c("Shox2","Hoxd13","CRE","EYFP","Msx1","Shh")

df3 <- df2[genes,]
df3 <- as.data.frame(t(df3))
df3$cells <- rownames(df3)

list2 <- lapply(colnames(df3), function(x) {
  df3[df3[[x]] > 0,"cells"]
})
names(list2) <- colnames(df3)
list3 <-list2[-5]

library("UpSetR")
upset(fromList(list3),sets = c("Shox2","Hoxd13","CRE","EYFP"),
      keep.order = TRUE, nintersects = NA)

upset(fromList(list3),sets = c("Shox2","Hoxd13","CRE","EYFP"),
      keep.order = TRUE, empty.intersections ="on", nintersects = NA)

### Shox2_Hoxd13 pos and neg classification 
genes2 <- c("Shox2","Hoxd13")
df4 <- df2[genes2,]
df4 <- as.data.frame(t(df4))
df4$cells <- rownames(df4)

list4 <- lapply(colnames(df4), function(x) {
  df4[df4[[x]] > 0,"cells"]
})
names(list4) <- colnames(df4)
list5 <-list4[-3]

library("UpSetR")
upset(fromList(list5),sets = c("Shox2","Hoxd13"),
      keep.order = TRUE, nintersects = NA)

upset(fromList(list5),sets = c("Shox2","Hoxd13"),
      keep.order = TRUE, empty.intersections ="on", nintersects = NA)

### Shox2_Hoxd13 pos and neg classification 
genes2 <- c("Shox2","Hoxd13")
df4 <- df2[genes2,]
df4 <- as.data.frame(t(df4))

create_label_column <- function(data) {
  # Check conditions for "double positive," "one_positive," and "two_positive"
  double_positive_condition <- (data$Shox2 > 0) & (data$Hoxd13 > 0)
  one_positive_condition <- (data$Shox2 > 0) & !(data$Hoxd13 > 0)
  two_positive_condition <- !(data$Shox2 > 0) & (data$Hoxd13 > 0)
  
  # Create a new column based on the conditions
  data$Shx_Hox <- ifelse(double_positive_condition, "Shx_Hox_pos",
                         ifelse(one_positive_condition, "Shx_pos",
                                ifelse(two_positive_condition, "Hox_pos", "negative")))
  
  # Return the modified data frame
  return(data)
}

df5 <- create_label_column(df4)

head(df5)

all.mesenchyme[[colnames(x = df5)]] <- df5
head(all.mesenchyme@meta.data)

Idents(all.mesenchyme) <- "Shx_Hox"
DimPlot(all.mesenchyme, reduction = "umap", label = TRUE)

Idents(all.mesenchyme) <- "merged.MES.clusters"
VlnPlot(object = all.mesenchyme, assay = "RNA", features = "Shox2", split.by = "Shx_Hox", 
        pt.size = 0, y.max = 4,add.noise = TRUE ) + ggtitle("Shox2") 


propotions_cells_percluster3 <- prop.table(table(all.mesenchyme$Shx_Hox, all.mesenchyme$merged.MES.clusters), margin = 2)
propotions_cells_percluster3.df <- as.data.frame(propotions_cells_percluster3)
t(propotions_cells_percluster3)

ggplot(propotions_cells_percluster3.df, aes(fill= Var1, x=Var2,y=Freq)) + 
  geom_bar(position="fill", stat = "identity") + 
  geom_text(aes(label = sprintf("%0.2f", round(Freq, 2))), color="black", size=2,position = position_stack(vjust = 0.5))

### EYFP Classification ###
all.mesenchyme

EYFP.expression <- all.mesenchyme@assays$RNA@data["EYFP",]
dfEYFP <- data.frame(EYFP.expression)
dfEYFP$stage <- row.names(dfEYFP)
dfEYFP$stage <- gsub("_[A-Z]+$","", dfEYFP$stage)
dfEYFP$stage <- gsub("HL_","", dfEYFP$stage)
table(dfEYFP$stage)

#Using y=..density.. scales the histograms so the area under each is 1
ggplot(dfEYFP, aes(EYFP.expression,fill=stage, color=stage)) + 
  geom_histogram(aes(y=..density..),position = "dodge", alpha=0.5) + 
  geom_density(alpha=0.05,trim=TRUE) +
  scale_y_continuous(trans = 'sqrt') +
  ggtitle("EYFP global")

ggplot(dfEYFP, aes(EYFP.expression,fill=stage, color=stage)) + 
  geom_density(alpha=0.05,trim=TRUE) +
  scale_y_continuous(trans = 'sqrt') +
  ggtitle("EYFP global")


dfEYFP_E105 <- dfEYFP[dfEYFP$stage == "E105",]

ggplot(dfEYFP_E105, aes(EYFP.expression,fill=stage, color=stage)) + 
  geom_density(alpha=0.05,trim=TRUE) +
  scale_y_continuous(trans = 'sqrt') +
  ggtitle("EYFP global")



### EYFP pos and neg classification 
df1 <- all.mesenchyme@assays$RNA@data[,]
df2 <- data.frame(df1)

genes2 <- c("EYFP")
df4 <- df2[genes2,]
df4 <- as.data.frame(t(df4))

create_label_column2 <- function(data) {
  # Check conditions for "double positive," "one_positive," and "two_positive"
  one_positive_condition <- (data$EYFP > 0) & (data$EYFP < 0.5) 
  two_positive_condition <- (data$EYFP >= 0.5) 
  
  # Create a new column based on the conditions
  data$EYFPclass <- ifelse(one_positive_condition, "low_EYFP",
                           ifelse(two_positive_condition, "High_EYFP", "negative"))
  
  # Return the modified data frame
  return(data)
}

df5 <- create_label_column2(df4)

head(df5)

all.mesenchyme[[colnames(x = df5)]] <- df5
head(all.mesenchyme@meta.data)

Idents(all.mesenchyme) <- "EYFPclass"
DimPlot(all.mesenchyme, reduction = "umap", label = TRUE)

Idents(all.mesenchyme) <- "merged.MES.clusters"
VlnPlot(object = all.mesenchyme, assay = "RNA", features = "EYFP", split.by = "EYFPclass", 
        pt.size = 0, y.max = 4,add.noise = TRUE ) + ggtitle("EYFP") 

Idents(all.mesenchyme) <- "gene_class"
VlnPlot(object = all.mesenchyme, assay = "RNA", features = "EYFP", split.by = "EYFPclass", 
        pt.size = 0, y.max = 4,add.noise = TRUE ) + ggtitle("EYFP") 


propotions_cells_percluster3 <- prop.table(table(all.mesenchyme$EYFPclass, all.mesenchyme$merged.MES.clusters), margin = 2)
propotions_cells_percluster3.df <- as.data.frame(propotions_cells_percluster3)
t(propotions_cells_percluster3)

ggplot(propotions_cells_percluster3.df, aes(fill= Var1, x=Var2,y=Freq)) + 
  geom_bar(position="fill", stat = "identity") + 
  geom_text(aes(label = sprintf("%0.2f", round(Freq, 2))), color="black", size=2,position = position_stack(vjust = 0.5))


Idents(all.mesenchyme) <- "EYFPclass"
lowEYFP <- subset(all.mesenchyme, idents = "low_EYFP")

low.EYFP.expression <- lowEYFP@assays$RNA@data["EYFP",]
dflowEYFP <- data.frame(low.EYFP.expression)
dflowEYFP$stage <- row.names(dflowEYFP)
dflowEYFP$stage <- gsub("_[A-Z]+$","", dflowEYFP$stage)
dflowEYFP$stage <- gsub("HL_","", dflowEYFP$stage)
table(dflowEYFP$stage)

g1 <- ggplot(dflowEYFP, aes(low.EYFP.expression,fill=stage, color=stage)) + 
  geom_density(alpha=0.05,trim=TRUE) +
  scale_y_continuous(trans = 'sqrt') +
  ggtitle("lowEYFP global")

g1


Idents(all.mesenchyme) <- "EYFPclass"
HighEYFP <- subset(all.mesenchyme, idents = "High_EYFP")

High.EYFP.expression <- HighEYFP@assays$RNA@data["EYFP",]
dfHighEYFP <- data.frame(High.EYFP.expression)
dfHighEYFP$stage <- row.names(dfHighEYFP)
dfHighEYFP$stage <- gsub("_[A-Z]+$","", dfHighEYFP$stage)
dfHighEYFP$stage <- gsub("HL_","", dfHighEYFP$stage)
table(dfHighEYFP$stage)

g2 <- ggplot(dfHighEYFP, aes(High.EYFP.expression,fill=stage, color=stage)) + 
  geom_density(alpha=0.05,trim=TRUE) +
  scale_y_continuous(trans = 'sqrt') +
  ggtitle("HighEYFP global")

g2

#### T.4 Some plots MERGED UMAP1.1: VlnPlot Plot ####
all.mesenchyme <- readRDS(file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.mesenchyme_UMAP1.1MERGED.rds")
all.mesenchyme
summary(all.mesenchyme@meta.data)
DimPlot(all.mesenchyme, reduction = "umap", label = TRUE)

## Violin Plot ##
DefaultAssay(all.mesenchyme) <- "RNA"
all.mesenchyme

#split by clusters

#These are the different Vlnplot options that we have tested
VlnPlot(object = all.mesenchyme, assay = "RNA", features = "Shox2", group.by = "merged.MES.clusters", 
        pt.size = 0, y.max = 4,add.noise = TRUE ) + ggtitle("Shox2") 

ggsave("Vlnplot_Shox2_clusters_noiseT.png", height = 5, width = 5)

#ggsave("Vlnplot_Shox2_clusters.eps", height = 3, width = 4)

VlnPlot(object = all.mesenchyme, assay = "RNA", features = "Shox2", group.by = "merged.MES.clusters", 
        pt.size = 0, y.max = 4,add.noise = FALSE ) + ggtitle("Shox2") 

ggsave("Vlnplot_Shox2_clusters_noiseF.png", height = 5, width = 5)

VlnPlot(object = all.mesenchyme, assay = "RNA", features = "Shox2", group.by = "merged.MES.clusters", 
        pt.size = 0, y.max = 4, add.noise = FALSE) + ggtitle("Shox2") + geom_boxplot()

ggsave("Vlnplot_Shox2_clusters_noiseF_boxplot.png", height = 5, width = 5)

VlnPlot(object = all.mesenchyme, assay = "RNA", features = "Shox2", group.by = "merged.MES.clusters", 
        pt.size = 0, y.max = 4, add.noise = TRUE) + ggtitle("Shox2") + geom_boxplot()

ggsave("Vlnplot_Shox2_clusters_noiseT_boxplot.png", height = 5, width = 5)

VlnPlot(object = all.mesenchyme, assay = "RNA", features = c("Shox2","Hoxd13"), group.by = "merged.MES.clusters", 
        pt.size = 0, y.max = 4, add.noise = TRUE, stack = TRUE) 

ggsave("Vlnplot_Shox2_Hoxd13_clusters_noiseT_STACK.png", height = 5, width = 5)
ggsave("Vlnplot_Shox2_Hoxd13_clusters_noiseT_STACK.eps", height = 5, width = 5)


## We finally choose to use the Noise option = TRUE
DefaultAssay(all.mesenchyme) <- "RNA"
all.mesenchyme

#average expression values are NOT normalized using the log1p. 
#Indeed what the function does is undo the log1p for each cell in the data slot
AverageExpression(all.mesenchyme, features = "CRE", group.by = "merged.MES.clusters",assays = "RNA", slot = "data")

avg.Shox2 <- AverageExpression(all.mesenchyme, features = "Shox2", group.by = "merged.MES.clusters",assays = "RNA", slot = "data")
write.csv(avg.Shox2,"Average_expression_Shox2_RNAdata_mesenchyme_nolog1p.csv")

avg.Hoxd13 <- AverageExpression(all.mesenchyme, features = "Hoxd13", group.by = "merged.MES.clusters",assays = "RNA", slot = "data")
write.csv(avg.Hoxd13,"Average_expression_Hoxd13_RNAdata_mesenchyme_nolog1p.csv")


#I can test that with the following functions
fun1 = function (x) {log1p(mean(x = expm1(x = x)))}
fun2 = function (x) {mean(x = expm1(x = x))}

#within the LPP population as an example
LPP <- subset(x = all.mesenchyme, idents = "LPP")
CRE.expression <- LPP@assays$RNA@data["CRE",]
LPP.CRE <-as.data.frame(CRE.expression)

#including log1p normalisation
fun1(LPP.CRE$CRE.expression)

#not doing normalisation
fun2(LPP.CRE$CRE.expression)

#for the violin plot I dediced to plot the average expression log1p normalised since the values plot a
VlnPlot(object = all.mesenchyme, assay = "RNA", features = "CRE", group.by = "merged.MES.clusters", 
        pt.size = 0, y.max = 4, slot = "data",add.noise = TRUE ) + ggtitle("CRE")  + stat_summary(fun= fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(object = all.mesenchyme, assay = "RNA", features = "CRE", group.by = "merged.MES.clusters", 
        pt.size = 0, y.max = 4,add.noise = TRUE) + ggtitle("CRE") + geom_boxplot()

VlnPlot(object = all.mesenchyme, assay = "RNA", features = "Shox2", group.by = "merged.MES.clusters", 
        pt.size = 0, y.max = 4, slot = "data",add.noise = TRUE ) + ggtitle("Shox2")  + stat_summary(fun= fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(object = all.mesenchyme, assay = "RNA", features = "Shox2", group.by = "merged.MES.clusters", 
        pt.size = 0, y.max = 4,add.noise = TRUE) + ggtitle("Shox2") + geom_boxplot()

VlnPlot(object = all.mesenchyme, assay = "RNA", features = "EYFP", group.by = "merged.MES.clusters", 
        pt.size = 0, y.max = 4, slot = "data",add.noise = TRUE ) + ggtitle("EYFP")  + stat_summary(fun= fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(object = all.mesenchyme, assay = "RNA", features = "EYFP", group.by = "merged.MES.clusters", 
        pt.size = 0, y.max = 4,add.noise = TRUE) + ggtitle("EYFP") + geom_boxplot()

VlnPlot(object = all.mesenchyme, assay = "RNA", features = "Hoxd13", group.by = "merged.MES.clusters", 
        pt.size = 0, y.max = 4, slot = "data",add.noise = TRUE ) + ggtitle("Hoxd13")  + stat_summary(fun= fun1, geom="text", aes(label=sprintf("%1.3f", ..y..))) 

VlnPlot(object = all.mesenchyme, assay = "RNA", features = "Hoxd13", group.by = "merged.MES.clusters", 
        pt.size = 0, y.max = 4,add.noise = TRUE) + ggtitle("Hoxd13") + geom_boxplot()
#

VlnPlot(object = all.mesenchyme, assay = "RNA", features = "Smarca2", group.by = "merged.MES.clusters", 
        pt.size = 0, y.max = 4,add.noise = TRUE ) + ggtitle("Smarca2") 

VlnPlot(object = all.mesenchyme, assay = "RNA", features = "Smarca4", group.by = "merged.MES.clusters", 
        pt.size = 0, y.max = 4,add.noise = TRUE ) + ggtitle("Smarca4") 

VlnPlot(object = all.mesenchyme, assay = "RNA", features = "Smarca2", group.by = "stage", 
        pt.size = 0, y.max = 4,add.noise = TRUE ) + ggtitle("Smarca2") 

VlnPlot(object = all.mesenchyme, assay = "RNA", features = "Smarca4", group.by = "stage", 
        pt.size = 0, y.max = 4,add.noise = TRUE ) + ggtitle("Smarca4") 

VlnPlot(object = all.mesenchyme, assay = "RNA", features = "Shox2", group.by = "stage", 
        pt.size = 0, y.max = 4,add.noise = TRUE ) + ggtitle("Shox2") 

ggsave("Vlnplot_Shox2_only_mesenchyme_noiseT_split_by_stage.png", height = 5, width = 5)
ggsave("Vlnplot_Shox2_only_mesenchyme_noiseT_split_by_stage.eps", height = 5, width = 5)

VlnPlot(object = all.mesenchyme, assay = "RNA", features = "CRE", group.by = "stage", 
        pt.size = 0, y.max = 4,add.noise = FALSE ) + ggtitle("CRE") 

ggsave("Vlnplot_CRE_only_mesenchyme_noiseF_split_by_stage.png", height = 5, width = 5)
ggsave("Vlnplot_CRE_only_mesenchyme_noiseF_split_by_stage.eps", height = 5, width = 5)

VlnPlot(object = all.mesenchyme, assay = "RNA", features = "Shox2", group.by = "stage", 
        pt.size = 0, y.max = 4,add.noise = FALSE ) + ggtitle("Shox2") 

ggsave("Vlnplot_Shox2_only_mesenchyme_noiseF_split_by_stage.png", height = 5, width = 5)
ggsave("Vlnplot_Shox2_only_mesenchyme_noiseF_split_by_stage.eps", height = 5, width = 5)

VlnPlot(object = all.mesenchyme, assay = "RNA", features = "EYFP", group.by = "stage", 
        pt.size = 0, y.max = 4,add.noise = FALSE ) + ggtitle("EYFP") 

ggsave("Vlnplot_EYFP_only_mesenchyme_noiseF_split_by_stage.png", height = 5, width = 5)
ggsave("Vlnplot_EYFP_only_mesenchyme_noiseF_split_by_stage.eps", height = 5, width = 5)


VlnPlot(object = all.mesenchyme, assay = "RNA", features = "Col1a1", group.by = "merged.MES.clusters", 
        pt.size = 0, add.noise = TRUE) + ggtitle("Col1a1")


#Violins split by gene_class
Idents(all.mesenchyme) <- "gene_class"
levels(all.mesenchyme)

VlnPlot(object = all.mesenchyme, assay = "RNA", features = "Shox2", split.by = "stage", 
        pt.size = 0, y.max = 4,add.noise = TRUE ) + ggtitle("Shox2") 

VlnPlot(object = all.mesenchyme, assay = "RNA", features = "EYFP", split.by = "stage", 
        pt.size = 0, y.max = 4,add.noise = TRUE ) + ggtitle("EYFP") 


Idents(all.mesenchyme) <- "merged.MES.clusters"
levels(all.mesenchyme)

VlnPlot(object = all.mesenchyme, assay = "RNA", features = "Col1a1", split.by = "stage", 
        pt.size = 0,add.noise = TRUE ) + ggtitle("Col1a1") 

Idents(all.mesenchyme) <- "stage"
levels(all.mesenchyme)

VlnPlot(object = all.mesenchyme, assay = "RNA", features = "Col1a1", 
        pt.size = 0,add.noise = TRUE ) + ggtitle("Col1a1") 


#### U. Scvelo MERGED UMAP1.1 ####
all.mesenchyme <- readRDS(file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.mesenchyme_UMAP1.1MERGED.rds")
all.mesenchyme
summary(all.mesenchyme@meta.data)
DimPlot(all.mesenchyme, reduction = "umap", label = TRUE)

#extract hex color codes for a plot with three elements in ggplot2 
hex <- hue_pal()(15)

#overlay hex color codes on actual colors
show_col(hex)

colors <- c("#F8766D","#E58700","#C99800","#A3A500","#6BB100","#00BA38","#00BF7D","#00C0AF","#00BCD8","#00B0F6","#619CFF",
            "#B983FF","#E76BF3","#FD61D1","#FF67A4")


##ScVelo MERGED Mesenchyme UMAP1_1##
DefaultAssay(all.mesenchyme) <- "RNA"
all.mesenchyme

SaveH5Seurat(all.mesenchyme,filename = "Mesenchyme.MERGED.UMAP1_1.h5Seurat")
Convert("Mesenchyme.MERGED.UMAP1_1.h5Seurat", dest="h5ad")


#Create an environment in conda and install scvelo
conda create -n en_scvelo
conda activate en_scvelo
pip install -U scvelo


#How to run scvelo
conda activate en_scvelo
python
import scvelo as scv
adata = scv.read("Mesenchyme.MERGED.UMAP1_1.h5ad")
adata

ident_colours = ['#F8766D','#E58700','#C99800','#A3A500','#6BB100','#00BA38','#00BF7D','#00C0AF','#00BCD8','#00B0F6','#619CFF','#B983FF','#E76BF3','#FD61D1','#FF67A4']
adata.uns['merged.MES.clusters']=ident_colours

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=10, n_neighbors=30)

scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)

scv.pl.velocity_embedding_stream(adata, basis="umap", color="merged.MES.clusters", palette = ident_colours, legend_loc='right margin')


#### V. Alexandre_Lucille velocity analisis ####
all.mesenchyme <- readRDS(file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.mesenchyme_UMAP1.1MERGED.rds")
all.mesenchyme
summary(all.mesenchyme@meta.data)
DimPlot(all.mesenchyme, reduction = "umap", label = TRUE)

levels(all.mesenchyme)

#this is the general velocity analysis
#safelyLoadAPackageInCRANorBioconductor("velocyto.R")
bm <- all.mesenchyme
rm(all.mesenchyme)
Sys.time()
bm <- RunVelocity(object = bm, deltaT = 1, kCells = 25, fit.quantile = 0.02 )
Sys.time()
DimPlot(bm, label = T)

ident.colors <- (scales::hue_pal())(n = length(x = levels(x = bm)))
names(x = ident.colors) <- levels(x = bm)
cell.colors <- ident.colors[Idents(object = bm)]
names(x = cell.colors) <- colnames(x = bm)
#velocity <- list()
#velocity[["stage"]] 
velocity <- show.velocity.on.embedding.cor(emb = Embeddings(object = bm, 
                                                            reduction = "umap"), 
                                           vel = Tool(object = bm, 
                                                      slot = "RunVelocity"), 
                                           n = 200, 
                                           scale = "sqrt", 
                                           cell.colors = ac(x = cell.colors, 
                                                            alpha = 0.5), 
                                           cex = 0.8, 
                                           arrow.scale = 3, 
                                           show.grid.flow = TRUE, 
                                           min.grid.cell.mass = 0.5, 
                                           grid.n = 40, 
                                           arrow.lwd = 1, 
                                           do.par = FALSE, 
                                           cell.border.alpha = 0.1, 
                                           min.arrow.size = 0.0, 
                                           n.cores = 5, return.details = T)
saveRDS(velocity, file = "velocity.rds")

#the gene specific dynamic is here
#bm.list <- list(bm)
velocity <- readRDS("/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/velocity.rds")

ggplot(as.data.frame(Embeddings(object = bm,reduction = "umap")), aes(x = UMAP_1, y = UMAP_2)) + 
  geom_point(aes(colour = velocity[["vel"]]["Shox2",])) +
  scale_colour_gradientn(colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu"))


DefaultAssay(bm) <- "RNA"
bm.list <- list(bm)
gene.wanted <- "Shox2"
ggplot(as.data.frame(Embeddings(object = bm.list[[1]],reduction = "umap")),
       aes(x = UMAP_1, y = UMAP_2)) + 
  geom_point(aes(colour = velocity[["vel"]][gene.wanted,])) +
  scale_colour_gradientn(colors = RColorBrewer::brewer.pal(n = 11, name = "RdBu"))+
  theme_classic() +
  NoLegend() +
  ggtitle(gene.wanted, subtitle = "Transcript dynamic")

#To check the Aggregate counts in each Assay 
AggregateExpression(all.mesenchyme, features  = "Shox2", slot = "counts")
AggregateExpression(all.mesenchyme, features  = "Hoxd13", slot = "counts")
AggregateExpression(all.mesenchyme, features  = "EYFP", slot = "counts")
AggregateExpression(all.mesenchyme, features  = "CRE", slot = "counts")
AggregateExpression(all.mesenchyme, features  = "CRE", slot = "counts")


AggregateExpression(all.mesenchyme, features  = c("Shox2", "Hoxd13"), slot = "counts")

AggregateExpression(bm, features  = "EYFP", slot = "counts")

#### W. monocle ####
all.mesenchyme <- readRDS(file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.mesenchyme_UMAP1.1MERGED.rds")
all.mesenchyme
summary(all.mesenchyme@meta.data)
Idents(all.mesenchyme)<-"merged.MES.clusters"
DimPlot(all.mesenchyme, reduction = "umap", label = TRUE)
DefaultAssay(all.mesenchyme) <- "RNA"
all.mesenchyme

#monocle3 uses a cell_data_set object, 
#the as.cell_data_set function from SeuratWrappers can be used to “convert” a Seurat object to Monocle object
all.mesenchyme.mono <- as.cell_data_set(all.mesenchyme)
all.mesenchyme.mono <- cluster_cells(all.mesenchyme.mono, cluster_method = 'louvain')

plot_cells(all.mesenchyme.mono, color_cells_by = "cluster", show_trajectory_graph = FALSE)
plot_cells(all.mesenchyme.mono, color_cells_by = "partition", show_trajectory_graph = FALSE)

#The first step in trajectory analysis is the learn_graph() function
integrated.sub <- subset(as.Seurat(all.mesenchyme.mono, assay = NULL), monocle3_partitions == 1)
all.mesenchyme.mono <- as.cell_data_set(integrated.sub)
all.mesenchyme.mono <- learn_graph(all.mesenchyme.mono, use_partition = TRUE, verbose = FALSE)

#After learning the graph, monocle can plot add the trajectory graph to the cell plot.
plot_cells(all.mesenchyme.mono,
           color_cells_by = "cluster",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)


#We can set the root to any one of our clusters by selecting the cells in that cluster to use as the root in the function order_cells.
all.mesenchyme.mono <- order_cells(all.mesenchyme.mono, root_cells = colnames(all.mesenchyme.mono[,clusters(all.mesenchyme.mono) == 2]))
plot_cells(all.mesenchyme.mono,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           )

#option 2
#max.avp <- which.max(unlist(FetchData(all.mesenchyme, "Shox2")))
#max.avp <- colnames(all.mesenchyme)[max.avp]
#all.mesenchyme.mono <- order_cells(all.mesenchyme.mono, root_cells = max.avp)
#plot_cells(all.mesenchyme.mono, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
#           label_branch_points = FALSE)

#We can export this data to the Seurat object and visualize
integrated.sub <- as.Seurat(all.mesenchyme.mono, assay = NULL)
summary(integrated.sub@meta.data)
FeaturePlot(integrated.sub, "monocle3_pseudotime")

#Monocle’s graph_test() function detects genes that vary over a trajectory.
all.mesenchyme.mono_graph_test_results <- graph_test(all.mesenchyme.mono,
                                     neighbor_graph = "principal_graph",
                                     cores = 8)

rowData(all.mesenchyme.mono)$gene_short_name <- row.names(rowData(all.mesenchyme.mono))

head(all.mesenchyme.mono_graph_test_results, error=FALSE, message=FALSE, warning=FALSE)

deg_ids <- rownames(subset(all.mesenchyme.mono_graph_test_results[order(all.mesenchyme.mono_graph_test_results$morans_I, decreasing = TRUE),], q_value < 0.05))

plot_cells(all.mesenchyme.mono,
           genes=head(deg_ids),
           show_trajectory_graph = FALSE,
           label_cell_groups = FALSE,
           label_leaves = FALSE)

gene_modules <- find_gene_modules(all.mesenchyme.mono[deg_ids,],
                                  resolution=c(10^seq(-6,-1)))
table(gene_modules$module)

#### X.1 Renaming gene_class creating final class ####
#all clusters including sattelites
all.merged.SCTCCSt.labelled <- readRDS("/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.merged.SCTCCSt.labelled.rds")
all.merged.SCTCCSt.labelled
summary(all.merged.SCTCCSt.labelled@meta.data)
DimPlot(all.merged.SCTCCSt.labelled, reduction = "umap",label = TRUE, label.size = 3)

Idents(all.merged.SCTCCSt.labelled) <- "gene_class"
Idents(all.merged.SCTCCSt.labelled)
levels(all.merged.SCTCCSt.labelled)
DimPlot(all.merged.SCTCCSt.labelled, reduction = "umap", label = FALSE)

all.merged.SCTCCSt.labelled@meta.data$finalclass <- all.merged.SCTCCSt.labelled@meta.data$gene_class
all.merged.SCTCCSt.labelled@meta.data$finalclass

all.merged.SCTCCSt.labelled@meta.data$finalclass[all.merged.SCTCCSt.labelled@meta.data$finalclass == "EYFPpos_Shox2pos_CREpos_cells"] <- "Maintained"
all.merged.SCTCCSt.labelled@meta.data$finalclass[all.merged.SCTCCSt.labelled@meta.data$finalclass == "EYFPpos_Shox2pos_CREneg_cells"] <- "Maintained"
all.merged.SCTCCSt.labelled@meta.data$finalclass[all.merged.SCTCCSt.labelled@meta.data$finalclass == "EYFPneg_Shox2pos_CREneg_cells"] <- "Initiating"
all.merged.SCTCCSt.labelled@meta.data$finalclass[all.merged.SCTCCSt.labelled@meta.data$finalclass == "EYFPneg_Shox2neg_CREneg_cells"] <- "Negative"
all.merged.SCTCCSt.labelled@meta.data$finalclass[all.merged.SCTCCSt.labelled@meta.data$finalclass == "EYFPneg_Shox2pos_CREpos_cells"] <- "Initiating"
all.merged.SCTCCSt.labelled@meta.data$finalclass[all.merged.SCTCCSt.labelled@meta.data$finalclass == "EYFPneg_Shox2neg_CREpos_cells"] <- "Other_cells"
all.merged.SCTCCSt.labelled@meta.data$finalclass[all.merged.SCTCCSt.labelled@meta.data$finalclass == "EYFPpos_Shox2neg_CREneg_cells"] <- "Decommissioned"
all.merged.SCTCCSt.labelled@meta.data$finalclass[all.merged.SCTCCSt.labelled@meta.data$finalclass == "EYFPpos_Shox2neg_CREpos_cells"] <- "Other_cells"

all.merged.SCTCCSt.labelled@meta.data$finalclass
Idents(all.merged.SCTCCSt.labelled) <- "finalclass"
Idents(all.merged.SCTCCSt.labelled)

DimPlot(all.merged.SCTCCSt.labelled, reduction = "umap", label = FALSE)
DimPlot(all.merged.SCTCCSt.labelled, reduction = "umap", label = FALSE, cols = c('Decommissioned'='#f39200','Maintained'='#be1622',
                                                                                 'Initiating'='#2d2e83','Negative' ='#000000',
                                                                                 'Other_cells'='#e5e5e5'), pt.size = 0.5)

saveRDS(all.merged.SCTCCSt.labelled, "/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.merged.SCTCCSt.labelled2.rds")

#only mesenchyme
all.mesenchyme <- readRDS(file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.mesenchyme_UMAP1.1MERGED.rds")
all.mesenchyme
summary(all.mesenchyme@meta.data)
DimPlot(all.mesenchyme, reduction = "umap", label = TRUE)

Idents(all.mesenchyme) <- "gene_class"
Idents(all.mesenchyme)
levels(all.mesenchyme)
DimPlot(all.mesenchyme, reduction = "umap", label = TRUE)

all.mesenchyme@meta.data$finalclass <- all.mesenchyme@meta.data$gene_class
all.mesenchyme@meta.data$finalclass

all.mesenchyme@meta.data$finalclass[all.mesenchyme@meta.data$finalclass == "EYFPpos_Shox2pos_CREpos_cells"] <- "Maintained"
all.mesenchyme@meta.data$finalclass[all.mesenchyme@meta.data$finalclass == "EYFPpos_Shox2pos_CREneg_cells"] <- "Maintained"
all.mesenchyme@meta.data$finalclass[all.mesenchyme@meta.data$finalclass == "EYFPneg_Shox2pos_CREneg_cells"] <- "Initiating"
all.mesenchyme@meta.data$finalclass[all.mesenchyme@meta.data$finalclass == "EYFPneg_Shox2neg_CREneg_cells"] <- "Negative"
all.mesenchyme@meta.data$finalclass[all.mesenchyme@meta.data$finalclass == "EYFPneg_Shox2pos_CREpos_cells"] <- "Initiating"
all.mesenchyme@meta.data$finalclass[all.mesenchyme@meta.data$finalclass == "EYFPneg_Shox2neg_CREpos_cells"] <- "Other_cells"
all.mesenchyme@meta.data$finalclass[all.mesenchyme@meta.data$finalclass == "EYFPpos_Shox2neg_CREneg_cells"] <- "Decommissioned"
all.mesenchyme@meta.data$finalclass[all.mesenchyme@meta.data$finalclass == "EYFPpos_Shox2neg_CREpos_cells"] <- "Other_cells"

all.mesenchyme@meta.data$finalclass
Idents(all.mesenchyme) <- "finalclass"
Idents(all.mesenchyme)

DimPlot(all.mesenchyme, reduction = "umap", label = FALSE)
DimPlot(all.mesenchyme, reduction = "umap", label = FALSE, cols = c('Decommissioned'='#f39200','Maintained'='#be1622',
                                                                                 'Initiating'='#2d2e83','Negative' ='#000000',
                                                                                 'Other_cells'='#e5e5e5'), pt.size = 0.5)

saveRDS(all.mesenchyme, "/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.mesenchyme_UMAP1.1MERGED2.rds")

#### X.2 Final class cell proportion ####
#all clusters inlcuding satellites
all.merged.SCTCCSt.labelled <- readRDS("/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.merged.SCTCCSt.labelled2.rds")
all.merged.SCTCCSt.labelled
summary(all.merged.SCTCCSt.labelled@meta.data)
DimPlot(all.merged.SCTCCSt.labelled, reduction = "umap",label = TRUE, label.size = 3)
DimPlot(all.merged.SCTCCSt.labelled, reduction = "umap", label = FALSE, cols = c('Decommissioned'='#F1D31E','Maintained'='#be1622',
                                                                                 'Initiating'='#6EC4EE','Negative' ='#000000',
                                                                                 'Other_cells'='#e5e5e5'), pt.size = 0.5)


propotions_cells_percluster4 <- prop.table(table(all.merged.SCTCCSt.labelled$finalclass, all.merged.SCTCCSt.labelled$stage), margin = 2)
propotions_cells_percluster4.df <- as.data.frame(propotions_cells_percluster4)
t(propotions_cells_percluster4)

ggplot(propotions_cells_percluster4.df, aes(fill= Var1, x=Var2,y=Freq)) + 
  geom_bar(position="fill", stat = "identity") + 
  geom_text(aes(label = sprintf("%0.2f", round(Freq, 2))), color="white", size=3,position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c('Decommissioned'='#F1D31E','Maintained'='#be1622',
                               'Initiating'='#6EC4EE','Negative' ='#000000',
                               'Other_cells'='#e5e5e5')) +theme_light()

write.csv(propotions_cells_percluster4.df,"UMAP_MERGED_all_clusters_cell_proportions_perstage_finalclass.csv")


propotions_cells_percluster5 <- prop.table(table(all.merged.SCTCCSt.labelled$finalclass, all.merged.SCTCCSt.labelled$merged.clusters), margin = 2)
propotions_cells_percluster5.df <- as.data.frame(propotions_cells_percluster5)
t(propotions_cells_percluster5)

write.csv(propotions_cells_percluster5.df,"UMAP_MERGED_all_clusters_cell_proportions_percluster_finalclass.csv")

ggplot(propotions_cells_percluster5.df, aes(fill= Var1, x=Var2,y=Freq)) + 
  geom_bar(position="fill", stat = "identity", col="black") + 
  geom_text(aes(label = sprintf("%0.2f", round(Freq, 2))), color="white", size=3,position = position_stack(vjust = 0.5))+
  scale_fill_manual(values = c('Decommissioned'='#F1D31E','Maintained'='#be1622',
                               'Initiating'='#6EC4EE','Negative' ='#000000',
                               'Other_cells'='#e5e5e5')) + theme_classic()

## We finally choose to use the Noise option = TRUE
DefaultAssay(all.merged.SCTCCSt.labelled) <- "RNA"
all.merged.SCTCCSt.labelled

#for the violin plot I dediced to plot the average expression log1p normalised since the values plot a
VlnPlot(object = all.merged.SCTCCSt.labelled, assay = "RNA", features = "CRE", group.by = "finalclass", 
        pt.size = 0, y.max = 4, slot = "data",add.noise = FALSE ) + ggtitle("CRE") 

VlnPlot(object = all.merged.SCTCCSt.labelled, assay = "RNA", features = "Shox2", group.by = "finalclass", 
        pt.size = 0, y.max = 4, slot = "data",add.noise = FALSE ) + ggtitle("Shox2")  


#only mesenchyme
all.mesenchyme <- readRDS(file="/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.mesenchyme_UMAP1.1MERGED2.rds")
all.mesenchyme
summary(all.mesenchyme@meta.data)
DimPlot(all.mesenchyme, reduction = "umap", label = FALSE)

DimPlot(all.mesenchyme, reduction = "umap", label = FALSE, cols = c('Decommissioned'='#F1D31E','Maintained'='#be1622',
                                                                    'Initiating'='#6EC4EE','Negative' ='#000000',
                                                                    'Other_cells'='#e5e5e5'), pt.size = 0.5)

propotions_cells_percluster4 <- prop.table(table(all.mesenchyme$finalclass, all.mesenchyme$merged.MES.clusters), margin = 2)
propotions_cells_percluster4.df <- as.data.frame(propotions_cells_percluster4)
t(propotions_cells_percluster4)

write.csv(propotions_cells_percluster4.df,"UMAP_Mesenchyme_all_cell_proportions_percluster_finalclass.csv")

positions <- c("LDC","EDC","IM","DP","DPP","Ms","LP", "EPP", "LPP", "PCT","ICT","TP", "EPC", "PGP", "PC")
propotions_cells_percluster4.df$Var2 <- factor(propotions_cells_percluster4.df$Var2, levels = positions)

ggplot(propotions_cells_percluster4.df, aes(fill= Var1, x=Var2,y=Freq)) + 
  geom_bar(position="fill", stat = "identity", col="black") + 
  geom_text(aes(label = sprintf("%0.2f", round(Freq, 2))), color="white", size=3,position = position_stack(vjust = 0.5))+
  scale_fill_manual(values = c('Decommissioned'='#F1D31E','Maintained'='#be1622',
                               'Initiating'='#6EC4EE','Negative' ='#000000',
                               'Other_cells'='#e5e5e5')) + theme_classic()

#### X.3 Markers by stage of final class groups in all cells including satellites #####
#all clusters inlcuding satellites
all.merged.SCTCCSt.labelled <- readRDS("/Volumes/RAQUEL-LAB3/Single_cell_58_data/8_single_cell_clone_58/scclone58_analisis_all_stages_with_velocyto/sc_58_all_stages_with_velocyto_regress_stage/all.merged.SCTCCSt.labelled2.rds")
all.merged.SCTCCSt.labelled
summary(all.merged.SCTCCSt.labelled@meta.data)
DimPlot(all.merged.SCTCCSt.labelled, reduction = "umap",label = TRUE, label.size = 3)
DimPlot(all.merged.SCTCCSt.labelled, reduction = "umap", label = FALSE, cols = c('Decommissioned'='#F1D31E','Maintained'='#be1622',
                                                                                 'Initiating'='#6EC4EE','Negative' ='#000000',
                                                                                 'Other_cells'='#e5e5e5'), pt.size = 0.5)

DefaultAssay(all.merged.SCTCCSt.labelled) <- "RNA"
all.merged.SCTCCSt.labelled
levels(all.merged.SCTCCSt.labelled)


all.merged.SCTCCSt.labelled@meta.data$finalclass_stage <- paste(all.merged.SCTCCSt.labelled@meta.data$finalclass, all.merged.SCTCCSt.labelled$stage, sep = "_")

Idents(all.merged.SCTCCSt.labelled) <- "finalclass_stage"
levels(all.merged.SCTCCSt.labelled)

Maintained_E105.markers <-  FindMarkers(all.merged.SCTCCSt.labelled, ident.1 = "Maintained_E105",
                            ident.2 = NULL,
                            only.pos =TRUE,
                            logfc.threshold = 0.7,
                            pseudocount.use = 0,
                            min.diff.pct=0.1)

write.csv(Maintained_E105.markers,file=("Maintained_E105.markers.csv"))

Maintained_E135.markers <-  FindMarkers(all.merged.SCTCCSt.labelled, ident.1 = "Maintained_E135",
                                        ident.2 = NULL,
                                        only.pos =TRUE,
                                        logfc.threshold = 0.7,
                                        pseudocount.use = 0,
                                        min.diff.pct=0.1)

write.csv(Maintained_E135.markers,file=("Maintained_E135.markers.csv"))

Decommissioned_E105.markers <-  FindMarkers(all.merged.SCTCCSt.labelled, ident.1 = "Decommissioned_E105",
                                        ident.2 = NULL,
                                        only.pos =TRUE,
                                        logfc.threshold = 0.7,
                                        pseudocount.use = 0,
                                        min.diff.pct=0.1)

write.csv(Decommissioned_E105.markers,file=("Decommissioned_E105.markers.csv"))

Decommissioned_E135.markers <-  FindMarkers(all.merged.SCTCCSt.labelled, ident.1 = "Decommissioned_E135",
                                        ident.2 = NULL,
                                        only.pos =TRUE,
                                        logfc.threshold = 0.7,
                                        pseudocount.use = 0,
                                        min.diff.pct=0.1)

write.csv(Decommissioned_E135.markers,file=("Decommissioned_E135.markers.csv"))


Negative_E105.markers <-  FindMarkers(all.merged.SCTCCSt.labelled, ident.1 = "Negative_E105",
                                            ident.2 = NULL,
                                            only.pos =TRUE,
                                            logfc.threshold = 0.7,
                                            pseudocount.use = 0,
                                            min.diff.pct=0.1)

write.csv(Negative_E105.markers,file=("Negative_E105.markers.csv"))

Negative_E135.markers <-  FindMarkers(all.merged.SCTCCSt.labelled, ident.1 = "Negative_E135",
                                            ident.2 = NULL,
                                            only.pos =TRUE,
                                            logfc.threshold = 0.7,
                                            pseudocount.use = 0,
                                            min.diff.pct=0.1)

write.csv(Negative_E135.markers,file=("Negative_E135.markers.csv"))






