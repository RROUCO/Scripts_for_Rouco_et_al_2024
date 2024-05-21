### Required for all steps ###
RNAseqFunctionPath<-"/Users/roucogar/Library/CloudStorage/Dropbox/Git_Hub_Lucille/rnaseq_rscripts-master/RNAseqFunctions.R"
samplesPlan<-"/Users/roucogar/Library/CloudStorage/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/7_RNA-seq/Final_analysis/RNA-seq_58_Forelimb/samplesPlan_example.txt" #This file should be a tabulated file with at least one column called "sample". Optionnaly, the paths to the counts tables and FPKM tables can be provided under the column called: htseq_count_file and cufflinks_file.


#### STEP 3 - PLOTS ###
#Required
  tableWithNormalizedExpression<-"/Users/roucogar/Library/CloudStorage/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/7_RNA-seq/Final_analysis/RNA-seq_58_Forelimb/outputs/DESeq2_DP58_E105vsE115/DESeq2Analysis_PerGeno.txt" #You can put here either the FPKM norm values or the count norm values obtained after DESeq2
useFPKM<-T #In case you are using a file with both raw counts and FPKM you need to choose which values you want to plot. If set to T, only columns called FPKM_sample will be used.
#Optional
outputFolder<-"/Users/roucogar/Library/CloudStorage/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/7_RNA-seq/Final_analysis/RNA-seq_58_Forelimb/outputs/DESeq2_DP58_E105vsE115/"
usePng<-T#By default pdf is used as output. If set to T, png will be used.
fixedColors<-list(Genotype=c('E105'="lightgreen",'E115'="red4")) #You can provide color for each value of each factor in you samples plan to have constistant graphs for PCA, correlation and genes.

### Common to PCA and clustering ###
restrictToNMoreVariantGenes<-500 #In DESeq2 they restrict to the 500 most variant genes. If you want to keep all genes, comment the line or put 1000000.

### PCA ###
nbOfPC<-3 #Put here the number of PC you want to see (0=do not perform PCA, 1=Only look at first component, 2=look at the 2 first etc...)
#If the nbOfPC is greater than 1, you will have a barchart of each PC and you may want to use different parameters to identify your samples using the column names of the samples plan.
PCA1D<-list(color="Stage", fill="Replicate", shape="Stage") 
#Possible personalizations are :
#fill is for the color of the bar, 
#alpha is for the transparency,
#color is for the color of the border of the bar
#linetype is for the type of border of the bar
#If you do not want to use one of the parameter, just remove it from the list.

#If the nbOfPC is greater than 2, you will have projection in 2 dimension and to identify your sample you may want to use the column names of the samples plan.
PCA2D<-list(color="Stage", fill="Replicate", shape="Stage")
#Possible personalizations are :
#color is for the color of the symbol, 
#alpha is for the transparency,
#shape is for the shape of the symbol.
#You can also choose 2 colors fill and color. If so, the fill will be inside and the color will be the border.
#If you do not want to use one of the parameter, just remove it from the list.

getGeneContributionToPCA<-T #Do you want to have the contribution of each gene to each PC (T=yes, F=no).

### Clustering ###

plotMatrixAndClustering<-T#Do you want to perform a correlation matrix and clustering (T=yes, F=no)

### Genes ###
fileWithGenes<-"/Users/roucogar/Library/CloudStorage/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/7_RNA-seq/Final_analysis/RNA-seq_58_Forelimb/relevant_genes.txt" #One gene per line. The first line of the gene file should correspond to a column in the expression file.
# geneIDToAdd<-"gene_short_name" #By default, the title of the plot is the id provided in the fileWithGenes but you can add a meaning full name like gene_short_name if it is provided in the tableWithNormalizedExpression.
useLogExpression<-T #By default, the values of expression plotted are log2(1+expression) (when T, if F the raw expression will be plotted.)
useSameYmaxForAllGenes<-T #By default, each gene is plotted on an adjusted scale. If useSameYmaxForAllGenes is T, all genes will be plotted with the same y axis.
xaxisForGenes<-"Genotype" #A factor which will be used as x axis.
plotGenesPara<-list(color="Stage",fill="Replicate", shape="Stage")
#Possible personalizations are :
#color is for the color of the symbol, 
#alpha is for the transparency,
#shape is for the shape of the symbol.
#You can also choose 2 colors fill and color. If so, the fill will be inside and the color will be the border.
#If you do not want to use one of the parameter, just remove it from the list.
doNotPlotGeneByGene<-F #If you only want a heatmap and not one gene per one gene. Put it to T.
addGlobalHeatmap<-T #If you want to have a heatmap with the genes provided in the list. All values of plotGenesPara will be used to annotate the samples.
keepGeneOrder<-T #By default, genes are clustered by euclidean distance and complete clustering. If you want to keep the original order. Put keepGeneOrder to T.
clusterSamples<-F #By default, samples are not clustered.
