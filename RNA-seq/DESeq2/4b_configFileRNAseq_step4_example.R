### Required for all steps ###
RNAseqFunctionPath<-"/Users/roucogar/Library/CloudStorage/Dropbox/Git_Hub_Lucille/rnaseq_rscripts-master/RNAseqFunctions.R"

#### STEP 4 - PLOTS RESULTS OF DESEQ2 ###
#Required
tableWithResultsOfDifferentialExpression<-"/Users/roucogar/Library/CloudStorage/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/7_RNA-seq/Final_analysis/RNA-seq_58_Forelimb/outputs/DESeq2_DP58_E105vsE115/DESeq2Analysis_PerGeno.txt" #You can put here the result of step2 or any table with at least 2 columns (padj and log2FoldChange) for Volcano and 3 columns (padj, log2FoldChange and baseMean) for MAP.

outputFolder<-"/Users/roucogar/Library/CloudStorage/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/7_RNA-seq/Final_analysis/RNA-seq_58_Forelimb/outputs/DESeq2_DP58_E105vsE115/"
usePng<-T#By default pdf is used as output. If set to T, png will be used.

### Common for Volcano and MAP ###
#significant thresholds:
maxPAdj<-0.05 #Put the maximum value of ajusted p-value to be significant (default 0.05)
minLFC<-1.5 #Put the minimum log2 fold change to be significant (default 1.5)
#colors
colOfNonSignificant<-"grey"
colOfSignificant<-"blue"
click<-T #Do you want to click on the plot to be able to identify some genes. (T=yes, F=no)
geneID<-"gene_name" #If you want to click, provide here the name of the column which can be used to label the gene.
fileWithGenes<-"/Users/roucogar/Library/CloudStorage/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/7_RNA-seq/Final_analysis/RNA-seq_58_Forelimb/relevant_genes.txt" #You can also provide a list of genes you want to identify on the plots. One gene per line. The first line of the gene file should correspond to a column in the tableWithResultsOfDifferentialExpression file.
colOfCircle<-"red"

### Volcano ###

maxYVolcano<-NA #If you want to zoom on bigger p-values because you have some genes with very low p-values, you can put here a value to restrict the plot to p-values higher than 10^-maxYVolcano. Put NA or comment line if you do not need to zoom the initial plot.

### MAP ###

ylimMAP<-NA #If you want to zoom on smaller log2 fold changes because you have some genes with very high log2 fold changes, you can put here minimum and maxiumum log2 fold change values to restrict the plot to these values. Put NA or comment line if you do not need to zoom the initial plot, put c(minValue,maxValue) if you want to restrict.
