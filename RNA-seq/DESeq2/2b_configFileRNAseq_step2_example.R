### Required for all steps ###
RNAseqFunctionPath<-"/Users/roucogar/Library/CloudStorage/Dropbox/Git_Hub_Lucille/rnaseq_rscripts-master/RNAseqFunctions.R"
samplesPlan<-"/Users/roucogar/Library/CloudStorage/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/7_RNA-seq/Final_analysis/RNA-seq_58_Forelimb/samplesPlan_example.txt" #This file should be a tabulated file with at least one column called "sample". Optionnaly, the paths to the counts tables and FPKM tables can be provided under the column called: htseq_count_file and cufflinks_file.


#### STEP 2 - DESEQ 2 ANALYSIS ###
#Required
tableWithCounts<-"/Users/roucogar/Library/CloudStorage/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/7_RNA-seq/Final_analysis/RNA-seq_58_Forelimb/outputs/mergedTables/AllHTSeqCounts_subset.txt"
geneIDColCounts<-"Ens_ID" #Specify here the name of the column which contains the gene IDs (they need to be unique).
#For the DESeq2 analysis you need to specify a factor on which you want to do the analysis:
factor<-"Stage" #This needs to be a name of a column of the samplesPlan file.

#Optional
tableWithAnnotations<-"/Users/roucogar/Library/CloudStorage/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/7_RNA-seq/Final_analysis/RNA-seq_58_Forelimb/outputs/mergedTables/AllCufflinks_Simplified_norm.txt" #This can be table from cufflinks or cuffdiff or Biomart to annotate genes. You will need to choose a file with at least one column with the Ensembl Gene IDs.
geneIDColInAnnotations<-"gene_id"  #Specify here the name of the column which contains the gene IDs (it must match with the content of the geneID from the table with counts).
#You can also provide a gtf:
gtfFile<-"/Users/roucogar/Library/CloudStorage/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/7_RNA-seq/Final_analysis/RNA-seq_58_Forelimb/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm39.104_ExonsOnly_UCSC_dsmCherry_P2A_CRE_EYFP.gtf"
changeTest<-F #Default test is Wald but you can change to likelihood ratio test (LRT) with reduced formula ~1. Put F to keep Wald and put T to use LRT.
outputDESeqTable<-"/Users/roucogar/Library/CloudStorage/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/7_RNA-seq/Final_analysis/RNA-seq_58_Forelimb/outputs/DESeq2_DP58_E105vsE115/DESeq2Analysis_PerGeno.txt"
outputSignificantTable<-T #If you want to have another table with only significant genes abs(l2FC)>1.5 and corrected p-value<0.05
