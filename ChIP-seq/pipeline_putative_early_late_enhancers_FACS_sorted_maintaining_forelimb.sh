######Run broadpeak MACS2 peak calling on the cluster
#https://github.com/macs3-project/MACS/blob/master/docs/callpeak.md

# --nomodel + --extsize 200 Whether or not to build the shifting model. If True, MACS will not build model. While --nomodel is set, MACS uses --extsize parameter to extend reads in 5'->3' direction to fix-sized fragments. For example, if the size of the binding region for your transcription factor is 200 bp, and you want to bypass the model building by MACS, this parameter can be set as 200.

# -B 2 Whether or not to save extended fragment pileup, and local lambda tracks (two files) at every bp into a bedGraph file

# --broad If set, MACS will try to call broad peaks using the --broad-cutoff setting. The maximum gap is expanded to 4 * MAXGAP (--max-gap parameter)As a result, MACS will output a 'gappedPeak' and a 'broadPeak' file instead of 'narrowPeak' file.

# --broad-cutoff 0.05 Cutoff for broad region. This option is not available unless --broad is set. If -p is set, this is a pvalue cutoff, otherwise,it's a qvalue cutoff. 

# --gsize mm parameter set to mm to fixed to the genome size that can be sequenced

#  --nolambda With this flag on, MACS will use the background lambda as local lambda. This means MACS will not consider the local bias at peak candidate regions.The whole point of the local lambda is to estimate the local background in the very situation that no input file is available. --nolambda turns this off and is imho a terrible choice for any peak calling situation as it then assumes uniform background across the genome. You will get plenty of false-positives with low enrichments.

# --call-summits it's incompatible with broad peak calling


# Explanation to remove the lambda: "Second, the estimation of dynamic background works well for histone modification ChIP-Seq data with control. However, if no control data is available for ChIP-Seq data with broad peaks, the local background estimation should be skipped via setting --nolambda in the command line.""
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3120977/
#trial3
macs2 callpeak -t ${sample}_mapped_sorted_q30.bam -n ${sample}_macs_SR200 --broad --nolambda --broad-cutoff 0.05 --nomodel --gsize mm --extsize 200 -B 2> ${pathResults}${sample}_macs_SR200.log


####### Run local
# build beds from broadPeak
cut -f 1-3 /Users/roucogar/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/6_Chip-seq/broad_trial3/H3K27ac-FL-E105-GA58_DP_Rep1_macs_SR200_peaks.broadPeak | uniq > H3K27ac-FL-E105-GA58_DP_Rep1_macs_SR200_broadpeaks.bed
cut -f 1-3 /Users/roucogar/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/6_Chip-seq/broad_trial3/H3K27ac-FL-E115-GA58_DP_Rep1_macs_SR200_peaks.broadPeak | uniq > H3K27ac-FL-E115-GA58_DP_Rep1_macs_SR200_broadpeaks.bed
cut -f 1-3 /Users/roucogar/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/6_Chip-seq/broad_trial3/H3K27ac-FL-E125-GA58_DP_Rep1_macs_SR200_peaks.broadPeak | uniq > H3K27ac-FL-E125-GA58_DP_Rep1_macs_SR200_broadpeaks.bed
cut -f 1-3 /Users/roucogar/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/6_Chip-seq/broad_trial3/H3K27ac-FL-E135-GA58_DP_Rep1_macs_SR200_peaks.broadPeak | uniq > H3K27ac-FL-E135-GA58_DP_Rep1_macs_SR200_broadpeaks.bed

# merge peaks that are closer than 600bp 
bedtools merge -header -i H3K27ac-FL-E105-GA58_DP_Rep1_macs_SR200_broadpeaks.bed -d 600 > H3K27ac-FL-E105-GA58_DP_macs_SR200_broadpeaks_600bp_merged.bed
bedtools merge -header -i H3K27ac-FL-E115-GA58_DP_Rep1_macs_SR200_broadpeaks.bed -d 600 > H3K27ac-FL-E115_GA58_DP_macs_SR200_broadpeaks_600bp_merged.bed
bedtools merge -header -i H3K27ac-FL-E125-GA58_DP_Rep1_macs_SR200_broadpeaks.bed -d 600 > H3K27ac-FL-E125_GA58_DP_macs_SR200_broadpeaks_600bp_merged.bed
bedtools merge -header -i H3K27ac-FL-E135-GA58_DP_Rep1_macs_SR200_broadpeaks.bed -d 600 > H3K27ac-FL-E135_GA58_DP_macs_SR200_broadpeaks_600bp_merged.bed

#using bedops to create a list of common broad H3K27Ac
#The --merge operation flattens all disjoint, overlapping, and adjoining element regions into contiguous, disjoint regions:
bedops --merge H3K27ac-FL-E105-GA58_DP_macs_SR200_broadpeaks_600bp_merged.bed H3K27ac-FL-E115_GA58_DP_macs_SR200_broadpeaks_600bp_merged.bed > H3K27ac-FL-E105_E115_GA58_DP_macs_SR200_broadpeaks_merged.bed
bedops --merge H3K27ac-FL-E125_GA58_DP_macs_SR200_broadpeaks_600bp_merged.bed H3K27ac-FL-E135_GA58_DP_macs_SR200_broadpeaks_600bp_merged.bed > H3K27ac-FL-E125_E135_GA58_DP_macs_SR200_broadpeaks_merged.bed
bedops --merge H3K27ac-FL-E105_E115_GA58_DP_macs_SR200_broadpeaks_merged.bed H3K27ac-FL-E125_E135_GA58_DP_macs_SR200_broadpeaks_merged.bed > H3K27ac-FL_GA58_DP_macs_SR200_broadpeaks_merged.bed

#Extend peaks up/down stream 300bp
bedtools slop -i H3K27ac-FL_GA58_DP_macs_SR200_broadpeaks_merged.bed -g /Users/roucogar/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/Annotations_Genomes/mm39_dsmCherry_P2A_CRE_EYFP/mm39_dsmCherry_P2A_CRE_EYFP.fa.fai -b 300 > H3K27ac-FL_GA58_DP_macs_SR200_broadpeaks_merged_extended.bed

#select only the TAD region coordinates in mm39 chr3:65,885,132-67,539,263
bedtools intersect -a /Users/roucogar/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/6_Chip-seq/Shox2_TAD_Coordinates.bed -b H3K27ac-FL_GA58_DP_macs_SR200_broadpeaks_merged_extended.bed > H3K27ac-FL_GA58_DP_macs_SR200_broadpeaks_merged_extended_Shox2_TAD.bed
bedtools intersect -a Shox2_TAD_Coordinates.bed -b H3K27ac-FL_GA58_DP_macs_SR200_broadpeaks_merged.bed > H3K27ac-FL_GA58_DP_macs_SR200_broadpeaks_merged_Shox2_TAD.bed

#conda install -c bioconda deeptools

deeptools multiBigwigSummary 
#multiBigwigSummary computes the average scores for each of the files in every genomic region.
#This analysis is performed for the entire genome by running the program in bins mode, or for certain user selected regions in BED-file mode.

BED-file --BED 
#the user provides a BED file that contains all regions that should be considered for the analysis. 

--smartLabels
# Instead of manually specifying labels for the input bigWig files, this causes deepTools to use the file name after removing the path and extension.
 
--outRawCounts	
#Save average scores per region for each bigWig file to a single tab-delimited file.

-out
#File name to save the compressed matrix file (npz format) needed by the “plotPCA” and “plotCorrelation” tools.

#conda create -n deeptools -c bioconda deeptools
#conda activate deeptools
#multiBigwigSummary BED-file --BED potential_enhancer_coordinates_mm39.bed -b  /Users/roucogar/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/6_Chip-seq/broad_trial3/H3K27ac-FL-E105-GA58_DP_Rep1_macs_SR200_norm.bw /Users/roucogar/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/6_Chip-seq/broad_trial3/H3K27ac-FL-E135-GA58_DP_Rep1_macs_SR200_norm.bw --smartLabels -out norm_coverage_enhancers_merged500bp.npz --outRawCounts norm_coverage_enhancers_merged500bp.bed
multiBigwigSummary BED-file --BED H3K27ac-FL_GA58_DP_macs_SR200_broadpeaks_merged_Shox2_TAD.bed -b H3K27ac-FL-E105-GA58_DP_Rep1_macs_SR200_norm.bw H3K27ac-FL-E115-GA58_DP_Rep1_macs_SR200_norm.bw H3K27ac-FL-E125-GA58_DP_Rep1_macs_SR200_norm.bw H3K27ac-FL-E135-GA58_DP_Rep1_macs_SR200_norm.bw --smartLabels -out norm_coverage_enhancers_merged_Shox2TAD.npz --outRawCounts norm_coverage_enhancers_merged_Shox2TAD.bed
#conda deactivate

#Go to R
#save norm_coverage_enhancers_merged_Shox2TAD.bed as xslm tipe of file to open it on R
#in R run:
library(readxl)
library(dplyr)
setwd("~/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/6_Chip-seq/broad_trial3")

norm_coverage_enhancers_merged_Shox2TAD_bed <- read_excel("/Users/roucogar/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/6_Chip-seq/broad_trial3/norm_coverage_enhancers_merged_Shox2TAD.xlsm")

enh <-data.frame(norm_coverage_enhancers_merged_Shox2TAD_bed)

	#filter by values of 0.3

enh_filter_0_3 <- filter(enh,enh$X.H3K27ac.FL.E105.GA58_DP_Rep1_macs_SR200_norm. >= 0.3 | 
                           enh$X.H3K27ac.FL.E115.GA58_DP_Rep1_macs_SR200_norm. >= 0.3 |
                           enh$X.H3K27ac.FL.E125.GA58_DP_Rep1_macs_SR200_norm. >= 0.3 |
                           enh$X.H3K27ac.FL.E135.GA58_DP_Rep1_macs_SR200_norm. >= 0.3)

write.table(enh_filter_0_3, "potential_Shox2_enhancers_filter0_3_with_values.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

coor_enh_filter_0_3 <- enh_filter_0_3[,1:3]

write.table(coor_enh_filter_0_3, "potential_Shox2_enhancers_filter0_3.bed", sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


#outside R filter de fragments that smaller than 600bp.
#We are considering that an enhancer should have at least 600bp size to be surrounded by at least two nucleosomes

awk '{if($3-$2 >= 600) print}' potential_Shox2_enhancers_filter0_3.bed > potential_Shox2_enhancers_filter0_3_600pinterval.bed

awk '{if($3-$2 >= 600) print}' potential_Shox2_enhancers_filter0_3_with_values.txt > potential_Shox2_enhancers_filter0_3_600bpintereval_with_values.txt











