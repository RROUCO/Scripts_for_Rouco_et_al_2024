RNA-seq FASTQ files from two replicates of entire forelimbs at E10.5 and E13.5 (Andrey et al., 2017) were re-analyzed following the RNA-seq pipeline described in the `RNA-seq` folder of this repository.
Output files after DESeq2 analysis were used as input files for the first step of this pipeline.

ChIP-seq H3K27Ac datasets of entire forelimbs at E10.5 and E13.5 (Andrey et al., 2017) were first reanalyzed following the ChIP-seq pipeline described in the `ChIP-seq` folder of this repository.
H3K27ac MACS2 narrowpeaks were used as input for the second step of this pipeline.

1st step: Filter genes by being stably expressed between E10.5 and E13.5
The input file "AllCufflinks_Simplified_norm.txt" obtained from the DESeq2 pipeline was used as input together with the manually curated list of genes "1b_Curated_Gene_list.txt" to filter developmental genes of interest FPKMs.
The average of normalized FPKM values was calculated and used to compute the ratio among E10.5 and E13.5 datasets. 
Since we were interested in genes having stable expression between E10.5 and E13.5, genes were filtered to keep those with FPKM values bigger than 5, at both stages.
Then, we excluded genes having a fold change larger than 3 between the two stages. 

2nd step: Filtering ChIP-seq peaks
H3K27ac MACS2 narrowpeaks were restricted within the interaction domain defined by promoter Capture-C (Andrey et al., 2017) using "2b_Filtered_Curated_gene_and_interaction_domain_mm39.bed" file generated afer manually assign the interaction domains from (Andrey et al., 2017) to the output list from the previous step "Filtered_Curated_Gene_list_avg_ratio.csv". 
Then H3K27ac peaks within protein-coding gene promoters were removed using the file "2c_exclusionAroundTSS_prot_coding.bed".
Remaining peaks were extended by +/- 300bp and merged.
H3K27Ac peaks were then classified as putative common enhancers when present in both E10.5 and E13.5 using bedtools intersect. 
H3K27Ac peaks present only in the E10.5 dataset were classified as putative early enhancers while H3K27Ac peaks present only in the E13.5 dataset were classified as putative late enhancers.

3rd step: Assign putative enhancers to gene interaction domain
Putative enhancers were then assigned to gene interaction domain. Run in R.

Out put file "domains_mm39.csv" from step3 was used to produce Figure 1A using "Figure1a_plot_enhancer_category_across_loci.R"
