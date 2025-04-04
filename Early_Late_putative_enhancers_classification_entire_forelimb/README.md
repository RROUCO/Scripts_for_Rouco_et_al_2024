`RNA-seq `FASTQ files from two replicates of entire forelimbs at E10.5 and E13.5 (Andrey et al., 2017) were re-analyzed following the `RNA-seq pipeline` described in the `RNA-seq folder` of this repository.
Output files, obtained from running the `DESeq2 analysis` following the pipeline described in the `RNA-seq folder`, were used as input files for the first step of this pipeline.

`ChIP-seq` FASTQ H3K27Ac datasets of entire forelimbs at E10.5 and E13.5 (Andrey et al., 2017) were first reanalyzed following the `ChIP-seq pipeline` described in the `ChIP-seq folder` of this repository.
H3K27ac MACS2 narrowpeaks were used as input for the second step of this pipeline.

1st step: Filter genes by being stably expressed between E10.5 and E13.5
The input file `AllCufflinks_Simplified_norm.txt` obtained from the DESeq2 pipeline was used as input together with the manually curated list of genes `1b_Curated_Gene_list.txt` to filter developmental genes of interest FPKMs using the `1a_Limb_wide_FPKM_Extraction_eFL_E105_E135.R`pipeline.
The average of normalized FPKM values was calculated and used to compute the ratio among E10.5 and E13.5 datasets. 
Since we were interested in genes having stable expression between E10.5 and E13.5, genes were filtered to keep those with FPKM values bigger than 5, at both stages.
Then, we excluded genes having a fold change larger than 3 between the two stages. 

2nd step: Filtering ChIP-seq peaks
H3K27ac MACS2 `.narrowpeaks` were restricted within the interaction domain defined by promoter Capture-C (Andrey et al., 2017) using `2b_Filtered_Curated_gene_and_interaction_domain_mm39.bed` file generated afer manually assign the interaction domains from (Andrey et al., 2017) to the output list from the previous step `Filtered_Curated_Gene_list_avg_ratio.csv`. 
Then H3K27ac peaks within protein-coding gene promoters were removed using the file `2c_exclusionAroundTSS_prot_coding.bed`. At the Shox2 locus, peak called at the alternative Veph1 alternative promoter was manually excluded. 
Remaining peaks were extended by +/- 300bp and merged.
H3K27Ac peaks were then classified as putative common enhancers when present in both E10.5 and E13.5 using bedtools intersect. 
H3K27Ac peaks present only in the E10.5 dataset were classified as putative early enhancers while H3K27Ac peaks present only in the E13.5 dataset were classified as putative late enhancers.

Output file  `H3K27ac_FL_classified_peaks.bed` from step 2 containing the list of 1'625 putative enhancers classified was used to produce the tornado plot from Figure 1A using the command lines described in `Figure1A_Tornado_plot_from_classified_enhancers.sh`.

3rd step: Assign putative enhancers to gene interaction domain
Putative enhancers were then assigned to gene interaction domain. Running in R `3_domains_enhancers_stage.R`

Output file `domains_mm39.csv` from step3 was used to produce Supplementary Figure 1A using `Supplementary_Figure1A_plot_enhancer_category_across_loci.R`

