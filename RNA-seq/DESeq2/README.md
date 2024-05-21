#`DESeq2` analysis
All the steps provided here were adapted from https://github.com/lldelisle/rnaseq_rscripts to analyze datasets for Rouco et al. 2024

1st step:
We run the `1a_step1-generateTables.R` providing the corresponding configuration file as the example here provided `1b_configFileRNAseq_step1_FL58_allstages_example.R`.
To configurate this file extra input files are provided: `RNAseqFunctions.R`, `samplesPlan_example.txt` , `genesInchrMfrommergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm39.104_ExonsOnly_UCSC_dsmCherry_P2A_CRE_EYFP.txt`.

2nd step:
We run the `2a_step2-DESeq2.R` providing the corresponding configuration file as the example here provided `2b_configFileRNAseq_step2_example.R`.
As input files we provided the `AllHTSeqCounts_subset.txt` obtained from the previous step containing the counts and the `AllCufflinks_Simplified_norm.txt` also obtained in the previous stept to provide as a table with annotations.
To configurate this file extra input files are provided: `RNAseqFunctions.R`, `samplesPlan_example.txt`.
The input file `mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm39.104_ExonsOnly_UCSC_dsmCherry_P2A_CRE_EYFP.gtf` can be downloaded from https://zenodo.org/records/11219861

3rd step:
We run the `3a_step3-graphClusteringPCAGenes.R` providing the corresponding configuration file as the example here provided `3b_configFileRNAseq_step3_example.R`.
As input files we provided the `DESeq2Analysis_PerGeno.txt` obtained from the previous step.
To configurate this file extra input files are provided: `RNAseqFunctions.R`, `samplesPlan_example.txt` ,`relevant_genes.txt`.

4th step:
We run the `4a_step4-graphVolcanoAndMAP.R` providing the corresponding configuration file as the example here provided `4b_configFileRNAseq_step4_example.R`.
As input files we provided the `DESeq2Analysis_PerGeno.txt` obtained from the previous step.
To configurate this file extra input files are provided: `RNAseqFunctions.R`, `relevant_genes.txt`.
