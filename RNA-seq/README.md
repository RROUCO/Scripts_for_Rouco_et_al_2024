To analyze` RNA-Seq` FASTQ files the following pipelines were adapted from https://github.com/lldelisle/myNGSanalysis/ and used depending on the type of dataset. 

For datasets from GEO that we wanted to reanalyze (Andrey et al., 2017) we used the sbach command line `1a_GEO_datasets_RNAseq_analysis.sh` providing the FASTQ files through an external table `1b_GEO_datasets_RNAseq_analysis_table.txt`.

For datasets produced in this study either from entire limb or from FACS sorted cells we used the sbach command line `2a_RRouco_datasets_RNAseq.sh` providing the FASTQ files through an external table like the one presented in this repository `2b_RRouco_datasets_RNAseq_table_example.txt`.
Custom genome used for processing this samples can be found in https://zenodo.org/records/11219861

`DESeq2` analysis:
Output files from previous steps were processed adapting the pipeline from https://github.com/lldelisle/rnaseq_rscripts
Rscripts with a representative example of the configuration files used and other required input files are provided