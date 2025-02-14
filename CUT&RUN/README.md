CUT & RUN paired FASTQ files of each sample were processed to subsample 78 M read pairs using seqtk toolkit as described in `1a_CR_downsampling.sh`

After subsampling, to analyze subsampled CUT&RUN FASTQ files the following pipeline was adapted from https://github.com/lldelisle/myNGSanalysis/blob/main/hpc/CUTandRUN/CUTandRUN.sh

Sbach command line `2a_CUTandRUN.sh` providing the subsampled FASTQ files through an external table like the one provided in `2b_CUT_RUN_table_example.txt`

Custom genome used for processing these samples can be found in https://zenodo.org/records/14865689.

