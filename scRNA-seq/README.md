1st Step:
scRNA-seq FASTQ files obtained from 10X experiments were processed with 10× Genomics Cell Ranger software (version 6.1.2) using the sbatch `Cellranger_analysis.sh` provided in this folder for each sample independently.
Custom genome used for processing this samples can be found in https://zenodo.org/records/14865689

2nd Step:
Cell Ranger output files for each dataset were further processed using the velocyto run10x command, as described below, from the velocyto.py tool (version 0.17.17) in Python (version 3.9.12) with our custom genome GTF (https://zenodo.org/records/11219861) and the UCSC genome browser repeat masker.gtf file to mask expressed repetitive elements to generate a loom file for each sample.
velocyto run10x -m /Volumes/RAQUEL-LAB3/genomes/velocyto/mm39_rmsk.gtf /Volumes/RAQUEL-LAB3/Single_cell_58_data/CellRanger_58/cr_HL_E105_58_REP1/ /Volumes/RAQUEL-LAB3/genomes/velocyto/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm39.104_ExonsOnly_UCSC_dsmCherry_P2A_CRE_EYFP.gtf\

3rd Step:
Each resulting loom matrix, comprising spliced/unspliced/ambiguous reads, was individually imported into R (version 4.1.2) using the Read Velocity function from the Seurat Wrappers package (version 0.3.0). 
Simultaneously, feature-filtered output matrices obtained from Cell Ranger were loaded into R separately through the Read10X function of the Seurat package (version 4.2.1).
Subsequently, the spliced, unspliced, ambiguous, and RNA feature data were combined into a single matrix for each dataset. 
Following this, each matrix was transformed into a Seurat object using the Seurat package. 
These steps and downstream scRNA-seq analysis can be found in the file `scRNA-seq_Seurat_Rouco_et_al_2024_updated04032024.R`.