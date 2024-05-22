#!/bin/sh
#SBATCH -J sc58_E105_1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem 64G
#SBATCH -p shared-cpu
#SBATCH -t 10:00:00

module load CellRanger/6.1.2

cellranger count --id=cr_HL_E105_58_REP1 --transcriptome=/home/share/andrey_lab/Genomes/cellranger/cr_mm39_dsmCherry_P2A_CRE_EYFP --fastqs=/home/roucogar/live/single_cell/fastq/processing --sample=HLREP1E105

