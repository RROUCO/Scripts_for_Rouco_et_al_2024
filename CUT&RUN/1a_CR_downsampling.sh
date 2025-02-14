#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out
#SBATCH -e slurm-%x-%A_%2a.err
#SBATCH --nodes 1
#SBATCH --mem 50G
#SBATCH --cpus-per-task 24
#SBATCH --time 03:00:00
#SBATCH --partition=shared-cpu
#SBATCH --array=1-8 #identifie les lignes/échantillons à traiter depuis Table.txt
#SBATCH --job-name downsampling_CR
#SBATCH --chdir /home/roucogar/scratch/CUT_RUN/

path="$PWD/"
pathForFastq="$path/"

pathForTable="${path}/1b_CR_downsampling_table_example.txt"

export PATH=$PATH:/home/roucogar/live/softwares/seqtk/

Fastq=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $1}')
#sra=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $2}') #if SRA then provide SRR in $2Table
#fastqFile=${sample}.fastq.gz #if sra
fastqFile=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $2}')  #use if fastq are stored locally

seqtk sample -s100 $fastqFile 78000000 > $Fastq
gzip $Fastq