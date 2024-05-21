#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out
#SBATCH -e slurm-%x-%A_%2a.err
#SBATCH --nodes 1
#SBATCH --mem 50G
#SBATCH --cpus-per-task 14
#SBATCH --time 08:00:00
#SBATCH --partition=shared-cpu
#SBATCH --array=25,27,38,41,46,53,60,66 #identifie les lignes/échantillons à traiter depuis Table.txt
#SBATCH --job-name CHIP
#SBATCH --chdir /home/users/d/darbellf/live/ChIPseq/fastq/processing

# This script remove adapters 
# make alignment with Bowtie2
# select MAPQ30 alignments
# Compute coverage SR200
# Normalize by nb of tags

path="$PWD/"
pathForFastq="$path/"
genome=mm39
#adapterSeq="TruSeq" #mask the non-relevant adapter type
adapterSeq="NextSeq"

pathForTable="${path}/${genome}/1b_Narrow_peak_MACS2_GEOdowloaded_table_example.txt"

pathForFasta="/home/users/d/darbellf/live/genomes/${genome}/"
pathForIndex="/home/users/d/darbellf/live/genomes/${genome}/bowtie2/bowtie2" #mind the pathbase, don't go with the index folder (eg. "/home/users/d/darbellf/live/genomes/${genome}/bowtie2/" is not recognized by bowtie2)

indexPath=${pathForIndex}

nbOfThreads=14

export PATH=$PATH:/home/users/d/darbellf/live/softwares/ucsc_utilities/ #to load bedGraphToBigWig
export PATH=$PATH:/home/users/d/darbellf/live/softwares/sratoolkit.2.11.2-centos_linux64/bin #to load fasterq-dump from sratoolkit/2.11.2

module purge
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load cutadapt/1.18-Python-3.6.6

sample=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $1}')
sra=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $2}') #if SRA then provide SRR in $2Table
fastqFile=${sample}.fastq.gz #if sra
#fastqFile=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $2}') #use if fastq are stored locally

mkdir -p ${path}/${genome}/${sample}

pathResults=${path}/${genome}/${sample}/
echo $sample
cd $pathResults


# Cutadapt
if [ ! -e ${path}/${genome}/allFinalFiles/reports/${sample}_report-cutadapt.txt ]; then
  fastq=${pathForFastq}/$fastqFile
  if [ ! -e $fastq ]; then
    mkdir -p $pathForFastq
    cd $pathForFastq
    fasterq-dump -o ${sample}.fastq ${sra}
    gzip ${sample}.fastq
    cd $pathResults
  fi
  if [ ! -s $fastq ]; then
    echo "FASTQ IS EMPTY"
    exit 1
  fi
  if [ $adapterSeq = "TruSeq" ]; then
    cutadapt -j $nbOfThreads -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 30 -m 15 -o ${pathResults}${sample}-cutadapt.fastq.gz $fastq > ${pathResults}${sample}_report-cutadapt.txt
  else
    if [ $adapterSeq = "NextSeq" ]; then
      cutadapt -j $nbOfThreads -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -q 30 -m 15 -o ${pathResults}${sample}-cutadapt.fastq.gz $fastq > ${pathResults}${sample}_report-cutadapt.txt
    else
      echo "YOU NEED TO WRITE THE CODE"
      exit 1
    fi
  fi
  
  mkdir -p ${path}/${genome}/allFinalFiles/reports/
  cp ${pathResults}${sample}_report-cutadapt.txt ${path}/${genome}/allFinalFiles/reports/${sample}_report-cutadapt.txt
fi

module purge
module load GCC/8.3.0 #required for samtools, bowtie2 and bedtools
module load OpenMPI/3.1.4
module load SAMtools/1.10
module load BEDTools/2.28.0
module load Bowtie2/2.3.5.1

# Mapping
if [ ! -e ${path}/${genome}/allFinalFiles/reports/${sample}_mapping_stats.txt ];then
  bowtie2 -p $nbOfThreads -x $indexPath -U ${pathResults}${sample}-cutadapt.fastq.gz 2> ${pathResults}${sample}_mapping_stats.txt  | samtools view --threads $nbOfThreads -Su - | samtools sort --threads $nbOfThreads -o ${pathResults}${sample}_mapped_sorted.bam
  mkdir -p ${path}/${genome}/allFinalFiles/reports/
  cp ${pathResults}${sample}_mapping_stats.txt ${path}/${genome}/allFinalFiles/reports/
fi

# MAPQ30
if [ ! -e ${pathResults}${sample}_mapped_sorted_q30.bam ]; then
  samtools view --threads $nbOfThreads -b ${pathResults}${sample}_mapped_sorted.bam -q 30 > ${pathResults}${sample}_mapped_sorted_q30.bam
fi

mkdir -p ${path}/${genome}/allFinalFiles/bedGraphs

if [ ! -e ${pathForFasta}${genome}.fa.fai ]; then
  samtools faidx ${pathForFasta}${genome}.fa
fi

module purge
module load GCC/4.9.3-2.25
module load OpenMPI/1.10.2
module load MACS2/2.1.1.20160309-Python-2.7.11

# Use macs2 with fixed fragment size of 200bp
if [ ! -e ${path}bedGraphs/${sample}_macs_SR200.bedGraph.gz ]; then
  macs2 callpeak -t ${sample}_mapped_sorted_q30.bam -n ${sample}_macs_SR200 --call-summits --nomodel --extsize 200 -B 2> ${pathResults}${sample}_macs_SR200.log
  cp ${pathResults}${sample}_macs_SR200.log ${path}/${genome}/allFinalFiles/reports/
  cp ${pathResults}${sample}_macs_SR200_peaks.narrowPeak ${path}/${genome}/allFinalFiles/bedGraphs
  # The output of macs2 can go over the size of chromosomes so a special care needs to be taken before coverting to bigwig
  bash /home/users/d/darbellf/live/softwares/fromMacs2BdgToSimplifiedBdgAndBw.sh ${pathResults}${sample}_macs_SR200_treat_pileup.bdg ${path}/${genome}/allFinalFiles/bedGraphs/${sample}_macs_SR200 "macs2 SR200 of ${sample}" ${pathForFasta}${genome}.fa.fai &
fi
wait

if [ ! -e ${path}bedGraphs/${sample}_macs_SR200_norm.bedGraph.gz ]; then
  nbtags=$(grep "total tags in treatment" ${pathResults}${sample}_macs_SR200.log | awk '{print $NF}')
  zcat ${path}/${genome}/allFinalFiles/bedGraphs/${sample}_macs_SR200.bedGraph.gz | awk -v s=$sample -v n=$nbtags -v OFS="\t" 'BEGIN{print "track type=bedGraph name=\""s" SR200 normalized by million tags\" visibility=full autoScale=on windowingFunction=mean"}NR>1{$4=$4/n*1e6; print}' > ${path}/${genome}/allFinalFiles/bedGraphs/${sample}_macs_SR200_norm.bedGraph
  bedGraphToBigWig ${path}/${genome}/allFinalFiles/bedGraphs/${sample}_macs_SR200_norm.bedGraph ${pathForFasta}/${genome}.fa.fai ${path}/${genome}/allFinalFiles/bedGraphs/${sample}_macs_SR200_norm.bw
  gzip ${path}/${genome}/allFinalFiles/bedGraphs/${sample}_macs_SR200_norm.bedGraph &
fi

wait

echo "Everything is done"
find . -size 0 -delete


