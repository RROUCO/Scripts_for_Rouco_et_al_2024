#!/bin/bash

#SBATCH -o slurm-%x-%A_%2a.out
#SBATCH -e slurm-%x-%A_%2a.err
#SBATCH --nodes 1
#SBATCH --mem 50G
#SBATCH --cpus-per-task 16
#SBATCH --time 04:00:00
#SBATCH --partition=shared-bigmem
#SBATCH --array=1-34 #identifie les lignes/échantillons à traiter depuis Table.txt
#SBATCH --job-name RNAseq
#SBATCH --chdir /srv/beegfs/scratch/users/r/roucogar/RNA-Seq/Forelimb/

# This script make alignment with STAR ENCODE parameters
# Evaluate FPKM with cufflinks
# Coverage with bedtools

path="$PWD/"
pathForFastq="/srv/beegfs/scratch/users/r/roucogar/RNA-Seq/Forelimb/"
genome=mm39_dsmCherry_P2A_CRE_EYFP
#adapterSeq="TruSeq"
adapterSeq="NextSeq"
ensemblVersion=104
versionOfGtf="Mus_musculus.GRCm39.$ensemblVersion"
gtfFile=/home/share/andrey_lab/Genomes/${genome}/filtered_gtf/mergeOverlapGenesOfFilteredTranscriptsOf${versionOfGtf}_ExonsOnly_UCSC_dsmCherry_P2A_CRE_EYFP.gtf

pathForTable="$path${genome}/2b_RRouco_datasets_RNAseq_table_example.txt"
pathForFasta="/home/share/andrey_lab/Genomes/${genome}/"
pathForIndex="/home/share/andrey_lab/Genomes/${genome}/STAR/"
indexPath=${pathForIndex}

# I tryed 36 and it failed:
# BAMoutput.cpp:27:BAMoutput: exiting because of *OUTPUT FILE* error: could not create output file ./_STARtmp//BAMsort/19/48
# SOLUTION: check that the path exists and you have write permission for this file. Also check ulimit -n and increase it to allow more open files.
nbOfThreads=16

export PATH=$PATH:/home/roucogar/live/softwares/ucsc_utilities/ #to load bedGraphToBigWig

module purge
module load GCC/7.3.0-2.30
module load OpenMPI/3.1.1
module load cutadapt/1.18-Python-3.6.6

sample=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $1}')
# sra=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $2}')
# fastqFile=${sample}.fastq.gz
fastqFile=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $2}') #complex
stranded=$(cat $pathForTable | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $3}') #write stranded or unstranded



if [ $SLURM_ARRAY_TASK_ID == 1 ] && [ ! -e ${path}/MTmouse.gtf ];then #mind the $SLURM_ARRAY_TASK_ID == 1 
  echo -e "chrM\tchrM_gene\texon\t0\t16299\t.\t+\t.\tgene_id \"chrM_gene_plus\"; transcript_id \"chrM_tx_plus\"; exon_id \"chrM_ex_plus\";">MTmouse.gtf
  echo -e "chrM\tchrM_gene\texon\t0\t16299\t.\t-\t.\tgene_id \"chrM_gene_minus\"; transcript_id \"chrM_tx_minus\"; exon_id \"chrM_ex_minus\";" >>MTmouse.gtf
fi

mkdir -p ${path}/${genome}/${sample}

pathResults=${path}/${genome}/${sample}/
echo $sample
cd $pathResults

if [ ! -e ${path}/${genome}/allFinalFiles/reports/${sample}_report-cutadapt.txt ]; then
  fastq=${pathForFastq}/$fastqFile
  if [ ! -s $fastq ]; then
    echo "FASTQ IS EMPTY"
    exit 1
  fi
  if [ $adapterSeq = "TruSeq" ]; then
    cutadapt -j $nbOfThreads -a GATCGGAAGAGCACACGTCTGAACTCCAGTCAC -q 30 -m 15 -o ${pathResults}${sample}-cutadapt.fastq.gz ${fastq} > ${pathResults}${sample}_report-cutadapt.txt
  else
    if [ $adapterSeq = "NextSeq" ]; then
    cutadapt -j $nbOfThreads -a CTGTCTCTTATACACATCTCCGAGCCCACGAGAC -q 30 -m 15 -o ${pathResults}${sample}-cutadapt.fastq.gz ${fastq} > ${pathResults}${sample}_report-cutadapt.txt #because_SR
     else
        echo "YOU NEED TO WRITE THE CODE"
        exit 1
    fi
  fi

  mkdir -p ${path}/${genome}/allFinalFiles/reports/
  cp ${pathResults}${sample}_report-cutadapt.txt ${path}/${genome}/allFinalFiles/reports/
fi

module purge
module load GCC/8.3.0 #required for samtools, star and bedtools
module load OpenMPI/3.1.4
module load OpenBLAS/0.3.7
module load R/3.6.2
module load STAR/2.7.2b
module load SAMtools/1.10
module load BEDTools/2.28.0
export PATH=$PATH:/home/roucogar/live/softwares/cufflinks-2.2.1.Linux_x86_64

if [ ! -e ${path}/${genome}/allFinalFiles/bam/${sample}_Aligned.sortedByCoord.out.bam ];then
  if [ "$stranded" = "unstranded" ]; then 
    #I need to add --outSAMstrandField intronMotif because it is not stranded library
    STAR --runThreadN $nbOfThreads --genomeDir ${indexPath} --readFilesIn ${pathResults}${sample}-cutadapt.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate  --outSAMstrandField intronMotif  --sjdbOverhang '99' --sjdbGTFfile $gtfFile  --quantMode GeneCounts  --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1
  else
    STAR --runThreadN $nbOfThreads --genomeDir ${indexPath} --readFilesIn ${pathResults}${sample}-cutadapt.fastq.gz --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate  --sjdbOverhang '99' --sjdbGTFfile $gtfFile  --quantMode GeneCounts  --outFilterType BySJout --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1
  fi
  mkdir -p ${path}/${genome}/allFinalFiles/reports/
  cp Log.final.out ${path}/${genome}/allFinalFiles/reports/${sample}_STAR_logFinal.txt
  mkdir -p ${path}/${genome}/allFinalFiles/bam/
  cp ${pathResults}Aligned.sortedByCoord.out.bam ${path}/${genome}/allFinalFiles/bam/${sample}_Aligned.sortedByCoord.out.bam
  samtools index ${path}/${genome}/allFinalFiles/bam/${sample}_Aligned.sortedByCoord.out.bam
fi

if [ -e ${pathForFasta}/${genome}.fa ] && [ -e ${path}MTmouse.gtf ] && [ -e $gtfFile ] && [ -s ${pathResults}Aligned.sortedByCoord.out.bam ];then
  if [ ! -e ${path}allFinalFiles/FPKM_${sample}_isoforms.txt ];then
    echo "export PATH=$PATH:/opt/ebsofts/MPI/GCC/4.9.3-2.25/OpenMPI/1.10.2/Cufflinks/2.2.1/bin" >cufflinks_${sample}.sh
    echo "mkdir -p ${pathResults}cufflinksWOMT" >>cufflinks_${sample}.sh #w/o MT
    echo "cufflinks -p 12 -o ${pathResults}cufflinksWOMT --max-bundle-length 10000000 --max-bundle-frags 100000000 --multi-read-correct --library-type \"fr-firststrand\" -b ${pathForFasta}/${genome}.fa  --no-effective-length-correction -M ${path}MTmouse.gtf -G $gtfFile ${pathResults}Aligned.sortedByCoord.out.bam" >>cufflinks_${sample}.sh
    echo "" >>cufflinks_${sample}.sh
    echo "mkdir -p ${path}/${genome}/allFinalFiles" >>cufflinks_${sample}.sh
    echo "cp ${pathResults}cufflinksWOMT/genes.fpkm_tracking ${path}/${genome}/allFinalFiles/FPKM_${sample}.txt" >>cufflinks_${sample}.sh
    echo "cp ${pathResults}cufflinksWOMT/isoforms.fpkm_tracking ${path}/${genome}/allFinalFiles/FPKM_${sample}_isoforms.txt" >>cufflinks_${sample}.sh
   
    if [ "$stranded" = "unstranded" ]; then 
      sed -i 's/fr-firststrand/fr-unstranded/g' cufflinks_${sample}.sh
    fi
    echo "Launching cufflinks"
    bash cufflinks_${sample}.sh &
  fi
else
  echo "cufflinks not launch because some files are missing."
fi

if { [ ! -e accepted_hits_unique_${sample}.bam ] && [ -s ${pathResults}Aligned.sortedByCoord.out.bam ] ;} || [ -e tmp.header ] ;then
  echo "Compute uniquely aligned"
  samtools view -H Aligned.sortedByCoord.out.bam >  tmp.header
  samtools view -@ 4 Aligned.sortedByCoord.out.bam | grep  -w "NH:i:1" | cat tmp.header - | samtools view -@ 5 -b > accepted_hits_unique_${sample}.bam
  rm tmp.header
fi

if [ ! -e ${path}/${genome}/allFinalFiles/htseqCount_${sample}.txt ] && [ -e ReadsPerGene.out.tab ];then 
  mkdir -p ${path}/${genome}/allFinalFiles
  echo "write htseqCount" #compile htseq counts-like from STAR counts
  if [ "$stranded" = "unstranded" ]; then 
    cat ReadsPerGene.out.tab | awk '{print $1"\t"$2}' > ${path}/${genome}/allFinalFiles/htseqCount_${sample}.txt
  else
    cat ReadsPerGene.out.tab | awk '{print $1"\t"$4}' > ${path}/${genome}/allFinalFiles/htseqCount_${sample}.txt
  fi
fi

if [ ! -e ${pathForFasta}${genome}.fa.fai ]; then
        samtools faidx ${pathForFasta}${genome}.fa
fi

#This is to make the bedgraph of coverage
mkdir -p ${path}/${genome}/allFinalFiles/bedGraphs
if { [ ! -e ${path}/${genome}/allFinalFiles/bedGraphs/${sample}_Uniq.bedGraph.gz ] && [ -e accepted_hits_unique_${sample}.bam ];} || [ -e tmp.header.u ] ;then
  echo "Building uniq reads bedGraph"
  echo "track type=bedGraph name=\"${sample} Uniq reads\"">tmp.header.u
  echo "bedtools genomecov -ibam accepted_hits_unique_${sample}.bam -bg -split | cat tmp.header.u - > ${path}/${genome}/allFinalFiles/bedGraphs/${sample}_Uniq.bedGraph">u.sh
  echo "cat ${path}/${genome}/allFinalFiles/bedGraphs/${sample}_Uniq.bedGraph | grep -v track | LC_ALL=C sort -k1,1 -k2,2n > ${path}/${genome}/allFinalFiles/bedGraphs/${sample}_Uniq_sorted.bedGraph">>u.sh
  echo "bedGraphToBigWig ${path}/${genome}/allFinalFiles/bedGraphs/${sample}_Uniq_sorted.bedGraph ${pathForFasta}/${genome}.fa.fai ${path}/${genome}/allFinalFiles/bedGraphs/${sample}_Uniq.bw">>u.sh
  echo "gzip ${path}/${genome}/allFinalFiles/bedGraphs/${sample}_Uniq.bedGraph">>u.sh
  echo "rm tmp.header.u">>u.sh
  echo "touch u.done">>u.sh
  bash u.sh &
fi
if [ "$stranded" != "unstranded" ]; then
  if { [ ! -e ${path}/${genome}/allFinalFiles/bedGraphs/${sample}_Uniq_fwd.bedGraph.gz ] && [ -e accepted_hits_unique_${sample}.bam ];} || [ -e tmp.header.uf ] ;then
  #strand + corresponds to reverse strand due to TruSeq
  echo "Building uniq fwd reads bedGraph"
  echo "track type=bedGraph name=\"${sample} Uniq reads forward\"">tmp.header.uf
  echo "bedtools genomecov -ibam accepted_hits_unique_${sample}.bam -bg -split -strand - | cat tmp.header.uf - > ${path}/${genome}/allFinalFiles/bedGraphs/${sample}_Uniq_fwd.bedGraph">uf.sh
  echo "cat ${path}/${genome}/allFinalFiles/bedGraphs/${sample}_Uniq_fwd.bedGraph | grep -v track | LC_ALL=C sort -k1,1 -k2,2n > ${path}/${genome}/allFinalFiles/bedGraphs/${sample}_Uniq_fwd_sorted.bedGraph">>uf.sh
  echo "bedGraphToBigWig ${path}/${genome}/allFinalFiles/${sample}_Uniq_fwd_sorted.bedGraph ${pathForFasta}/${genome}.fa.fai ${path}/${genome}/allFinalFiles/bedGraphs/${sample}_Uniq_fwd.bw">>uf.sh
  echo "gzip ${path}/${genome}/allFinalFiles/bedGraphs/${sample}_Uniq_fwd.bedGraph">>uf.sh
  echo "rm tmp.header.uf">>uf.sh
  echo "touch uf.done">>uf.sh
  bash uf.sh &
  fi
  if { [ ! -e ${path}/${genome}/allFinalFiles/bedGraphs/${sample}_Uniq_rev.bedGraph.gz ] && [ -e accepted_hits_unique_${sample}.bam ];} || [ -e tmp.header.ur ];then
  echo "Building uniq rev reads bedGraph"
  echo "track type=bedGraph name=\"${sample} Uniq reads reverse\"">tmp.header.ur
  echo "bedtools genomecov -ibam accepted_hits_unique_${sample}.bam -bg -split -strand + | cat tmp.header.ur - > ${path}/${genome}/allFinalFiles/bedGraphs/${sample}_Uniq_rev.bedGraph">ur.sh
  echo "cat ${path}/${genome}/allFinalFiles/bedGraphs/${sample}_Uniq_rev.bedGraph | grep -v track | LC_ALL=C sort -k1,1 -k2,2n > ${path}/${genome}/allFinalFiles/bedGraphs/${sample}_Uniq_rev_sorted.bedGraph">>ur.sh
  echo "bedGraphToBigWig ${path}/${genome}/allFinalFiles/${sample}_Uniq_rev_sorted.bedGraph ${pathForFasta}/${genome}.fa.fai ${path}/${genome}/allFinalFiles/bedGraphs/${sample}_Uniq_rev.bw">>ur.sh
  echo "gzip ${path}/${genome}/allFinalFiles/bedGraphs/${sample}_Uniq_rev.bedGraph">>ur.sh
  echo "rm tmp.header.ur">>ur.sh
  echo "touch ur.done">>ur.sh
  bash ur.sh &
  fi
fi

wait
echo "Everything is done"
find . -size 0 -delete