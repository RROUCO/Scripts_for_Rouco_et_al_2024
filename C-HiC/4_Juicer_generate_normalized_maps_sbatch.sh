#!/bin/sh

#SBATCH -J CHIC-58Dp
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -p shared-cpu
#SBATCH -t 02:00:00
#SBATCH --mem-per-cpu=10000

# stop immediately, when a program return non-zero results (i.e. an error)
set -e

EXPORTER="java -jar /home/roucogar/Software/C-HiC_workflow_Robert/CHIC_exporter_java/dist/lib/CaptureC_exporter.jar"
JUICER_TOOLS="java -jar /home/roucogar/Software/juicer_tools.1.9.9_jcuda.0.8.jar"

# 1 based coordinates of the enriched region, in round down to next complete bin, if not multiple of the binsize
CHROM="chr3"
START=65000001
END=68500000
GENOME="capture_chrom.sizes"
BINSIZES="5000 10000"
# Minimal MAPQ value, that is required (di-tags below are ignored for the create of the Hi-C map)
MAPQ="30"
CHROM_SIZE=$(($END-$START+1))

echo $CHROM"	"$CHROM_SIZE > $GENOME

DPNII_STARTS="/home/roucogar/Genomes/Juicerfragmented/DpnII_starts_mm9.txt"


BAMS="/srv/beegfs/scratch/users/r/roucogar/CHiC_Robert_pipeline/C58ESCs2/CHiC_ESCs2_GA58_S1_L001_R1_2_001.hicup.bam"

WASHU_DIR=Maps_WashU

mkdir -p txt
mkdir -p $WASHU_DIR
mkdir -p hic

for BAM in $BAMS; do

	###
	# Writing of BAM fite to Text format
	###
	TXT=`basename $BAM .bam`".txt"
	echo "Export to Juicebox intermediate file format"
	# Line 41 of the header has an inproper key value structure for Picard-lib
	# Due to option LENIENT in SAMReader construction it is ignored
	$EXPORTER $BAM $DPNII_STARTS txt/$TXT 2> txt/$TXT.sorted_out

	###
	# Shifting of the coordinates AND filtering of the chromosome region
	###
	SHIFTED=`basename $TXT .txt`".shifted.txt"
	echo "Shift coordinates, filter entries"
	# kick out line 41 of the header because has a inproper key value structure for Picard-lib and is not needed
	less txt/$TXT | awk -v chrom=$CHROM -v start=$START -v end=$END '
		BEGIN {
			offset = start-1;
			# 1 based intervals
			range = end - start + 1; 
		}
		{
		if (NF != 11) {
		    print "Wrong input format, line:", NR
		    exit 1
		}
		# readname str1 chr1 pos1 frag1 str2 chr2 pos2 frag2 mapq1 mapq2
		if($3==chrom && $7==chrom && 
		   $4 >= start && $4 <= end && 
		   $8 >= start && $8 <= end) {
			# substract offset from the position coordinates
			# -1 because Juicebox coordinates start at 0
			pos1 = $4-offset-1
			pos2 = $8-offset-1
			print $1, $2, $3, pos1, $5, $6, $7, pos2, $9, $10, $11;

			# sanity checks
			if(pos1 < 0 || pos1 >= range || pos2 < 0 || pos2 >= range ) {
				print "ERROR in input line:", NR
			    exit 1
			}
		}
	}' > txt/$SHIFTED


	####
	# Converting the text file into a HiC map
	####
	HIC=`basename $TXT .txt`".MAPQ$MAPQ.hic"
	# Juicebox hic maps
	# mapping quality has to be at least MAPQ
	$JUICER_TOOLS pre -q $MAPQ txt/$SHIFTED hic/$HIC $GENOME


	for BINSIZE in $BINSIZES; do

		KB=$(($BINSIZE/1000))"kb"
		mkdir -p $WASHU_DIR/$KB
		####
		# Exporting of the balanced HiC map 
		####
		BALANCED=`basename $HIC .hic`".KR_$(($BINSIZE/1000))kb.txt"
		echo "Exporting balanced map."
		# Juicebox hic maps
		$JUICER_TOOLS dump observed KR hic/$HIC $CHROM $CHROM BP $BINSIZE hic/$BALANCED

		####
		# Exporting of the raw HiC map 
		####
		RAWCOUNT=`basename $HIC .hic`".Raw_$(($BINSIZE/1000))kb.txt"
		echo "Exporting raw count map."
		# Juicebox hic maps
		$JUICER_TOOLS dump observed NONE hic/$HIC $CHROM $CHROM BP $BINSIZE hic/$RAWCOUNT


		####
		# Shifting back the coordinates to original
		####
		FINALBALANCED=`basename $BALANCED .txt`.WashU.txt
		less hic/$BALANCED | awk -v chrom=$CHROM -v start=$START -v binsize=$BINSIZE ' BEGIN {
			offset = start-1;
		}
		{
			print chrom "," $1+offset+1 "," $1+offset+binsize "\t" chrom ":" $2+offset+1 "-" $2+offset+binsize "\t" $3 
		}' > "$WASHU_DIR"/"$KB"/"$FINALBALANCED"


		####
		# Shifting back the coordinates to original
		####
		FINALRAW=`basename $RAWCOUNT .txt`.WashU.txt
		less hic/$RAWCOUNT | awk -v chrom=$CHROM -v start=$START -v binsize=$BINSIZE ' BEGIN {
			offset = start-1
		}
		{
			print chrom "," $1+offset+1 "," $1+offset+binsize "\t" chrom ":" $2+offset+1 "-" $2+offset+binsize "\t" $3 
		}' > "$WASHU_DIR"/"$KB"/"$FINALRAW"
	done

done

# cleanup
rm hic/*.txt
rm txt/*.hicup.txt
rm txt/*.hicup.shifted.txt



