# stop upon occurrence of first error (non-negative return value of a program)
set -e

#WORKDIR="/path/to/output/profiles"

FRAG="/home/roucogar/Software/virtual-4C-Robert/profiles/DpnII_starts_mm9.txt"
VIEWPOINTS="/home/roucogar/Software/virtual-4C-Robert/profiles/Virtual_viewpoints_5kb.txt"
CHROMSIZES="/home/roucogar/Software/virtual-4C-Robert/profiles/mm9.chrom.sizes"

PROFILE_JAR="/home/roucogar/Software/virtual-4C-Robert/Virtual-Capture-C-Exporter/dist/lib/Virtual_Capture_C_Exporter.jar"

# 1 based coordinates for filtering the output of the exporter only for the enriched region
CHROM_ENRICHED="chr3"
START_ENRICHED="65000001"
END_ENRICHED="68500000"


# Generate Capture-C profiles
GENERATE_PROFILE() {
	BAM=$1
	OUT=$2
	
	cd $WORKDIR

	#########
	# Prepare directory
	#########

    rm -rf $OUT	
	mkdir -p $OUT
	cd $OUT

	# create files, clean in case old files are in it
	mkdir -p bam
	mkdir -p frag
	mkdir -p frag_norm
	mkdir -p binned
	mkdir -p smoothed_3kb
	mkdir -p smoothed_3kb_norm
	mkdir -p smoothed_5kb
	mkdir -p smoothed_5kb_norm

	###
	# Create Profiles
	# java -jar <infile.bam> <fragStartFile.txt> <chromSizes.txt> <viewpoints.txt> <outdir> <trackprefix> <excludedMargin>
	###
	java -jar -Xmx40024m $PROFILE_JAR $BAM $FRAG $CHROMSIZES $VIEWPOINTS "." $OUT 5000 $CHROM_ENRICHED $START_ENRICHED $END_ENRICHED 
                                                                                                  

	###
	# Process BAM Files
	# Sort and indexing
	###

	cd bam
	
	echo "--Sort BAM--"
	for BAM_OUT in $(ls *.bam); do
		
		BASE=$(basename $BAM_OUT .bam)
		samtools sort $BAM_OUT -o $BASE.sort.bam
		samtools index $BASE.sort.bam
		# remove the unsorted version of the bam file
		rm $BAM_OUT
	done

}


GENERATE_PROFILE /srv/beegfs/scratch/users/r/roucogar/CHiC_Robert_pipeline/CFL58DPrep3E115/CFL58D3115_S1_L001_R1_2_001.hicup.bam /srv/beegfs/scratch/users/r/roucogar/CHiC_Robert_pipeline/CFL58DPrep3E115/virtual4C_FL58D3115






