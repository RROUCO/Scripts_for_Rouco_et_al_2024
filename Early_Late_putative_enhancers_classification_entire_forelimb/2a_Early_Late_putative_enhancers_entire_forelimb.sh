# build beds from narrowPeak
cut -f 1-3 H3K27Ac-FL-E105-Wt-Mm-Rep1-L4288_macs_SR200_peaks.narrowPeak | uniq > H3K27Ac-FL-E105-Wt-Rep1_peaks.bed
cut -f 1-3 H3K27Ac-FL-E135-Wt-Mm-Rep2-L7110_macs_SR200_peaks.narrowPeak | uniq > H3K27Ac-FL-E135-Wt-Rep2_peaks.bed

# restrict H3K27ac peaks to selected interaction domains (liftedOver online from mm9 to mm39)
bedtools intersect -a H3K27Ac-FL-E105-Wt-Rep1_peaks.bed -b ~/Dropbox/Trajectories_Rouco_et_al/Limb_wide_early-late/Filtered_Curated_gene_and_interaction_domain_mm39.bed -wa | uniq > H3K27Ac-FL-E105-Wt-Rep1_id_peaks.bed
bedtools intersect -a H3K27Ac-FL-E135-Wt-Rep2_peaks.bed -b ~/Dropbox/Trajectories_Rouco_et_al/Limb_wide_early-late/Filtered_Curated_gene_and_interaction_domain_mm39.bed -wa | uniq > H3K27Ac-FL-E135-Wt-Rep2_id_peaks.bed

### #exclusion of protein coding gene promoters
# write 1st exon of protein_coding genes
#grep "exon_number \"1\"" /Users/nccrgenetics/Dropbox/analysis/132_Col2a1_sensor/peaks_analysis/mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm39.104_ExonsOnly_UCSC.gtf | grep protein_coding > first_exons_prot_coding_mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm39.104_ExonsOnly_UCSC.gtf
# write left flanking coordinates of protein_coding gene 1st exon (TSS)
#bedtools flank -s -l 10 -r 0 -i first_exons_prot_coding_mergeOverlapGenesOfFilteredTranscriptsOfMus_musculus.GRCm39.104_ExonsOnly_UCSC.gtf -g /Users/nccrgenetics/igv/genomes/seq/mm39.fa.fai > tss_prot_coding.gtf
# reformat .gtf to .bed
#awk -F "\t" -v OFS="\t" '{split($9, a, ";"); for(i in a){if(a[i]~/gene_name/){split(a[i], b, "\"")}};print $1, $4, $5, b[2]}' input.gtf > output.bed
#awk -v OFS="\t" '{print $1, $4, $5, ".", ".", $7}' tss_prot_coding.gtf > tss_prot_coding.bed
# exclude promoter region (+ is downstream (right),- is upstream (left), convention is -2kb/+500bp)
#bedtools slop -i tss_prot_coding.bed -g /Users/nccrgenetics/igv/genomes/seq/mm39.fa.fai -l 2000 -r 500 -s > exclusionAroundTSS_prot_coding.bed
###

# remove H3K27ac peaks within protein-coding gene promoters (+500bp/-2kb)
bedtools intersect -a H3K27Ac-FL-E105-Wt-Rep1_id_peaks.bed -b exclusionAroundTSS_prot_coding.bed -v > H3K27Ac-FL-E105-Wt-Rep1_id_peaks_outsideTSS_prot_coding.bed
bedtools intersect -a H3K27Ac-FL-E135-Wt-Rep2_id_peaks.bed -b exclusionAroundTSS_prot_coding.bed -v > H3K27Ac-FL-E135-Wt-Rep2_id_peaks_outsideTSS_prot_coding.bed

# extended H3K27ac peaks by +/-300bp
bedtools slop -i H3K27Ac-FL-E105-Wt-Rep1_id_peaks_outsideTSS_prot_coding.bed -g /Users/nccrgenetics/igv/genomes/seq/mm39.fa.fai -b 300 > H3K27Ac-FL-E105-Wt-Rep1_id_peaks_outsideTSS_prot_coding_300bp_extended.bed
bedtools slop -i H3K27Ac-FL-E135-Wt-Rep2_id_peaks_outsideTSS_prot_coding.bed -g /Users/nccrgenetics/igv/genomes/seq/mm39.fa.fai -b 300 > H3K27Ac-FL-E135-Wt-Rep2_id_peaks_outsideTSS_prot_coding_300bp_extended.bed

# merge extended H3K27ac peaks, for each samples
bedtools merge -i H3K27Ac-FL-E105-Wt-Rep1_id_peaks_outsideTSS_prot_coding_300bp_extended.bed > H3K27Ac-FL-E105-Wt-Rep1_id_peaks_outsideTSS_prot_coding_300bp_extended_merged.bed
bedtools merge -i H3K27Ac-FL-E135-Wt-Rep2_id_peaks_outsideTSS_prot_coding_300bp_extended.bed > H3K27Ac-FL-E135-Wt-Rep2_id_peaks_outsideTSS_prot_coding_300bp_extended_merged.bed

# print common extended merged H3K27ac peaks
bedtools intersect -a H3K27Ac-FL-E105-Wt-Rep1_id_peaks_outsideTSS_prot_coding_300bp_extended_merged.bed -b H3K27Ac-FL-E135-Wt-Rep2_id_peaks_outsideTSS_prot_coding_300bp_extended_merged.bed | uniq | awk -v OFS="\t" '{ print $1,$2,$3,"C","0",".",$2,$3,"135,206,235"}' > H3K27Ac_FL_common_peaks_300bp_extended_merged.bed

# print stage-specific extended H3K27ac peaks 
bedtools intersect -a H3K27Ac-FL-E105-Wt-Rep1_id_peaks_outsideTSS_prot_coding_300bp_extended_merged.bed -b H3K27Ac_FL_common_peaks_300bp_extended_merged.bed -v | uniq | awk -v OFS="\t" '{ print $1,$2,$3, "E","0",".",$2,$3,"80,200,120"}' > H3K27Ac_FL_E105_peaks_300bp_extended_merged.bed
bedtools intersect -a H3K27Ac-FL-E135-Wt-Rep2_id_peaks_outsideTSS_prot_coding_300bp_extended_merged.bed -b H3K27Ac_FL_common_peaks_300bp_extended_merged.bed -v | uniq | awk -v OFS="\t" '{ print $1,$2,$3, "L","0",".",$2,$3,"191,64,191"}' > H3K27Ac_FL_E135_peaks_300bp_extended_merged.bed

# assemble into a final bed file for graphical visualisation
cat H3K27Ac_FL_common_peaks_300bp_extended_merged.bed H3K27Ac_FL_E105_peaks_300bp_extended_merged.bed H3K27Ac_FL_E135_peaks_300bp_extended_merged.bed |  LC_ALL=C sort -k1,1 -k2,2n > H3K27ac_FL_classified_peaks.bed

# move into R