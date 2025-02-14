conda create --name ChIP
conda install -c conda-forge -c bioconda deeptools

#Split BED file with the 1'625 enhancers by putative enhancer type  Common, Early, Late
awk '$4 == "C" {print > "Common_H3K27ac_FL_classified_peaks.bed"}' H3K27ac_FL_classified_peaks.bed
awk '$4 == "E" {print > "Early_H3K27ac_FL_classified_peaks.bed"}' H3K27ac_FL_classified_peaks.bed
awk '$4 == "L" {print > "Late_H3K27ac_FL_classified_peaks.bed"}' H3K27ac_FL_classified_peaks.bed

#Matrix score + tornado plot
conda activate ChIP

computeMatrix reference-point \
--referencePoint center \
--missingDataAsZero \
-S /Users/roucogar/Library/CloudStorage/Dropbox/Trajectories_Rouco_et_al/Limb_wide_early-late/ChIP/H3K27Ac-FL-E105-Wt-Mm-Rep1-L4288_macs_SR200_norm.bw \
/Users/roucogar/Library/CloudStorage/Dropbox/Trajectories_Rouco_et_al/Limb_wide_early-late/ChIP/H3K27Ac-FL-E135-Wt-Mm-Rep2-L7110_macs_SR200_norm.bw \
--samplesLabel H3K27ac_eFL_E105_GA H3K27ac_eFL_E135_GA \
-R /Users/roucogar/Library/CloudStorage/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/6_Chip-seq_CUT\&RUN/entierFL_Guillame_data/Early_H3K27ac_FL_classified_peaks.bed Common_H3K27ac_FL_classified_peaks.bed Late_H3K27ac_FL_classified_peaks.bed \
--beforeRegionStartLength 3000 \
--afterRegionStartLength 3000 \
-o /Users/roucogar/Library/CloudStorage/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/6_Chip-seq_CUT\&RUN/entierFL_Guillame_data/all_matrix_H3K27ac_eFL_GA_3000.mat.gz


plotHeatmap -m /Users/roucogar/Library/CloudStorage/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/6_Chip-seq_CUT\&RUN/entierFL_Guillame_data/all_matrix_H3K27ac_eFL_GA_3000.mat.gz \
--colorMap Blues \
--regionsLabel  Early Common Late \
--heatmapWidth 10 \
-o /Users/roucogar/Library/CloudStorage/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/6_Chip-seq_CUT\&RUN/entierFL_Guillame_data/all_tornado_H3K27ac_eFL_GA_3000.pdf

conda deactivate
