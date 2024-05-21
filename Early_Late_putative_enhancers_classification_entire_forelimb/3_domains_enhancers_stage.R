####
# Note that 4 interaction domains have some overlap with each others (Tbx3 and Tbx5, Hoxa13 and Hoxa11). 
# These overlaps contain some H3K27ac peaks (12 and 30 H3K27ac peaks, respectively).
# This is why the number of enhancers/peaks is different than the total number of plotted enhancers in the final boxplot.

# move into R
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(ggplot2)
library(plyranges)
library(ggrepel)
library(reshape2)

wd <- "/Users/nccrgenetics/Library/CloudStorage/Dropbox/Trajectories_Rouco_et_al/Limb_wide_early-late/"
domains <- "Filtered_Curated_gene_and_interaction_domain_mm39.bed"
enhancers <- "ChIP/H3K27ac_FL_classified_peaks.bed"

# domains and stages-classified FL enhancers/H3K27ac peaks
setwd(wd)

domains.gr <- import.bed(domains)
enhancers.gr <- import.bed(enhancers)

ov.enh <- as.data.frame(findOverlaps(domains.gr, enhancers.gr))
ov.enh$cate <- enhancers.gr$name[ov.enh$subjectHits]
ov.enh$genes <- domains.gr$name[ov.enh$queryHits]

# add to the enhancer.gr the name of the interaction domain gene
enhancers.gr$name <- paste0(enhancers.gr$name, "_", ov.enh$genes[match(1:length(enhancers.gr),ov.enh$subjectHits)])
export.bed(enhancers.gr, "H3K27ac_FL_classified_peaks_with_genes_mm39.bed")

ov.enh.summary <- summarise(group_by(ov.enh, queryHits),
                            genes = paste(unique(sort(genes))),
                            enhancers = paste(unique(sort(subjectHits)), collapse = ", "),
                            # categories = paste(unique(sort(cate)), collapse = ", "),
                            categories = paste(unique(sort(cate)), collapse = "_"),
                            nenhancers = length(unique(subjectHits)),
                            ncat = length(unique(cate)),
                            n_C = sum(cate == "C"),
                            n_E = sum(cate == "E"),
                            n_L = sum(cate == "L"))

domains.gr$enh <- NA
domains.gr$enh[ov.enh.summary$queryHits] <- ov.enh.summary$enhancers

domains.gr$n.enh <- 0
domains.gr$n.enh[ov.enh.summary$queryHits] <- ov.enh.summary$nenhancers

domains.gr$enh.category <- "No_enhancer"
domains.gr$enh.category[ov.enh.summary$queryHits] <- ov.enh.summary$categories

domains.gr$n.enh.category <- 0
domains.gr$n.enh.category[ov.enh.summary$queryHits] <- ov.enh.summary$ncat

order.enh.cat = c("C_E_L", "C_E", "C_L", "E_L", "C", "E", "L")
domains.gr$enh.category <- factor(domains.gr$enh.category, levels = order.enh.cat)

domains.gr$n_C_enhancers <- 0
domains.gr$n_C_enhancers <- ov.enh.summary$n_C

domains.gr$n_E_enhancers <- 0
domains.gr$n_E_enhancers <- ov.enh.summary$n_E

domains.gr$n_L_enhancers <- 0
domains.gr$n_L_enhancers <- ov.enh.summary$n_L

# transform domains.gr in a dataframe and reshape the variables with melt to run ggplot (one row per point)
# mcols extract metadata of the gr object into a DataFrame object
frame_domains <- as.data.frame(mcols(domains.gr))
write.table(select(frame_domains, name, n.enh, enh.category, n.enh.category, n_C_enhancers, n_E_enhancers, n_L_enhancers), "domains_mm39.csv", sep=",", quote = FALSE, col.names = TRUE, row.names = FALSE)
