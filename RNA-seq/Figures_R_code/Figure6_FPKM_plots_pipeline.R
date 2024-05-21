library(reshape2)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(usefulLDfunctions)


prefix <- "/Users/roucogar/Library/CloudStorage/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/7_RNA-seq/Final_analysis/DESeq2_251_vs_58//"

# Get all FPKMs
fpkm.table <- read.delim(file.path(prefix, "All_251_vs_58/outputs/mergedTables/AllCufflinks_Simplified_norm.txt"))

# Get samples plan
samples.plan.df <- read.delim(file.path(prefix, "All_251_vs_58/samplesPlan.txt"))

# Select genes
my.genes <- c("Shox2","EYFP","CRE","dsmCherry","Rsrc1","Shh","Irx3","Hoxd13","Msx1","Meis2")

# Extract FPKM of interest
big.df <- melt(
  fpkm.table[fpkm.table$gene_short_name %in% my.genes, ],
  value.name = "FPKM"
)

# Extract sample column
big.df$sample <- gsub("FPKM_", "", big.df$variable)

# Merge with samples plan
big.df <- merge(big.df, samples.plan.df, all = TRUE)

# Set the Control as ref
big.df$Genotype <- factor(big.df$Genotype,
  levels = c("control", setdiff(unique(big.df$Genotype), "control"))
)

# Get log2fc and p-values from DESeq2

# First find the files
deseq2.files <- list.files(
  path = prefix,
  recursive = TRUE,
  pattern = "^DESeq2Analysis_PerGeno.txt$",
  full.names = TRUE
)

# Read them
deseq2.results <- lapply(deseq2.files, read.delim)

# Extract rows and cols of interest
deseq2.results.subset <- lapply(
  deseq2.results,
  subset,
  subset = gene_name %in% my.genes,
  select = c("gene_name", "pvalue", "padj", "log2FoldChange")
)

# Add name of file before merging
deseq2.results.subset.df <- do.call(
  rbind,
  lapply(1:length(deseq2.files), function(i) {
    df <- deseq2.results.subset[[i]]
    df$file <- deseq2.files[[i]]
    return(df)
  })
)

# Deduce metadata from file
deseq2.results.subset.df <- deseq2.results.subset.df %>%
  mutate(
    relPath = gsub("^/+", "", gsub(prefix, "", file)),
    Stage = gsub("/.*", "", relPath),
  )

deseq2.results.subset.df <- merge(
  deseq2.results.subset.df,
  unique(samples.plan.df[, c("Genotype", "Stage")])
)

for (my.gene in my.genes) {
  for (my.name in setdiff(unique(samples.plan.df$Genotype), "control")) {
    print(paste(my.gene, "in", my.name))
    subset.df <- subset(
      big.df,
      gene_short_name == my.gene & Genotype %in% c("control", my.name)
    )

    stat.test <- subset(
      deseq2.results.subset.df,
      Genotype == my.name & gene_name == my.gene
    ) %>%
      rename(group1 = Genotype) %>%
      mutate(
        group2 = "control",
        y.position = max(subset.df$FPKM) * 1.05
      )

    g <- ggplot(
      subset.df,
      aes(x = Genotype, y = FPKM, color = Genotype)
    ) +
      geom_boxplot(outlier.color = NA) +
      # geom_jitter() +
      geom_point(position=position_jitterdodge())+
      facet_grid(. ~ Stage) +
      theme_minimal() +
      xlab("") +
      scale_y_continuous(
        expand = expansion(mult = c(0, .1)),
        limits = c(0, NA)
      ) +
      theme(axis.text.x = element_blank()) +
      scale_color_manual(values = c("black", "red")) +
      stat_pvalue_manual(
        stat.test,
        label = "padj = {format(padj, digits = 2)}\nL2FC = {format(log2FoldChange, digits = 2)}",
        label.size = 2.5
      ) +
      ggtitle(my.gene)
    ggsave(file.path(prefix, "FPKM_plots", paste0(my.gene, "_", my.name, ".pdf")), g, width = 5, height = 5)
  }
}

