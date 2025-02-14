# Load necessary libraries
library(ggplot2)
library(reshape2)
library(ggpubr)
library(readxl)

#load table
qPCR_58_242_all_genes_of_interest <- read_excel("~/Library/CloudStorage/Dropbox/Andrey_LAB_EXP_DESIGN_RESULTS/3_qPCRs/29_58_242_CD1_RT-qPCR_2024/qPCR_58_CD1_all_genes_of_interest.xlsx")
View(qPCR_58_242_all_genes_of_interest)       

data_matrix <- qPCR_58_242_all_genes_of_interest

# Assume your matrix is named `data_matrix`
# Create a data frame for better manipulation
genes <- data_matrix$Gene  # Store gene names if not already
samples <- colnames(data_matrix)  # Store sample names if not already
data_long <- melt(data_matrix, varnames = c("Gene", "Sample"), value.name = "Expression")

data_long$Group <- gsub("[0-9]$", "", data_long$variable)

#plot  with statistics 

p2 <- ggplot(data_long, aes(x = Group, y = Expression, fill = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.5) + 
  geom_jitter(width = 0.2, size = 1.5) +  # Jitter points for replicates
  facet_wrap(~Gene, scales = "free") +  # Separate plots for each gene
  labs(title = "RT-qPCR Expression Data",
       x = "Sample Group", y = "Relative Expression") +
  theme_minimal() + scale_y_continuous(expand = c(.004, .004)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  stat_summary(fun=mean, geom="crossbar", linetype = "dashed", linewidth = 0.2, width=0.5) 

my_comparisons1 <- list(c("E135-#242-FL","E135-CD1-FL"))

p2 + stat_compare_means(comparisons = my_comparisons1, method = "t.test") +
  ggtitle("RT-qPCR relative expression 58 & CD1, t-test")

ggsave("RT-qPCR relative expression 58 & CD1 t-test.pdf", width = 12, height = 8)
