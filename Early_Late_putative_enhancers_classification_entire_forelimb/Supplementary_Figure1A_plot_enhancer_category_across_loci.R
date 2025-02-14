library(reshape2)
library(ggalluvial)
library(ggplot2)
library(ggrepel)

domains_mm39 <- read.csv("~/Library/CloudStorage/Dropbox/Trajectories_Rouco_et_al/Limb_wide_early-late/domains_mm39.csv")

df2 <- domains_mm39[,c(1,3,5:7)]

df2.wide <- melt(df2, idvar = "name", direction = "wide")

ggplot(df2.wide, aes(fill= variable, x=enh.category, y= value))+ 
  geom_bar(position="stack", stat="identity") + scale_fill_manual(values=c('#1c75bc', '#006838', '#662d91'))+
  labs(x = "Enhancer category distribtuion across loci", y = "Number of enhancers")

#006838 early
#662d91 late
#1c75bc common

