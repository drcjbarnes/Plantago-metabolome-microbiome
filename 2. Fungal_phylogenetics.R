
library(tidyr)
library(tidyverse)
library(reshape2)
library(dplyr)
library(plyr)
library(tibble)
library(taRifx)
library(stringr)

library(vegan)
library(mvabund)
library(class) 
library(lme4)
library(dplyr)
library(MASS)
library(data.table)
library(RColorBrewer)
library(sjmisc)

library(ggplot2); theme_set(theme_bw(base_size = 16))
library(ggrepel)
library(RColorBrewer)
library(corrplot)
library(grid)
library(gridExtra)
library(patchwork)
library(venn)

library(ape)
library(phytools)
library(caper)
library(picante)

library(mclust)
library(cluster)

library(geosphere)
library(mapdata)
library(maps)
library(googleway)
library(ggmap)

library(dada2)
library(funrar) #Converting into relative abundances 
library(metacoder)
library(taxa)
library(stringr)

library(usedist)
library(phyloseq)
library(sjmisc)
library(pracma)
library(FSA)

library(mixOmics)
library(mixKernel)
library(igraph)
library(qgraph)


library(msa)
library(ape)
library(ips)
library(phangorn)
library(DECIPHER)
library("seqinr")
library(phytools)

setwd("/Users/christopherbarnes/Dropbox/Danish_Plantago_microbiomes/Fungi")

rel_data_mean_leaves <- readRDS(file = "rel_data_mean_leaves.RDS")
targeted_chem_mean_leaves <- readRDS(file = "targeted_chem_mean_leaves.RDS")
untargeted_chem_mean_leaves <- readRDS(file = "untargeted_chem_mean_leaves.RDS")
metadata_mean_leaves <- readRDS(file = "metadata_mean_leaves.RDS")

tree <- read.tree(file = "fungal_tree.tree")


a <- which(colnames(rel_data_mean_leaves) %in% tree$tip.label)
tree_rel_data <- rel_data_mean_leaves[ , (a) ]

tips <- c(colnames(tree_rel_data))
tree_trimmed <- keep.tip(tree, tips)

length(tree_trimmed$tip.label)
length(which(colnames(tree_rel_data) %in% tree_trimmed$tip.label))

t_rel_data_reduced <- as.data.frame(t(tree_rel_data))




pd_result <- pd(tree_rel_data, tree_trimmed, include.root = FALSE)
saveRDS(pd_result, file = "pd_results_combined_leaves.RDS")
aggregate(pd_result$PD, by = list(metadata_mean_leaves$Stratum), FUN = mean)

phydist <- cophenetic(tree_trimmed)
comdist.result <- comdist(tree_rel_data, phydist)
saveRDS(comdist.result, file = "comdist_result_combined_leaves.RDS")

comdistnt.data <- comdistnt(tree_rel_data, phydist, 
  abundance.weighted = FALSE, exclude.conspecifics = FALSE)


#ses.mpd.result <- 
# ses.mpd(tree_rel_data, phydist, null.model="taxa.labels", abundance.weighted=FALSE, runs=99)
#saveRDS(ses.mpd.result, file = "ses.mpd.result_combined_leaves.RDS")
ses.mpd.result <- readRDS(file = "ses.mpd.result_combined_leaves.RDS")

aggregate(ses.mpd.result$mpd.obs.z, by = list(metadata_mean_leaves$Stratum), FUN = mean)
#mpd.obs.z is equivalent to -NRI

ses.mntd.result <- 
  ses.mntd(tree_rel_data, phydist, null.model="taxa.labels", abundance.weighted=FALSE, runs=99)
saveRDS(ses.mntd.result, file = "ses.mntd.result_combined_leaves.RDS")

aggregate(ses.mntd.result$mntd.obs.z, by = list(metadata_mean_leaves$Stratum), FUN = mean)
aggregate(ses.mntd.result$mntd.obs.p, by = list(metadata_mean_leaves$Stratum), FUN = mean)

mntd_results <- mntd(tree_rel_data, phydist, abundance.weighted=FALSE)    
aggregate(mntd_results, by = list(metadata_mean_leaves$Stratum), FUN = mean)
saveRDS(mntd_results, file = "mntd_results_combined_leaves.RDS")




ses.mpd.result.weighted <- 
  ses.mpd(tree_rel_data, phydist, null.model="taxa.labels", abundance.weighted = TRUE, runs=99)
saveRDS(ses.mpd.result.weighted, file = "ses.mpd.result.weighted_combined_leaves.RDS")

aggregate(ses.mpd.result.weighted$mpd.obs.z, by = list(metadata_mean_leaves$Stratum), FUN = mean)
#mpd.obs.z is equivalent to -NRI

ses.mntd.result.weighted <- 
  ses.mntd(tree_rel_data, phydist, null.model="taxa.labels", abundance.weighted = TRUE, runs=99)
saveRDS(ses.mntd.result.weighted, file = "ses.mntd.result.weighted_combined_leaves.RDS")

aggregate(ses.mntd.result.weighted$mntd.obs.z, by = list(metadata_mean_leaves$Stratum), FUN = mean)
aggregate(ses.mntd.result.weighted$mntd.obs.p, by = list(metadata_mean_leaves$Stratum), FUN = mean)

mntd_results.weighted <- mntd(tree_rel_data, phydist, abundance.weighted = TRUE)    
aggregate(mntd_results.weighted, by = list(metadata_mean_leaves$Stratum), FUN = mean)
saveRDS(mntd_results.weighted, file = "mntd_results.weighted_combined_leaves.RDS")

mntd_results.weighted <- readRDS(file = "mntd_results.weighted_combined_leaves.RDS")



Fungal_phylo_results <- data.frame("PD" = pd_result$PD, "MPD" = ses.mpd.result$mpd.obs.z, "MNTD" = ses.mntd.result$mntd.obs.z, "Stratum" = metadata_mean_leaves$Stratum_simplified)
saveRDS(Fungal_phylo_results, file = "Fungal_phylo_results.RDS")







