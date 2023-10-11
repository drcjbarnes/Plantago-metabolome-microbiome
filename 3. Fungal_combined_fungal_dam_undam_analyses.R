
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


setwd("/Users/christopherbarnes/Dropbox/Danish_Plantago_microbiomes/Fungi")

rel_data <- readRDS(file = "final_data.RDS")
metadata <- readRDS(file = "final_metadata.RDS")
untargeted_chem <- readRDS(file = "final_untargeted_chem.RDS")
targeted_chem <- readRDS(file = "final_targeted_chem.RDS")
Taxonomy <- readRDS(file = "Fungal_taxonomy")
dim(rel_data)

metadata$Stratum_simplified <- str_replace(metadata$Stratum, "Damaged leaf", "Leaf")
metadata$Stratum_simplified <- str_replace(metadata$Stratum_simplified, "Undamaged leaf", "Leaf")
metadata$merging_var <- paste(metadata$PlantID, metadata$Stratum_simplified, sep = "_")

rel_data_mean_leaves <- aggregate(rel_data, by = list(Category = metadata$merging_var), FUN = mean)
rownames(rel_data_mean_leaves) <- rel_data_mean_leaves[ , 1]
rel_data_mean_leaves <- rel_data_mean_leaves[ , -c(1)]
colnames(rel_data_mean_leaves) <- colnames(rel_data)


metadata_mean_leaves <- metadata[!duplicated(metadata$merging_var),]
dim(metadata_mean_leaves)

metadata_mean_leaves <- metadata_mean_leaves[order(metadata_mean_leaves$merging_var) , ]
dim(metadata_mean_leaves)
rel_data_mean_leaves <- rel_data_mean_leaves[order(rownames(rel_data_mean_leaves)),]
dim(rel_data_mean_leaves)

dim(untargeted_chem)
dim(targeted_chem)
dim(metadata)

untargeted_chem_mean_leaves <- aggregate(untargeted_chem, by = list(Category = metadata$merging_var), FUN = mean)
rownames(untargeted_chem_mean_leaves) <- untargeted_chem_mean_leaves[ , 1]
untargeted_chem_mean_leaves <- untargeted_chem_mean_leaves[ , -c(1)]
colnames(untargeted_chem_mean_leaves) <- colnames(untargeted_chem)
untargeted_chem_mean_leaves <- untargeted_chem_mean_leaves[order(rownames(untargeted_chem_mean_leaves)),]
untargeted_chem_mean_leaves[is.na(untargeted_chem_mean_leaves)] <- 0

targeted_chem_mean_leaves <- aggregate(targeted_chem, by = list(Category = metadata$merging_var), FUN = mean)
rownames(targeted_chem_mean_leaves) <- targeted_chem_mean_leaves[ , 1]
targeted_chem_mean_leaves <- targeted_chem_mean_leaves[ , -c(1)]
targeted_chem_mean_leaves <- targeted_chem_mean_leaves[order(rownames(targeted_chem_mean_leaves)),]
targeted_chem_mean_leaves[is.na(targeted_chem_mean_leaves)] <- 0

#saveRDS(rel_data_mean_leaves, file = "rel_data_mean_leaves.RDS")
#saveRDS(targeted_chem_mean_leaves, file = "targeted_chem_mean_leaves.RDS")
#saveRDS(untargeted_chem_mean_leaves, file = "untargeted_chem_mean_leaves.RDS")
#saveRDS(metadata_mean_leaves, file = "metadata_mean_leaves.RDS")


rownames(rel_data_mean_leaves)
rownames(targeted_chem_mean_leaves)
rownames(untargeted_chem_mean_leaves)
metadata_mean_leaves$merging_var








#####################Leaf data analysis#######################

leaf_data <- subset(rel_data_mean_leaves, metadata_mean_leaves$Stratum_simplified == "Leaf")
leaf_metadata <- subset(metadata_mean_leaves, metadata_mean_leaves$Stratum_simplified == "Leaf")
leaf_untargeted_chem <- subset(untargeted_chem_mean_leaves, metadata_mean_leaves$Stratum_simplified == "Leaf")
leaf_targeted_chem <- subset(targeted_chem_mean_leaves, metadata_mean_leaves$Stratum_simplified == "Leaf")

leaf_data <- with(leaf_data, leaf_data[order(leaf_metadata$PlantID),])
leaf_untargeted_chem <- with(leaf_untargeted_chem, leaf_untargeted_chem[order(leaf_metadata$PlantID),])
leaf_targeted_chem <- with(leaf_targeted_chem, leaf_targeted_chem[order(leaf_metadata$PlantID),])
leaf_metadata <- with(leaf_metadata, leaf_metadata[order(leaf_metadata$PlantID),])

#adonis2(leaf_data ~ leaf_metadata$Site)
adonis2(leaf_data ~ leaf_metadata$Habitat)
#adonis2(leaf_data ~ leaf_metadata$Light)

adonis2(leaf_data ~ leaf_metadata$AnnPrec)
adonis2(leaf_data ~ leaf_metadata$ColdestQuartMeanTemp)
adonis2(leaf_data ~ leaf_metadata$MeanAnnTemp)

adonis2(leaf_data ~ leaf_metadata$pH)
adonis2(leaf_data ~ leaf_metadata$EC)
adonis2(leaf_data ~ leaf_metadata$BulkDensity)

#adonis2(leaf_data ~ leaf_metadata$PCNM1)
#adonis2(leaf_data ~ leaf_metadata$PCNM2)
adonis2(leaf_data ~ leaf_metadata$Island)
adonis2(leaf_data ~ leaf_metadata$Geography3)

adonis2(leaf_data ~ leaf_metadata$Habitat + leaf_metadata$Island +
          leaf_metadata$ColdestQuartMeanTemp, 
        group = leaf_metadata$Site)


kruskal.test(leaf_metadata$Richness ~ leaf_metadata$Site)
kruskal.test(leaf_metadata$Richness ~ leaf_metadata$Habitat)
kruskal.test(leaf_metadata$Richness ~ leaf_metadata$Light)

cor.test(leaf_metadata$Richness, leaf_metadata$AnnPrec, method = "spearman")
cor.test(leaf_metadata$Richness, leaf_metadata$ColdestQuartMeanTemp, method = "spearman")
cor.test(leaf_metadata$Richness, leaf_metadata$MeanAnnTemp, method = "spearman")

cor.test(leaf_metadata$Richness, leaf_metadata$pH, method = "spearman")
cor.test(leaf_metadata$Richness, leaf_metadata$EC, method = "spearman")
cor.test(leaf_metadata$Richness, leaf_metadata$BulkDensity, method = "spearman")

cor.test(leaf_metadata$Richness, leaf_metadata$PCNM1, method = "spearman")
cor.test(leaf_metadata$Richness, leaf_metadata$PCNM2, method = "spearman")
kruskal.test(leaf_metadata$Richness ~ leaf_metadata$Island)
kruskal.test(leaf_metadata$Richness ~ leaf_metadata$Geography3)

leaf_bray_data <- vegdist(leaf_data, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
leaf_untargeted_bray_data <- vegdist(leaf_untargeted_chem, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
#leaf_untargeted_bray_data <- vegdist((leaf_untargeted_chem - (min(leaf_untargeted_chem))), Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
leaf_targeted_bray_data <- vegdist(leaf_targeted_chem, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
#leaf_targeted_bray_data <- vegdist((leaf_targeted_chem - (min(leaf_targeted_chem))), Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence

mantel(leaf_bray_data, leaf_untargeted_bray_data)
mantel(leaf_bray_data, leaf_targeted_bray_data)

leaf_bray_vector <- as.vector(vegdist(leaf_data, Type = "bray", binary = T))
label = rep("Damaged leaf", times = length(leaf_bray_vector))
leaf_sim_data <- data.frame("Similarity" = leaf_bray_vector, label)


adonis2(leaf_bray_data ~ leaf_targeted_chem$GeniposideHeight)
adonis2(leaf_bray_data ~ leaf_targeted_chem$GardosideHeight)
adonis2(leaf_bray_data ~ leaf_targeted_chem$MelittosideHeight)
adonis2(leaf_bray_data ~ leaf_targeted_chem$X10.acetoxymajorHeight)
adonis2(leaf_bray_data ~ leaf_targeted_chem$VerbascosideHeight)
adonis2(leaf_bray_data ~ leaf_targeted_chem$AucubinHeight)

leaf_MDS_data <- monoMDS(leaf_bray_data) # Performs multidimensional scaling for the data
plot(leaf_MDS_data)

Nmds1 <- leaf_MDS_data$points[,1] #Extracts dimension 1 as a vector
Nmds2 <- leaf_MDS_data$points[,2] #Extracts dimension 2 as a vector

Nmds <-cbind(Nmds1, Nmds2) #Makes a new sheet with dimension1 and 2
Nmds <- as.data.frame(Nmds) #Creates a new dataframe with dimensions 1 and 2

ggplot(data = Nmds, aes(y = Nmds2, x = Nmds1, Type="p", 
  colour = leaf_metadata$Habitat)) + 
  geom_point(size = 4, alpha = 0.5) + labs(x="Dimension 1", y = "Dimension 2", colour = NULL, shape = NULL) + 
  theme(plot.title=element_text(hjust=0)) + 
  ggtitle("D") +
  theme(legend.position = "right") +
  scale_colour_brewer(palette = "Set1") 


a <- prcomp(leaf_untargeted_chem)
summary(a)
plot(a)
leaf_untar_chem_PC1 <- a$x[,1]
leaf_untar_chem_PC2 <- a$x[,2]
leaf_untar_chem_PC3 <- a$x[,3]
leaf_untar_chem_PC4 <- a$x[,4]

adonis2(leaf_bray_data ~ leaf_untar_chem_PC1)
adonis2(leaf_bray_data ~ leaf_untar_chem_PC2)
adonis2(leaf_bray_data ~ leaf_untar_chem_PC3)
adonis2(leaf_bray_data ~ leaf_untar_chem_PC4)















########################Cross strata comparisons########################
###################Separating the other datasets########################
soil_data <- subset(rel_data_mean_leaves, metadata_mean_leaves$Stratum_simplified == "Soil")
soil_metadata <- subset(metadata_mean_leaves, metadata_mean_leaves$Stratum_simplified == "Soil")
soil_untargeted_chem <- subset(untargeted_chem_mean_leaves, metadata_mean_leaves$Stratum_simplified == "Soil")
soil_targeted_chem <- subset(targeted_chem_mean_leaves, metadata_mean_leaves$Stratum_simplified == "Soil")

root_data <- subset(rel_data_mean_leaves, metadata_mean_leaves$Stratum_simplified == "Root")
root_metadata <- subset(metadata_mean_leaves, metadata_mean_leaves$Stratum_simplified == "Root")
root_untargeted_chem <- subset(untargeted_chem_mean_leaves, metadata_mean_leaves$Stratum_simplified == "Root")
root_targeted_chem <- subset(targeted_chem_mean_leaves, metadata_mean_leaves$Stratum_simplified == "Root")




#soil-root comparisons
a <- which(table(c(soil_metadata$PlantID, root_metadata$PlantID)) > 1)
a <- names(a)

root_data_soil_root <- subset(root_data, subset = root_metadata$PlantID %in% c(a))
dim(root_data_soil_root)
root_metadata_soil_root <- subset(root_metadata, subset = root_metadata$PlantID %in% c(a))
dim(root_metadata_soil_root)

soil_data_soil_root <- subset(soil_data, subset = soil_metadata$PlantID %in% c(a))
dim(soil_data_soil_root)
soil_metadata_soil_root <- subset(soil_metadata, subset = soil_metadata$PlantID %in% c(a))
dim(soil_metadata_soil_root)

soil_data_soil_root_bray <- vegdist(soil_data_soil_root, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
root_data_soil_root_bray <- vegdist(root_data_soil_root, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence

mantel(soil_data_soil_root_bray, root_data_soil_root_bray)
cor.test(root_metadata_soil_root$Richness, soil_metadata_soil_root$Richness, method = "spearman")


root_data_soil_root_binary <- (root_data_soil_root > 0) # logical, or
root_data_soil_root_binary <- (root_data_soil_root_binary > 0)*1L # integer 01
root_data_soil_root_binary <- as.data.frame(root_data_soil_root_binary)

soil_data_soil_root_binary <- (soil_data_soil_root > 0) # logical, or
soil_data_soil_root_binary <- (soil_data_soil_root_binary > 0)*1L # integer 01
soil_data_soil_root_binary <- as.data.frame(soil_data_soil_root_binary)

both_data_soil_root_binary <- list()
just_soil_data_soil_root_binary <- list()
just_root_data_soil_root_binary <- list()

for (i in 1:nrow(soil_data_soil_root_binary)){
  both_data_soil_root_binary[[i]] <- length(which(root_data_soil_root_binary[i,] == 1 & soil_data_soil_root_binary[i,] == 1))
  just_root_data_soil_root_binary[[i]] <- length(which(root_data_soil_root_binary[i,] == 1 & soil_data_soil_root_binary[i,] == 0))
  just_soil_data_soil_root_binary[[i]] <- length(which(root_data_soil_root_binary[i,] == 0 & soil_data_soil_root_binary[i,] == 1))
}

both_data_soil_root_binary <- unlist(both_data_soil_root_binary, use.names = FALSE)
just_root_data_soil_root_binary <- unlist(just_root_data_soil_root_binary, use.names = FALSE)
just_soil_data_soil_root_binary <- unlist(just_soil_data_soil_root_binary, use.names = FALSE)
soil_root_binary_comparisons <- data.frame(both_data_soil_root_binary, just_root_data_soil_root_binary, just_soil_data_soil_root_binary)
soil_root_binary_comparisons

soil_root_binary_comparisons_means <- colMeans(soil_root_binary_comparisons)





#soil-undamaged comparisons
a <- which(table(c(soil_metadata$PlantID, leaf_metadata$PlantID)) > 1)
a <- names(a)

soil_data_leaf_soil <- subset(soil_data, subset = soil_metadata$PlantID %in% c(a))
dim(soil_data_leaf_soil)
soil_metadata_leaf_soil <- subset(soil_metadata, subset = soil_metadata$PlantID %in% c(a))
dim(soil_metadata_leaf_soil)

leaf_data_leaf_soil <- subset(leaf_data, subset = leaf_metadata$PlantID %in% c(a))
dim(leaf_data_leaf_soil)
leaf_metadata_leaf_soil <- subset(leaf_metadata, subset = leaf_metadata$PlantID %in% c(a))
dim(leaf_metadata_leaf_soil)

leaf_targeted_chem_leaf_soil <- subset(leaf_targeted_chem, subset = leaf_metadata$PlantID %in% c(a))
dim(leaf_targeted_chem_leaf_soil)
leaf_untargeted_chem_leaf_soil <- subset(leaf_untargeted_chem, subset = leaf_metadata$PlantID %in% c(a))
dim(leaf_untargeted_chem_leaf_soil)

leaf_data_leaf_soil_bray <- vegdist(leaf_data_leaf_soil, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
soil_data_leaf_soil_bray <- vegdist(soil_data_leaf_soil, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
leaf_targeted_chem_leaf_soil_bary <- vegdist(leaf_targeted_chem_leaf_soil, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
leaf_untargeted_chem_leaf_soil_bary <- vegdist(leaf_untargeted_chem_leaf_soil, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence

mantel(leaf_data_leaf_soil_bray, soil_data_leaf_soil_bray)
cor.test(soil_metadata_leaf_soil$Richness, leaf_metadata_leaf_soil$Richness, method = "spearman")

mantel(soil_data_leaf_soil_bray, leaf_targeted_chem_leaf_soil_bary)
mantel(soil_data_leaf_soil_bray, leaf_untargeted_chem_leaf_soil_bary)

adonis2(soil_data_leaf_soil_bray ~ leaf_targeted_chem_leaf_soil$GeniposideHeight)
adonis2(soil_data_leaf_soil_bray ~ leaf_targeted_chem_leaf_soil$GardosideHeight)
adonis2(soil_data_leaf_soil_bray ~ leaf_targeted_chem_leaf_soil$MelittosideHeight)
adonis2(soil_data_leaf_soil_bray ~ leaf_targeted_chem_leaf_soil$X10.acetoxymajorHeight)
adonis2(soil_data_leaf_soil_bray ~ leaf_targeted_chem_leaf_soil$VerbascosideHeight)
adonis2(soil_data_leaf_soil_bray ~ leaf_targeted_chem_leaf_soil$AucubinHeight)


leaf_data_leaf_soil_binary <- (leaf_data_leaf_soil > 0) # logical, or
leaf_data_leaf_soil_binary <- (leaf_data_leaf_soil_binary > 0)*1L # integer 01
leaf_data_leaf_soil_binary <- as.data.frame(leaf_data_leaf_soil_binary)

soil_data_leaf_soil_binary <- (soil_data_leaf_soil > 0) # logical, or
soil_data_leaf_soil_binary <- (soil_data_leaf_soil_binary > 0)*1L # integer 01
soil_data_leaf_soil_binary <- as.data.frame(soil_data_leaf_soil_binary)

both_leaf_soil_binary <- list()
just_soil_data_leaf_soil_binary <- list()
just_leaf_data_leaf_soil_binary <- list()

for (i in 1:nrow(soil_data_leaf_soil_binary)){
  both_leaf_soil_binary[[i]] <- length(which(leaf_data_leaf_soil_binary[i,] == 1 & soil_data_leaf_soil_binary[i,] == 1))
  just_leaf_data_leaf_soil_binary[[i]] <- length(which(leaf_data_leaf_soil_binary[i,] == 1 & soil_data_leaf_soil_binary[i,] == 0))
  just_soil_data_leaf_soil_binary[[i]] <- length(which(leaf_data_leaf_soil_binary[i,] == 0 & soil_data_leaf_soil_binary[i,] == 1))
}

both_leaf_soil_binary <- unlist(both_leaf_soil_binary, use.names = FALSE)
just_leaf_data_leaf_soil_binary <- unlist(just_leaf_data_leaf_soil_binary, use.names = FALSE)
just_soil_data_leaf_soil_binary <- unlist(just_soil_data_leaf_soil_binary, use.names = FALSE)
leaf_soil_binary_comparisons <- data.frame(both_leaf_soil_binary, just_leaf_data_leaf_soil_binary, just_soil_data_leaf_soil_binary)
leaf_soil_binary_comparisons

leaf_soil_binary_comparisons_means <- colMeans(leaf_soil_binary_comparisons)





#root-undamaged comparisons
a <- which(table(c(root_metadata$PlantID, leaf_metadata$PlantID)) > 1)
a <- names(a)

root_data_leaf_soil <- subset(root_data, subset = root_metadata$PlantID %in% c(a))
dim(root_data_leaf_soil)
root_metadata_leaf_soil <- subset(root_metadata, subset = root_metadata$PlantID %in% c(a))
dim(root_metadata_leaf_soil)

leaf_data_leaf_soil <- subset(leaf_data, subset = leaf_metadata$PlantID %in% c(a))
dim(leaf_data_leaf_soil)
leaf_metadata_leaf_soil <- subset(leaf_metadata, subset = leaf_metadata$PlantID %in% c(a))
dim(leaf_metadata_leaf_soil)

leaf_targeted_chem_leaf_soil <- subset(leaf_targeted_chem, subset = leaf_metadata$PlantID %in% c(a))
dim(leaf_targeted_chem_leaf_soil)
leaf_untargeted_chem_leaf_soil <- subset(leaf_untargeted_chem, subset = leaf_metadata$PlantID %in% c(a))
dim(leaf_untargeted_chem_leaf_soil)

leaf_data_leaf_root_bray <- vegdist(leaf_data_leaf_soil, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
root_data_leaf_root_bray <- vegdist(root_data_leaf_soil, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
leaf_targeted_chem_leaf_root_bary <- vegdist(leaf_targeted_chem_leaf_soil, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
leaf_untargeted_chem_leaf_root_bary <- vegdist(leaf_untargeted_chem_leaf_soil, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence

mantel(leaf_data_leaf_root_bray, root_data_leaf_root_bray)
cor.test(root_metadata_leaf_soil$Richness, leaf_metadata_leaf_soil$Richness, method = "spearman")

mantel(root_data_leaf_root_bray, leaf_targeted_chem_leaf_root_bary)
mantel(root_data_leaf_root_bray, leaf_untargeted_chem_leaf_root_bary)

adonis2(root_data_leaf_root_bray ~ leaf_targeted_chem_leaf_soil$GeniposideHeight)
adonis2(root_data_leaf_root_bray ~ leaf_targeted_chem_leaf_soil$GardosideHeight)
adonis2(root_data_leaf_root_bray ~ leaf_targeted_chem_leaf_soil$MelittosideHeight)
adonis2(root_data_leaf_root_bray ~ leaf_targeted_chem_leaf_soil$X10.acetoxymajorHeight)
adonis2(root_data_leaf_root_bray ~ leaf_targeted_chem_leaf_soil$VerbascosideHeight)
adonis2(root_data_leaf_root_bray ~ leaf_targeted_chem_leaf_soil$AucubinHeight)


leaf_data_leaf_root_binary <- (leaf_data_leaf_soil > 0) # logical, or
leaf_data_leaf_root_binary <- (leaf_data_leaf_root_binary > 0)*1L # integer 01
leaf_data_leaf_root_binary <- as.data.frame(leaf_data_leaf_root_binary)

root_data_leaf_root_binary <- (root_data_leaf_soil > 0) # logical, or
root_data_leaf_root_binary <- (root_data_leaf_root_binary > 0)*1L # integer 01
root_data_leaf_root_binary <- as.data.frame(root_data_leaf_root_binary)

both_leaf_root_binary <- list()
just_root_data_leaf_root_binary <- list()
just_leaf_data_leaf_root_binary <- list()

for (i in 1:nrow(root_data_leaf_root_binary)){
  both_leaf_root_binary[[i]] <- length(which(leaf_data_leaf_root_binary[i,] == 1 & root_data_leaf_root_binary[i,] == 1))
  just_leaf_data_leaf_root_binary[[i]] <- length(which(leaf_data_leaf_root_binary[i,] == 1 & root_data_leaf_root_binary[i,] == 0))
  just_root_data_leaf_root_binary[[i]] <- length(which(leaf_data_leaf_root_binary[i,] == 0 & root_data_leaf_root_binary[i,] == 1))
}

both_leaf_root_binary <- unlist(both_leaf_root_binary, use.names = FALSE)
just_leaf_data_leaf_root_binary <- unlist(just_leaf_data_leaf_root_binary, use.names = FALSE)
just_root_data_leaf_root_binary <- unlist(just_root_data_leaf_root_binary, use.names = FALSE)
leaf_root_binary_comparisons <- data.frame(both_leaf_root_binary, just_leaf_data_leaf_root_binary, just_root_data_leaf_root_binary)
leaf_root_binary_comparisons

leaf_root_binary_comparisons_means <- colMeans(leaf_root_binary_comparisons)
















################Defining Jaccard Index calculations################
Jaccard = function (x, y) {
  M.11 = sum(x == 1 & y == 1)
  M.10 = sum(x == 1 & y == 0)
  M.01 = sum(x == 0 & y == 1)
  return (M.11 / (M.11 + M.10 + M.01))
}

#################Calculating Jaccard coefficients##################

binary_data <- (rel_data > 0) # logical, or
binary_data <- (binary_data > 0)*1L # integer 01
t_binary_data <- t(binary_data)
input.variables = data.frame(t_binary_data)

m = matrix(data = NA, nrow = length(input.variables), ncol = length(input.variables))
for (r in 1:length(input.variables)) {
  for (c in 1:length(input.variables)) {
    if (c == r) {
      m[r,c] = 1
    } else if (c > r) {
      m[r,c] = Jaccard(input.variables[,r], input.variables[,c])
    }
  }
}

variable.names = sapply(input.variables, attr, "label")
colnames(m) = variable.names
rownames(m) = variable.names   
jaccards = m
jaccards
dim(jaccards)


makeSymm <- function(m) {
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  return(m)
}

jaccards <- as.data.frame(makeSymm(jaccards))
jaccards <- as.data.frame(makeSymm(jaccards))
jaccards[is.na(jaccards)] <- 0

row.names(jaccards) <- row.names(rel_data)
colnames(jaccards) <- row.names(rel_data)

#saveRDS(jaccards, "JaccardsRDS")
jaccards <- readRDS(file = "JaccardsRDS")



all_bray_data <- as.data.frame(makeSymm(as.matrix(all_bray_data)))
all_bray_data <- as.data.frame(makeSymm(all_bray_data))
all_bray_data[is.na(all_bray_data)] <- 0
dim(jaccards)

long_status <- rep(metadata$Stratum, times = nrow(jaccards))
long_PlantID <- rep(metadata$PlantID, times = nrow(jaccards))
long_jaccards <- as.vector(jaccards)

long_data <- data.frame("Stratum" = long_status, "PlantID" = long_PlantID, "jaccard" = long_jaccards)
long_data$Stratum_simplified <- str_replace(long_data$Stratum, "Damaged leaf", "Leaf")
long_data$Stratum_simplified <- str_replace(long_data$Stratum_simplified, "Undamaged leaf", "Leaf")

jaccard_Stratums <- dist_groups(jaccards, metadata$Stratum)
jaccard_PlantIDs <- dist_groups(jaccards, metadata$PlantID)
jaccard_Stratum_simplified <- dist_groups(jaccards, metadata$Stratum_simplified)
jaccard_Sites <- dist_groups(jaccards, metadata$Site)

Similarity_Stratums <- dist_groups(all_bray_data, metadata$Stratum)
Similarity_PlantIDs <- dist_groups(all_bray_data, metadata$PlantID)
Similarity_Sites <- dist_groups(all_bray_data, metadata$Site)


ECDifferences <- outer(metadata$EC, metadata$EC, '-')
ECDifferences <- abs(ECDifferences)
dim(ECDifferences)

ECDifferences_Stratums <- dist_groups(ECDifferences, metadata$Stratum)
ECDifferences_PlantIDs <- dist_groups(ECDifferences, metadata$PlantID)
ECDifferences_Sites <- dist_groups(ECDifferences, metadata$Site)

pHDifferences <- outer(metadata$pH, metadata$pH, '-')
pHDifferences <- abs(pHDifferences)
dim(pHDifferences)

pHDifferences_Stratums <- dist_groups(pHDifferences, metadata$Stratum)
pHDifferences_PlantIDs <- dist_groups(pHDifferences, metadata$PlantID)
pHDifferences_Sites <- dist_groups(pHDifferences, metadata$Site)


Big.data <- data.frame("Jaccards" = jaccard_Stratums$Distance, "Similarity" = Similarity_Stratums$Distance,
                       "Stratum" = jaccard_Stratums$Label,
                       "Stratum_simplified" = jaccard_Stratum_simplified$Label,
                       "Site" = jaccard_Sites$Label, "PlantID" = jaccard_PlantIDs$Label, 
                       "EC" = ECDifferences_Sites$Distance,
                       "pH" = pHDifferences_Sites$Distance)


#######Within sites excluding self#########
Big_data_within_plant <- Big.data[grep("^Within", Big.data$PlantID),]
a1 <- aggregate(Big_data_within_plant$Jaccards, by = list(Category = Big_data_within_plant$Stratum), FUN = mean)
a2 <- data.frame("Category" = rep("NA", times = 4), "x" = rep("NA", times = 4))
a <- rbind (a1, a2)
colnames(a) <- c("Comparison", "Within plants")

Big_data_no_within_plant <- Big.data[grep("^Within", Big.data$Site),]
Big_data_no_within_plant <- Big_data_no_within_plant[grep("^Between", Big_data_no_within_plant$PlantID),]
b <- aggregate(Big_data_no_within_plant$Jaccards, by = list(Category = Big_data_no_within_plant$Stratum), FUN = mean)

Big_data_between_sites <- Big.data[grep("^Between", Big.data$Site),]
c <- aggregate(Big_data_between_sites$Jaccards, by = list(Category = Big_data_between_sites$Stratum), FUN = mean)
d <- data.frame(a, "Within sites" = b$x, "Between sites" = c$x)
d$Comparison <- c$Category

Big_data_within_plant$Group <- rep("Within plants", nrow(Big_data_within_plant))
Big_data_no_within_plant$Group <- rep("Within sites", nrow(Big_data_no_within_plant))
Big_data_between_sites$Group <- rep("Between sites", nrow(Big_data_between_sites))

Big_data_within_plant$Group2 <- rep("Within stratum", nrow(Big_data_within_plant))
Big_data_no_within_plant$Group2 <- rep("Between strata", nrow(Big_data_no_within_plant))
Big_data_between_sites$Group2 <- rep("Between strata", nrow(Big_data_between_sites))

long.order.data <- melt(d, id.vars = c("Comparison"), variable.name = "Group", value.name = "Jaccards")
long.order.data$Group <- factor(long.order.data$Group, 
                                levels = c("Within.plants", "Within.sites", "Between.sites"),
                                labels = c("Within.plants", "Within.sites", "Between.sites"))
long.order.data$Jaccards <- as.numeric(as.character(long.order.data$Jaccards))



Big_data <- rbind(Big_data_within_plant, Big_data_no_within_plant, Big_data_between_sites)
Big_data$Group <- factor(Big_data$Group, 
                         levels = c("Within plants", "Within sites", "Between sites"),
                         labels = c("Within plants", "Within sites", "Between sites"))

Big_data$Group2 <- ifelse(str_detect(as.vector(Big_data$Stratum),"^Within"),
                          "Within strata","Between strata")

Big_data$Group2 <- factor(Big_data$Group2, 
                          levels = c("Within strata","Between strata"),
                          labels = c("Within strata","Between strata"))

Big_data$Similarity <- (1 - Big_data$Similarity)

Big_data$Barrier <- Big_data$Stratum

Big_data$Barrier <- fct_recode(Big_data$Barrier, 
                               Within_tissue = "Within Undamaged leaf",
                               Within_tissue = "Within Damaged leaf",
                               Within_tissue =  "Within Root",
                               Within_tissue =  "Within Soil",
                               Aboveground_transfer = "Between Undamaged leaf and Damaged leaf",
                               Above_belowground_transfer = "Between Undamaged leaf and Root",
                               Above_belowground_transfer = "Between Undamaged leaf and Soil",
                               Above_belowground_transfer = "Between Damaged leaf and Root",
                               Above_belowground_transfer = "Between Damaged leaf and Soil",
                               Belowground_transfer = "Between Root and Soil")

Big_data$Barrier <- factor(Big_data$Barrier, 
                           levels = c("Within_tissue", "Aboveground_transfer", 
                                      "Belowground_transfer", "Above_belowground_transfer"), 
                           labels = c("Within tissue", "Aboveground transfer", 
                                      "Belowground transfer", "Above-belowground transfer"))


ggplot(subset(Big_data, Big_data$Barrier == "Above-belowground transfer"), 
       aes(x = Group, y = Jaccards, fill = Stratum_simplified))  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_boxplot(alpha = 0.8) + labs(x = NULL, y = "Shared ASVs", fill = NULL) + 
  scale_fill_brewer(palette = "Accent") #+ facet_grid(metadata$Plant ~ metadata$Site) 

saveRDS(Big_data, file = "Fungal_big_data.RDS")









dim(Taxonomy)
dim(leaf_data)

###########################Network analyses#################################

colnames(leaf_data) <- gsub("SeqVar", "fOTU", colnames(rel_data_mean_leaves))
genus_descriptors <- paste(Taxonomy$Genus, colnames(leaf_data), sep = "_")
colnames(leaf_data) <- genus_descriptors
rownames(Taxonomy) <- colnames(leaf_data)

leaf_data_filtered <- leaf_data[,c(which(colSums(leaf_data > 0) > 4))]
x <- nearZeroVar(leaf_untargeted_chem)
leaf_untargeted_chem_filtered <- leaf_untargeted_chem[,-c(x$Position)]

x <- nearZeroVar(leaf_untargeted_chem)
leaf_untargeted_chem_filtered <- leaf_untargeted_chem[,-c(x$Position)]
dim(leaf_untargeted_chem_filtered)

rownames(leaf_untargeted_chem_filtered) <- rownames(leaf_data_filtered)
rownames(leaf_targeted_chem) <- rownames(leaf_data_filtered)

imgCor(leaf_data_filtered, leaf_untargeted_chem_filtered, sideColors = c("purple", "green")) 
imgCor(leaf_data_filtered, leaf_targeted_chem, sideColors = c("purple", "green")) 

# set grid search values for each regularisation parameter
grid1 <- seq(0.001, 0.2, length = 10) 
grid2 <- seq(0.001, 0.2, length = 10)


#leaf_data_untar_rcc_tune <- tune.rcc(leaf_data_filtered, leaf_untargeted_chem_filtered, grid1 = grid1, grid2 = grid2, validation = "loo") 
#leaf_data_untar_rcc_tune <- saveRDS(leaf_data_untar_rcc_tune, file = "leaf_data_untar_rcc_tune")
#leaf_data_untar_rcc_tune <- rcc(leaf_data_filtered, leaf_untargeted_chem_filtered, method = 'shrinkage') 
leaf_data_untar_rcc_tune <- readRDS(file = "leaf_data_untar_rcc_tune")

#leaf_data_tar_rcc_tune <- tune.rcc(leaf_data_filtered, leaf_targeted_chem, grid1 = grid1, grid2 = grid2, validation = "loo") 
#leaf_data_tar_rcc_tune <- saveRDS(leaf_data_tar_rcc_tune, file = "leaf_data_tar_rcc_tune")
leaf_data_tar_rcc_tune <- readRDS(file = "leaf_data_tar_rcc_tune")



# formed optimised CV rCCA. Need to adjust for this for all four comparisons. 
leaf_data_untar_rcc_tune_opt.l1 <- leaf_data_untar_rcc_tune$opt.lambda1 # extract the optimal lambda values
leaf_data_untar_rcc_tune_opt.l2 <- leaf_data_untar_rcc_tune$opt.lambda2

leaf_data_untar_CV_rcc_tune <- rcc(leaf_data_filtered, leaf_untargeted_chem_filtered, method = "ridge", lambda1 = leaf_data_untar_rcc_tune_opt.l1, lambda2 = leaf_data_untar_rcc_tune_opt.l2) 

colnames(leaf_data_untar_CV_rcc_tune$X) <- colnames(leaf_data_filtered)
rownames(leaf_data_untar_CV_rcc_tune$loadings$X) <- colnames(leaf_data_filtered)
leaf_data_untar_CV_rcc_tune$names$colnames$X <- colnames(leaf_data_filtered)

colnames(leaf_data_untar_CV_rcc_tune$Y) <- paste("Untar_metab", seq(1: length(colnames(leaf_data_untar_CV_rcc_tune$Y))), sep = "_")
rownames(leaf_data_untar_CV_rcc_tune$loadings$Y) <- paste("Untar_metab", seq(1: length(colnames(leaf_data_untar_CV_rcc_tune$Y))), sep = "_")
leaf_data_untar_CV_rcc_tune$names$colnames$Y <- paste("Untar_metab", seq(1: length(colnames(leaf_data_untar_CV_rcc_tune$Y))), sep = "_")
plot(leaf_data_untar_CV_rcc_tune, type = "barplot", main = "Cross Validation") 

plotIndiv(leaf_data_untar_CV_rcc_tune, comp = 1:2, 
          ind.names = leaf_metadata$Habitat,
          group = leaf_metadata$Habitat, rep.space = "XY-variate", 
          legend = TRUE, title = '(a) Nutrimouse, rCCA CV XY-space')


plotArrow(leaf_data_untar_CV_rcc_tune, group = leaf_metadata$Habitat, 
          col.per.group = color.mixo(1:5),
          title = '(a) Nutrimouse, CV method')


plotVar(leaf_data_untar_CV_rcc_tune, var.names = c(TRUE, TRUE),
        cex = c(4, 4), cutoff = 0.4,
        title = '(a) Nutrimouse, rCCA CV comp 1 - 2')

network_leaf_data_untar_CV_rcc_tune <- network(leaf_data_untar_CV_rcc_tune, comp = 1:2, interactive = FALSE, lwd.edge = 1, cutoff = 0.30, symkey = FALSE, block.var.names = FALSE, alpha.node = 0.85)
#write.graph(network_leaf_data_untar_CV_rcc_tune$gR, file = "network_leaf_data_untar_CV_rcc_tune.glm")
#write.graph(network_leaf_data_untar_CV_rcc_tune$gR, fil = "network_leaf_data_untar_CV_rcc_tune.glm")
cim(leaf_data_untar_CV_rcc_tune, comp = 1:2, cutoff = 0.3, margins = c(10, 17))

leaf_untar_network <- network_leaf_data_untar_CV_rcc_tune$gR
class(leaf_untar_network)
leaf_untar_network_matrix <- network_leaf_data_untar_CV_rcc_tune$M
dim(leaf_untar_network_matrix)

saveRDS(leaf_untar_network, file = "leaf_untar_network.RDS")

Taxonomy[c(which(row.names(Taxonomy) == "SeqVar143")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar151")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar806")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar836")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar910")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar1115")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar1947")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar2342")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar3093")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar3306")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar3337")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar3447")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar3665")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar3745")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar4568")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar4631")), ] 


leaf_untar_fungal_taxa <- 
  which(row.names(Taxonomy) == "SeqVar143" | row.names(Taxonomy) == "SeqVar151" | row.names(Taxonomy) == "SeqVar806"
        | row.names(Taxonomy) == "SeqVar836" | row.names(Taxonomy) == "SeqVar910" | row.names(Taxonomy) == "SeqVar1115"
        | row.names(Taxonomy) == "SeqVar1947" | row.names(Taxonomy) == "SeqVar2342" | row.names(Taxonomy) == "SeqVar3093"
        | row.names(Taxonomy) == "SeqVar3306" | row.names(Taxonomy) == "SeqVar3337" | row.names(Taxonomy) == "SeqVar3447"
        | row.names(Taxonomy) == "SeqVar3665" | row.names(Taxonomy) == "SeqVar3745" | row.names(Taxonomy) == "SeqVar4568"
        | row.names(Taxonomy) == "SeqVar4631")

write.table (Taxonomy[c(leaf_untar_fungal_taxa) , ], file = "leaf_untar_fungal_taxa.txt")


key_taxa <- c("SeqVar143", "SeqVar151", "SeqVar806", "SeqVar836", "SeqVar910", "SeqVar1115", "SeqVar1947", "SeqVar2342", 
              "SeqVar3093", "SeqVar3306", "SeqVar3337", "SeqVar3447", "SeqVar3665", "SeqVar3745", "SeqVar4568", "SeqVar4631")

Fungal_Guilds <- readRDS(file = "Fungal_Guilds.RDS")

Fungal_Guilds$trophicMode[which(rownames(Fungal_Guilds) == "SeqVar1721")]

leaf_untar_fungal_guilds <- Fungal_Guilds[ which(rownames(Fungal_Guilds) %in% key_taxa) , ]
write.table (leaf_untar_fungal_guilds, file = "leaf_untar_fungal_guilds.txt")





# formed optimised CV rCCA. Need to adjust for this for all four comparisons. 
leaf_data_tar_rcc_tune_opt.l1 <- leaf_data_tar_rcc_tune$opt.lambda1 # extract the optimal lambda values
leaf_data_tar_rcc_tune_opt.l2 <- leaf_data_tar_rcc_tune$opt.lambda2

leaf_data_tar_CV_rcc_tune <- rcc(leaf_data_filtered, leaf_targeted_chem, method = "ridge", lambda1 = leaf_data_tar_rcc_tune_opt.l1, lambda2 = leaf_data_tar_rcc_tune_opt.l2) 
plot(leaf_data_tar_CV_rcc_tune, type = "barplot", main = "Cross Validation") 

plotIndiv(leaf_data_tar_CV_rcc_tune, comp = 1:2, 
          ind.names = leaf_metadata$Habitat,
          group = leaf_metadata$Habitat, rep.space = "XY-variate", 
          legend = TRUE, title = '(a) Nutrimouse, rCCA CV XY-space')

plotArrow(leaf_data_tar_CV_rcc_tune, group = leaf_metadata$Habitat, 
          col.per.group = color.mixo(1:5),
          title = '(a) Nutrimouse, CV method')

plotVar(leaf_data_tar_CV_rcc_tune, var.names = c(TRUE, TRUE),
        cex = c(4, 4), cutoff = 0.35,
        title = '(a) Nutrimouse, rCCA CV comp 1 - 2')

network_leaf_data_tar_CV_rcc_tune <- network(leaf_data_tar_CV_rcc_tune, comp = 1:2, interactive = FALSE, lwd.edge = 1, cutoff = 0.30, symkey = FALSE, block.var.names = FALSE, alpha.node = 0.85)
plot(igraph)
write.graph(network_leaf_data_tar_CV_rcc_tune$gR, "Targeted_fungi_undamaged.graphml", format = "graphml")

saveRDS(network_leaf_data_tar_CV_rcc_tune$gR, file = "leaf_tar_network.RDS")
#write.graph(network_leaf_data_tar_CV_rcc_tune$gR, fil = "network_leaf_data_tar_CV_rcc_tune.glm")
#write.graph(network_leaf_data_tar_CV_rcc_tune$gR, fil = "network_leaf_data_tar_CV_rcc_tune.glm")
cim(leaf_data_tar_CV_rcc_tune, comp = 1:2, cutoff = 0.3, margins = c(10, 17))


Taxonomy[c(which(row.names(Taxonomy) == "SeqVar138")), ] #Too poorly assigned
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar151")), ] #Too poorly defined
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar549")), ] #Leaf spot fungus
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar910")), ] #Yeast - saprotroph
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar2342")), ] #Yeast
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar2979")), ] #Phyllosphere yeast
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar3327")), ] #Too poorly assigned
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar3325")), ] #Cryptococcus (fungi?)
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar3331")), ] #Ergot!!!
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar3337")), ] #Leaf spot
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar3660")), ] #Yeast
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar4568")), ] #Too poorly assigned

leaf_tar_fungal_taxa <- 
  which(row.names(Taxonomy) == "SeqVar138" | row.names(Taxonomy) == "SeqVar151" | row.names(Taxonomy) == "SeqVar549"
        | row.names(Taxonomy) == "SeqVar910" | row.names(Taxonomy) == "SeqVar2342" | row.names(Taxonomy) == "SeqVar2979"
        | row.names(Taxonomy) == "SeqVar3327" | row.names(Taxonomy) == "SeqVar3325" | row.names(Taxonomy) == "SeqVar3331"
        | row.names(Taxonomy) == "SeqVar3337" | row.names(Taxonomy) == "SeqVar3660" | row.names(Taxonomy) == "SeqVar4568")

write.table (Taxonomy[c(leaf_tar_fungal_taxa) , ], file = "leaf_tar_fungal_taxa.txt")


#plot(leaf_data_filtered[(which(colnames(leaf_data_filtered) == "SeqVar1721")), ])

key_taxa <- c("SeqVar138", "SeqVar151", "SeqVar549", "SeqVar910", "SeqVar2342", "SeqVar2979", "SeqVar3327", "SeqVar3325", 
  "SeqVar3331", "SeqVar3337", "SeqVar3660", "SeqVar4568")

Fungal_Guilds <- readRDS(file = "Fungal_Guilds.RDS")

Fungal_Guilds$trophicMode[which(rownames(Fungal_Guilds) == "SeqVar1721")]

leaf_tar_fungal_guilds.txt <- Fungal_Guilds[ which(rownames(Fungal_Guilds) %in% key_taxa) , ]
write.table (leaf_tar_fungal_guilds.txt, file = "leaf_tar_fungal_guilds.txt")




Taxonomy[3740,]




fungi_untar_network <- readRDS(file = "leaf_untar_network.RDS")
#plot(fungi_untar_network, weighted=TRUE, mode="lower")


V(fungi_untar_network)$label
V(fungi_untar_network)$label[[4]] <- "Microbotryales_fOTU836"
V(fungi_untar_network)$label[[5]] <- "Fungi_fOTU910"
V(fungi_untar_network)$label[[14]] <- "Fungi_fOTU3745"
V(fungi_untar_network)$label
#plot(fungi_untar_network, weighted = TRUE, mode="lower")

V(fungi_untar_network)$color <- c(rep("cyan", times = 16), rep("light green", times = 20))
V(fungi_untar_network)$shape <- c(rep("circle", times = 16), rep("circle", times = 20))
V(fungi_untar_network)$shape <- rep("circle", times = 36)
V(fungi_untar_network)$size <- rep(10, times = 36)
#plot(fungi_untar_network, layout = layout.circle, vertex.label.dist=1.5)

V(fungi_untar_network)$color[2] <- c("blue")
V(fungi_untar_network)$shape[2] <- c("square")

V(fungi_untar_network)$color[5] <- c("blue")
V(fungi_untar_network)$shape[5] <- c("square")

V(fungi_untar_network)$color[8] <- c("blue")
V(fungi_untar_network)$shape[8] <- c("square")

V(fungi_untar_network)$color[15] <- c("blue")
V(fungi_untar_network)$shape[15] <- c("square")

V(fungi_untar_network)$color[11] <- c("blue")
V(fungi_untar_network)$shape[11] <- c("square")

#plot(fungi_untar_network, layout = layout.circle, vertex.label.dist=1.5)
E(fungi_untar_network)$color <- ifelse(E(fungi_untar_network)$weight  > 0, "black", "red")
#plot(fungi_untar_network, layout = layout.circle, vertex.label.dist=1.5)

V(fungi_untar_network)$label.cex = 1.2
plot(fungi_untar_network, layout = layout.circle, vertex.label.dist=1.2)




fungi_tar_network <- readRDS(file = "leaf_tar_network.RDS")
#plot(fungi_tar_network, weighted=TRUE, mode="lower")


V(fungi_tar_network)$label
V(fungi_tar_network)$label[[13]] <- "Gardoside"
V(fungi_tar_network)$label[[14]] <- "Melittoside"
V(fungi_tar_network)$label[[15]] <- "10-acetoxy major"
#V(fungi_tar_network)$label
#plot(fungi_tar_network, weighted = TRUE, mode="lower")

V(fungi_tar_network)$color <- c(rep("cyan", times = 12), rep("green3", times = 3))
V(fungi_tar_network)$shape <- c(rep("circle", times = 12), rep("circle", times = 3))
V(fungi_tar_network)$shape <- rep("circle", times = 15)
V(fungi_tar_network)$size <- rep(10, times = 15)
#plot(fungi_tar_network, layout = layout.circle, vertex.label.dist=1.5)

V(fungi_tar_network)$color[2] <- c("blue")
V(fungi_tar_network)$shape[2] <- c("square")

V(fungi_tar_network)$color[4] <- c("blue")
V(fungi_tar_network)$shape[4] <- c("square")

V(fungi_tar_network)$color[5] <- c("blue")
V(fungi_tar_network)$shape[5] <- c("square")

V(fungi_tar_network)$color[10] <- c("blue")
V(fungi_tar_network)$shape[10] <- c("square")

V(fungi_tar_network)$color[12] <- c("blue")
V(fungi_tar_network)$shape[12] <- c("square")


#plot(fungi_tar_network, layout = layout.circle, vertex.label.dist=1.5)
E(fungi_tar_network)$color <- ifelse(E(fungi_tar_network)$weight  > 0, "black", "red")
#plot(fungi_tar_network, layout = layout.circle, vertex.label.dist=1.5)

#V(fungi_tar_network)$label.cex = 1.2
plot(fungi_tar_network, layout = layout.circle, vertex.label.dist=1.5)





















###############################Guild analyses###################################



#devtools::install_github("brendanf/FUNGuildR")
library(FUNGuildR)

string_taxonomy <- str_c(Taxonomy[,1], "; ", Taxonomy[,2], "; ", Taxonomy[,3], "; ", Taxonomy[,4], "; ", Taxonomy[,5], "; ", Taxonomy[,6], "; ", Taxonomy[,7])
string_taxonomy <- data.frame("Taxonomy" = string_taxonomy)
rownames(string_taxonomy) <- rownames(Taxonomy)

#Fungal_Guilds <- funguild_assign(string_taxonomy)
#Fungal_Guilds$guild <- as.factor(Fungal_Guilds$guild)
#saveRDS(Fungal_Guilds, file = "Fungal_Guilds.RDS")
Fungal_Guilds <- readRDS(file = "Fungal_Guilds.RDS")
Fungal_Guilds$guild
dim(Fungal_Guilds)
dim(Taxonomy)

Fungal_Guilds$trophicMode <- as.factor(Fungal_Guilds$trophicMode)
Fungal_Guilds$trophicMode
rownames(Fungal_Guilds) <- rownames(Taxonomy)


dim(rel_data_mean_leaves)
dim(metadata_mean_leaves)
dim(Fungal_Guilds)

Pathogens <- which(Fungal_Guilds$trophicMode == "Pathotroph")
pathogens_guilds <- Fungal_Guilds[Pathogens,]
pathogens_rel_data <- rel_data_mean_leaves[,Pathogens]
pathogens_totals <- rowSums(pathogens_rel_data)
aggregate(pathogens_totals, by = list(
  Category = metadata_mean_leaves$Stratum), FUN = mean)

Saprophytes <- which(Fungal_Guilds$trophicMode == "Saprotroph")
Saprophyte_guilds <- Fungal_Guilds[Saprophytes,]
Saprophyte_rel_data <- rel_data_mean_leaves[,Saprophytes]
Saprophyte_totals <- rowSums(Saprophyte_rel_data)
aggregate(Saprophyte_totals, by = list(
  Category = metadata_mean_leaves$Stratum), FUN = mean)

Symbiotic <- which(Fungal_Guilds$trophicMode == "Symbiotroph")
Symbiotic_guilds <- Fungal_Guilds[Symbiotic,]
Symbiotic_rel_data <- rel_data_mean_leaves[,Symbiotic]
Symbiotic_totals <- rowSums(Symbiotic_rel_data)
aggregate(Symbiotic_totals, by = list(
  Category = metadata_mean_leaves$Stratum), FUN = mean)


x <- (which(Fungal_Guilds$trophicMode == "Symbiotroph" | Fungal_Guilds$trophicMode == "Saprotroph" | Fungal_Guilds$trophicMode == "Pathotroph"))
Fungal_Guilds$trophicMode[-x] <- "Ambigious"

guild_means <- aggregate(t(rel_data_mean_leaves), by = list(Fungal_Guilds$trophicMode), FUN = sum)







#root_undam_paired_means <- aggregate(root_undam_paired_data, by = list(Category = root_undam_paired_metadata$Stratum), mean)
rownames(guild_means) <- guild_means[,1]
guild_means <- guild_means[,-c(1)]
guild_means <- as.data.frame(t(guild_means))
dim(guild_means)


wide_guild_means <- data.frame("Stratum_simplified" = metadata_mean_leaves$Stratum_simplified, "Habitat" = metadata_mean_leaves$Habitat, guild_means)
long_guild_means <-  long.order.data <- melt(wide_guild_means, id.vars = c("Stratum_simplified", "Habitat"), variable.name = "Guild", value.name = "Abundance")
#long_guild_means <- long_guild_means[-c(which(long_guild_means$Guild == "Ambigious")) , ]

ggplot(long_guild_means, aes(x = Guild, y = Abundance, 
  fill = Stratum_simplified)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

Overall_guild_plots <- ggplot(long_guild_means, aes(x = Habitat, y = Abundance, 
  fill = Stratum_simplified)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_grid(. ~ Guild) + scale_fill_manual(values = brewer.pal(4, "Accent")) +
  ggtitle("A")
Overall_guild_plots



t_rel_data <- as.data.frame(t(rel_data_mean_leaves))
t_guild_sums <- aggregate(t_rel_data, by = list(Fungal_Guilds$trophicMode), FUN = sum)
guild_sums <- as.data.frame(t(t_guild_sums))
colnames(guild_sums) <- guild_sums[1,]
guild_sums <- guild_sums[-1,]
#guild_sums <- guild_sums[,-1] #Removes ambigious
guild_sums <- data.frame(lapply(guild_sums, function(x) as.numeric(as.character(x))))

#long_guild_means <- long_guild_means[-c(which(long_guild_means$variable == "Ambigious")), ]
guild_sums <- (make_relative(as.matrix(guild_sums))*100)

wide_guild <- data.frame("Stratum_simplified" = metadata_mean_leaves$Stratum_simplified, "Habitat" = metadata_mean_leaves$Habitat, "PlantID" = metadata_mean_leaves$PlantID, "Site" = metadata_mean_leaves$Site, guild_sums)
long_guild_means <-  long.order.data <- melt(wide_guild, id.vars = c("Stratum_simplified", "Habitat", "PlantID", "Site"), value.name = "Abundance")
#long_guild_means$Guilds <-  rep(Fungal_Guilds$trophicMode, each = nrow(wide_guild))

long_guild_means$Abundance <- as.numeric(long_guild_means$Abundance)
long_guild_means$PlantID <- as.factor(long_guild_means$PlantID)

Figure4 <- ggplot(long_guild_means, aes(x = PlantID, y = Abundance, fill = variable)) + geom_bar(stat="identity", colour="black", size = 0.1) + 
  #facet_grid(Method ~., labeller = as_labeller(new_labels)) + 
  theme(axis.title.x=element_blank(), legend.title = element_blank(), text = element_text(size = 9, angle = 90)) +
  labs(y = "Relative abundance") +
  theme(legend.position = "bottom") +
  theme(strip.text.x = element_text(size = 12))  +
  theme(strip.text.x = element_text(angle = 0))  +
  facet_grid(Stratum_simplified ~ Habitat, space = 'free', scales = 'free') +
  theme(strip.text.y = element_text(size = 12))

colourCount = length(unique(long_guild_means$variable))
getPalette = colorRampPalette(brewer.pal(9, "Accent"))

Figure4 + scale_fill_manual(values = getPalette(colourCount))










#####################Leaf guild analyses #######################

rownames(rel_data_mean_leaves)
rownames(targeted_chem_mean_leaves)
rownames(untargeted_chem_mean_leaves)
metadata_mean_leaves$merging_var
dim(Taxonomy)
dim(rel_data_mean_leaves)
dim(Fungal_Guilds)

Pathogens <- which(Fungal_Guilds$trophicMode == "Pathotroph")
pathogens_Taxonomy <- Taxonomy[ Pathogens, ]

leaf_data <- subset(rel_data_mean_leaves, metadata_mean_leaves$Stratum_simplified == "Leaf")
leaf_metadata <- subset(metadata_mean_leaves, metadata_mean_leaves$Stratum_simplified == "Leaf")
leaf_untargeted_chem <- subset(untargeted_chem_mean_leaves, metadata_mean_leaves$Stratum_simplified == "Leaf")
leaf_targeted_chem <- subset(targeted_chem_mean_leaves, metadata_mean_leaves$Stratum_simplified == "Leaf")

leaf_data <- with(leaf_data, leaf_data[order(leaf_metadata$PlantID),])
leaf_untargeted_chem <- with(leaf_untargeted_chem, leaf_untargeted_chem[order(leaf_metadata$PlantID),])
leaf_targeted_chem <- with(leaf_targeted_chem, leaf_targeted_chem[order(leaf_metadata$PlantID),])
leaf_metadata <- with(leaf_metadata, leaf_metadata[order(leaf_metadata$PlantID),])




#####################Pathogens#######################

leaf_pathogens_data <- leaf_data[ , c(Pathogens)]

x <- which(rowSums(leaf_pathogens_data) == 0)
leaf_pathogens_metadata <- leaf_metadata[ -c(x), ]
leaf_pathogens_untargeted_chem <- leaf_untargeted_chem[ -c(x), ]
leaf_pathogens_targeted_chem <- leaf_targeted_chem[ -c(x), ]
leaf_pathogens_data <- leaf_pathogens_data[ -c(x), ]


leaf_pathogens_bray <- vegdist(leaf_pathogens_data, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
leaf_pathogens_MDS <- monoMDS(leaf_pathogens_bray) # Performs multidimensional scaling for the data
plot(leaf_pathogens_MDS)
rowSums(leaf_pathogens_data)
Nmds1 <- leaf_pathogens_MDS$points[,1] #Extracts dimension 1 as a vector
Nmds2 <- leaf_pathogens_MDS$points[,2] #Extracts dimension 2 as a vector

Nmds <-cbind(Nmds1, Nmds2) #Makes a new sheet with dimension1 and 2
Nmds <- as.data.frame(Nmds) #Creates a new dataframe with dimensions 1 and 2

ggplot(data = Nmds, aes(y = Nmds2, x = Nmds1, Type="p", 
  colour = leaf_pathogens_metadata$Habitat)) + 
  geom_point(size = 4, alpha = 0.5) + labs(x="Dimension 1", y = "Dimension 2", colour = NULL, shape = NULL) + 
  theme(plot.title=element_text(hjust=0)) + 
  ggtitle("D") +
  theme(legend.position = "right") +
  scale_colour_brewer(palette = "Set1") 


adonis2(leaf_pathogens_data ~ leaf_pathogens_metadata$Habitat)

adonis2(leaf_pathogens_data ~ leaf_pathogens_metadata$AnnPrec)
adonis2(leaf_pathogens_data ~ leaf_pathogens_metadata$ColdestQuartMeanTemp)
adonis2(leaf_pathogens_data ~ leaf_pathogens_metadata$MeanAnnTemp)

adonis2(leaf_pathogens_data ~ leaf_pathogens_metadata$pH)
adonis2(leaf_pathogens_data ~ leaf_pathogens_metadata$EC)
adonis2(leaf_pathogens_data ~ leaf_pathogens_metadata$BulkDensity)

#adonis2(leaf_pathogens_data ~ leaf_pathogens_metadata$PCNM1)
#adonis2(leaf_pathogens_data ~ leaf_pathogens_metadata$PCNM2)
adonis2(leaf_pathogens_data ~ leaf_pathogens_metadata$Island)
adonis2(leaf_pathogens_data ~ leaf_pathogens_metadata$Geography3)

adonis2(leaf_pathogens_data ~ leaf_pathogens_metadata$Habitat + leaf_pathogens_metadata$MeanAnnTemp, 
        group = leaf_pathogens_metadata$Site)


pathogens_richness <- rowSums(leaf_pathogens_data > 0)

kruskal.test(pathogens_richness ~ leaf_pathogens_metadata$Site)
kruskal.test(pathogens_richness ~ leaf_pathogens_metadata$Habitat)
kruskal.test(pathogens_richness ~ leaf_pathogens_metadata$Light)

cor.test(pathogens_richness, leaf_pathogens_metadata$AnnPrec, method = "spearman")
cor.test(pathogens_richness, leaf_pathogens_metadata$ColdestQuartMeanTemp, method = "spearman")
cor.test(pathogens_richness, leaf_pathogens_metadata$MeanAnnTemp, method = "spearman")

cor.test(pathogens_richness, leaf_pathogens_metadata$pH, method = "spearman")
cor.test(pathogens_richness, leaf_pathogens_metadata$EC, method = "spearman")
cor.test(pathogens_richness, leaf_pathogens_metadata$BulkDensity, method = "spearman")

cor.test(pathogens_richness, leaf_pathogens_metadata$PCNM1, method = "spearman")
cor.test(pathogens_richness, leaf_pathogens_metadata$PCNM2, method = "spearman")
kruskal.test(pathogens_richness ~ leaf_pathogens_metadata$Island)
kruskal.test(pathogens_richness ~ leaf_pathogens_metadata$Geography3)

leaf_pathogens_bray_data <- vegdist(leaf_pathogens_data, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
leaf_pathogens_untargeted_chem_bray <- vegdist(leaf_pathogens_untargeted_chem, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
leaf_pathogens_targeted_chem_bray <- vegdist(leaf_pathogens_targeted_chem, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence

mantel(leaf_pathogens_bray_data, leaf_pathogens_untargeted_chem_bray)
mantel(leaf_pathogens_bray_data, leaf_pathogens_targeted_chem_bray)

adonis2(leaf_pathogens_bray_data ~ leaf_pathogens_targeted_chem$GeniposideHeight)
adonis2(leaf_pathogens_bray_data ~ leaf_pathogens_targeted_chem$GardosideHeight)
adonis2(leaf_pathogens_bray_data ~ leaf_pathogens_targeted_chem$MelittosideHeight)
adonis2(leaf_pathogens_bray_data ~ leaf_pathogens_targeted_chem$X10.acetoxymajorHeight)
adonis2(leaf_pathogens_bray_data ~ leaf_pathogens_targeted_chem$VerbascosideHeight)
adonis2(leaf_pathogens_bray_data ~ leaf_pathogens_targeted_chem$AucubinHeight)










Saprogens <- which(Fungal_Guilds$trophicMode == "Saprotroph")
saprogens_Taxonomy <- Taxonomy[ Saprogens, ]

leaf_saprogens_data <- leaf_data[ , c(Saprogens)]


x <- which(rowSums(leaf_saprogens_data) == 0)
leaf_saprogens_metadata <- leaf_metadata[ -c(x), ]
leaf_saprogens_untargeted_chem <- leaf_untargeted_chem[ -c(x), ]
leaf_saprogens_targeted_chem <- leaf_targeted_chem[ -c(x), ]
leaf_saprogens_data <- leaf_saprogens_data[ -c(x), ]

rowSums(leaf_saprogens_data)
rowSums(leaf_saprogens_data > 0)

leaf_saprogens_bray <- vegdist(leaf_saprogens_data, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
leaf_saprogens_MDS <- monoMDS(leaf_saprogens_bray) # Performs multidimensional scaling for the data
plot(leaf_saprogens_MDS)
rowSums(leaf_saprogens_data)
Nmds1 <- leaf_saprogens_MDS$points[,1] #Extracts dimension 1 as a vector
Nmds2 <- leaf_saprogens_MDS$points[,2] #Extracts dimension 2 as a vector

Nmds <-cbind(Nmds1, Nmds2) #Makes a new sheet with dimension1 and 2
Nmds <- as.data.frame(Nmds) #Creates a new dataframe with dimensions 1 and 2

ggplot(data = Nmds, aes(y = Nmds2, x = Nmds1, Type="p", 
  colour = leaf_saprogens_metadata$Habitat)) + 
  geom_point(size = 4, alpha = 0.5) + labs(x="Dimension 1", y = "Dimension 2", colour = NULL, shape = NULL) + 
  theme(plot.title=element_text(hjust=0)) + 
  ggtitle("D") +
  theme(legend.position = "right") +
  scale_colour_brewer(palette = "Set1") 


adonis2(leaf_saprogens_data ~ leaf_saprogens_metadata$Habitat)

adonis2(leaf_saprogens_data ~ leaf_saprogens_metadata$AnnPrec)
adonis2(leaf_saprogens_data ~ leaf_saprogens_metadata$ColdestQuartMeanTemp)
adonis2(leaf_saprogens_data ~ leaf_saprogens_metadata$MeanAnnTemp)

adonis2(leaf_saprogens_data ~ leaf_saprogens_metadata$pH)
adonis2(leaf_saprogens_data ~ leaf_saprogens_metadata$EC)
adonis2(leaf_saprogens_data ~ leaf_saprogens_metadata$BulkDensity)

#adonis2(leaf_saprogens_data ~ leaf_saprogens_metadata$PCNM1)
#adonis2(leaf_saprogens_data ~ leaf_saprogens_metadata$PCNM2)
adonis2(leaf_saprogens_data ~ leaf_saprogens_metadata$Island)
adonis2(leaf_saprogens_data ~ leaf_saprogens_metadata$Geography3)

adonis2(leaf_saprogens_data ~ leaf_saprogens_metadata$Habitat + leaf_saprogens_metadata$Island + 
        leaf_saprogens_metadata$ColdestQuartMeanTemp, 
        group = leaf_saprogens_metadata$Site)


saprogens_richness <- rowSums(leaf_saprogens_data > 0)

kruskal.test(saprogens_richness ~ leaf_saprogens_metadata$Site)
kruskal.test(saprogens_richness ~ leaf_saprogens_metadata$Habitat)
kruskal.test(saprogens_richness ~ leaf_saprogens_metadata$Light)

cor.test(saprogens_richness, leaf_saprogens_metadata$AnnPrec, method = "spearman")
cor.test(saprogens_richness, leaf_saprogens_metadata$ColdestQuartMeanTemp, method = "spearman")
cor.test(saprogens_richness, leaf_saprogens_metadata$MeanAnnTemp, method = "spearman")

cor.test(saprogens_richness, leaf_saprogens_metadata$pH, method = "spearman")
cor.test(saprogens_richness, leaf_saprogens_metadata$EC, method = "spearman")
cor.test(saprogens_richness, leaf_saprogens_metadata$BulkDensity, method = "spearman")

cor.test(saprogens_richness, leaf_saprogens_metadata$PCNM1, method = "spearman")
cor.test(saprogens_richness, leaf_saprogens_metadata$PCNM2, method = "spearman")
kruskal.test(saprogens_richness ~ leaf_saprogens_metadata$Island)
kruskal.test(saprogens_richness ~ leaf_saprogens_metadata$Geography3)

leaf_saprogens_bray_data <- vegdist(leaf_saprogens_data, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
leaf_saprogens_untargeted_chem_bray <- vegdist(leaf_saprogens_untargeted_chem, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
leaf_saprogens_targeted_chem_bray <- vegdist(leaf_saprogens_targeted_chem, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence

mantel(leaf_saprogens_bray_data, leaf_saprogens_untargeted_chem_bray)
mantel(leaf_saprogens_bray_data, leaf_saprogens_targeted_chem_bray)

adonis2(leaf_saprogens_bray_data ~ leaf_saprogens_targeted_chem$GeniposideHeight)
adonis2(leaf_saprogens_bray_data ~ leaf_saprogens_targeted_chem$GardosideHeight)
adonis2(leaf_saprogens_bray_data ~ leaf_saprogens_targeted_chem$MelittosideHeight)
adonis2(leaf_saprogens_bray_data ~ leaf_saprogens_targeted_chem$X10.acetoxymajorHeight)
adonis2(leaf_saprogens_bray_data ~ leaf_saprogens_targeted_chem$VerbascosideHeight)
adonis2(leaf_saprogens_bray_data ~ leaf_saprogens_targeted_chem$AucubinHeight)














Symbionts <- which(Fungal_Guilds$trophicMode == "Symbiotroph")
symbiotrophs_Taxonomy <- Taxonomy[ Symbionts, ]

leaf_symbiotrophs_data <- leaf_data[ , c(Symbionts)]


x <- which(rowSums(leaf_symbiotrophs_data) == 0)
leaf_symbiotrophs_metadata <- leaf_metadata[ -c(x), ]
leaf_symbiotrophs_untargeted_chem <- leaf_untargeted_chem[ -c(x), ]
leaf_symbiotrophs_targeted_chem <- leaf_targeted_chem[ -c(x), ]
leaf_symbiotrophs_data <- leaf_symbiotrophs_data[ -c(x), ]

rowSums(leaf_symbiotrophs_data)
rowSums(leaf_symbiotrophs_data > 0)

leaf_symbiotrophs_bray <- vegdist(leaf_symbiotrophs_data, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
leaf_symbiotrophs_MDS <- monoMDS(leaf_symbiotrophs_bray) # Performs multidimensional scaling for the data
plot(leaf_symbiotrophs_MDS)
rowSums(leaf_symbiotrophs_data)
Nmds1 <- leaf_symbiotrophs_MDS$points[,1] #Extracts dimension 1 as a vector
Nmds2 <- leaf_symbiotrophs_MDS$points[,2] #Extracts dimension 2 as a vector

Nmds <-cbind(Nmds1, Nmds2) #Makes a new sheet with dimension1 and 2
Nmds <- as.data.frame(Nmds) #Creates a new dataframe with dimensions 1 and 2

ggplot(data = Nmds, aes(y = Nmds2, x = Nmds1, Type="p", 
                        colour = leaf_symbiotrophs_metadata$Habitat)) + 
  geom_point(size = 4, alpha = 0.5) + labs(x="Dimension 1", y = "Dimension 2", colour = NULL, shape = NULL) + 
  theme(plot.title=element_text(hjust=0)) + 
  ggtitle("D") +
  theme(legend.position = "right") +
  scale_colour_brewer(palette = "Set1") 


adonis2(leaf_symbiotrophs_data ~ leaf_symbiotrophs_metadata$Habitat)

adonis2(leaf_symbiotrophs_data ~ leaf_symbiotrophs_metadata$AnnPrec)
adonis2(leaf_symbiotrophs_data ~ leaf_symbiotrophs_metadata$ColdestQuartMeanTemp)
adonis2(leaf_symbiotrophs_data ~ leaf_symbiotrophs_metadata$MeanAnnTemp)

adonis2(leaf_symbiotrophs_data ~ leaf_symbiotrophs_metadata$pH)
adonis2(leaf_symbiotrophs_data ~ leaf_symbiotrophs_metadata$EC)
adonis2(leaf_symbiotrophs_data ~ leaf_symbiotrophs_metadata$BulkDensity)

#adonis2(leaf_symbiotrophs_data ~ leaf_symbiotrophs_metadata$PCNM1)
#adonis2(leaf_symbiotrophs_data ~ leaf_symbiotrophs_metadata$PCNM2)
adonis2(leaf_symbiotrophs_data ~ leaf_symbiotrophs_metadata$Island)
adonis2(leaf_symbiotrophs_data ~ leaf_symbiotrophs_metadata$Geography3)

adonis2(leaf_symbiotrophs_data ~ leaf_symbiotrophs_metadata$Habitat + leaf_symbiotrophs_metadata$ColdestQuartMeanTemp, 
        group = leaf_symbiotrophs_metadata$Site)


symbiotrophs_richness <- rowSums(leaf_symbiotrophs_data > 0)

kruskal.test(symbiotrophs_richness ~ leaf_symbiotrophs_metadata$Site)
kruskal.test(symbiotrophs_richness ~ leaf_symbiotrophs_metadata$Habitat)
kruskal.test(symbiotrophs_richness ~ leaf_symbiotrophs_metadata$Light)

cor.test(symbiotrophs_richness, leaf_symbiotrophs_metadata$AnnPrec, method = "spearman")
cor.test(symbiotrophs_richness, leaf_symbiotrophs_metadata$ColdestQuartMeanTemp, method = "spearman")
cor.test(symbiotrophs_richness, leaf_symbiotrophs_metadata$MeanAnnTemp, method = "spearman")

cor.test(symbiotrophs_richness, leaf_symbiotrophs_metadata$pH, method = "spearman")
cor.test(symbiotrophs_richness, leaf_symbiotrophs_metadata$EC, method = "spearman")
cor.test(symbiotrophs_richness, leaf_symbiotrophs_metadata$BulkDensity, method = "spearman")

cor.test(symbiotrophs_richness, leaf_symbiotrophs_metadata$PCNM1, method = "spearman")
cor.test(symbiotrophs_richness, leaf_symbiotrophs_metadata$PCNM2, method = "spearman")
kruskal.test(symbiotrophs_richness ~ leaf_symbiotrophs_metadata$Island)
kruskal.test(symbiotrophs_richness ~ leaf_symbiotrophs_metadata$Geography3)

leaf_symbiotrophs_bray_data <- vegdist(leaf_symbiotrophs_data, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
leaf_symbiotrophs_untargeted_chem_bray <- vegdist(leaf_symbiotrophs_untargeted_chem, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
leaf_symbiotrophs_targeted_chem_bray <- vegdist(leaf_symbiotrophs_targeted_chem, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence

mantel(leaf_symbiotrophs_bray_data, leaf_symbiotrophs_untargeted_chem_bray)
mantel(leaf_symbiotrophs_bray_data, leaf_symbiotrophs_targeted_chem_bray)

adonis2(leaf_symbiotrophs_bray_data ~ leaf_symbiotrophs_targeted_chem$GeniposideHeight)
adonis2(leaf_symbiotrophs_bray_data ~ leaf_symbiotrophs_targeted_chem$GardosideHeight)
adonis2(leaf_symbiotrophs_bray_data ~ leaf_symbiotrophs_targeted_chem$MelittosideHeight)
adonis2(leaf_symbiotrophs_bray_data ~ leaf_symbiotrophs_targeted_chem$X10.acetoxymajorHeight)
adonis2(leaf_symbiotrophs_bray_data ~ leaf_symbiotrophs_targeted_chem$VerbascosideHeight)
adonis2(leaf_symbiotrophs_bray_data ~ leaf_symbiotrophs_targeted_chem$AucubinHeight)









