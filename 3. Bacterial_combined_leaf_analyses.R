

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

setwd("/Users/christopherbarnes/Dropbox/Danish_Plantago_microbiomes/Bacteria/DADA2_files")


rel_data <- readRDS(file = "data_data_no_mito.RDS")
metadata <- readRDS(file = "metadata_data_no_mito.RDS")
untargeted_chem <- readRDS(file = "final_untargeted_chem_no_mito.RDS")
targeted_chem <- readRDS(file = "final_targeted_chem_no_mito.RDS")
Taxonomy <- readRDS(file = "final_taxonomy_no_mito.RDS")

metadata$Stratum_simplified <- str_replace(metadata$Stratum, "Damaged leaf", "Leaf")
metadata$Stratum_simplified <- str_replace(metadata$Stratum_simplified, "Undamaged leaf", "Leaf")
metadata$merging_var <- paste(metadata$PlantID, metadata$Stratum_simplified, sep = "_")

rel_data_mean_leaves <- aggregate(rel_data, by = list(Category = metadata$merging_var), FUN = mean)
rownames(rel_data_mean_leaves) <- rel_data_mean_leaves[ , 1]
rel_data_mean_leaves <- rel_data_mean_leaves[ , -c(1)]
colnames(rel_data_mean_leaves) <- colnames(rel_data)


rel_data_mean_leaves <- readRDS(file = "rel_data_mean_leaves_no_mito.RDS")
metadata_mean_leaves <- metadata[!duplicated(metadata$merging_var),]
dim(metadata_mean_leaves)

metadata_mean_leaves <- metadata_mean_leaves[order(metadata_mean_leaves$merging_var) , ]
dim(metadata_mean_leaves)
rel_data_mean_leaves <- rel_data_mean_leaves[order(rownames(rel_data_mean_leaves)),]
dim(rel_data_mean_leaves)

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

#saveRDS(jaccards, "Jaccards_no_mito.RDS")
jaccards <- readRDS(file = "Jaccards_no_mito.RDS")



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

saveRDS(Big_data, file = "Bacterial_big_data_no_mito.RDS")















###########################Network analyses#################################



leaf_data <- subset(rel_data_mean_leaves, metadata_mean_leaves$Stratum_simplified == "Leaf")
leaf_metadata <- subset(metadata_mean_leaves, metadata_mean_leaves$Stratum_simplified == "Leaf")
leaf_untargeted_chem <- subset(untargeted_chem_mean_leaves, metadata_mean_leaves$Stratum_simplified == "Leaf")
leaf_targeted_chem <- subset(targeted_chem_mean_leaves, metadata_mean_leaves$Stratum_simplified == "Leaf")

leaf_data <- with(leaf_data, leaf_data[order(leaf_metadata$PlantID),])
leaf_untargeted_chem <- with(leaf_untargeted_chem, leaf_untargeted_chem[order(leaf_metadata$PlantID),])
leaf_targeted_chem <- with(leaf_targeted_chem, leaf_targeted_chem[order(leaf_metadata$PlantID),])
leaf_metadata <- with(leaf_metadata, leaf_metadata[order(leaf_metadata$PlantID),])

colnames(rel_data_mean_leaves) <- gsub("SeqVar", "ASV", colnames(rel_data_mean_leaves))
genus_descriptors <- paste(Taxonomy$Genus, colnames(rel_data_mean_leaves), sep = "_")
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
#leaf_data_untar_rcc_tune <- saveRDS(leaf_data_untar_rcc_tune, file = "leaf_data_untar_rcc_tune_no_mito")
#leaf_data_untar_rcc_tune <- rcc(leaf_data_filtered, leaf_untargeted_chem_filtered, method = 'shrinkage') 
leaf_data_untar_rcc_tune <- readRDS(file = "leaf_data_untar_rcc_tune_no_mito")

#leaf_data_tar_rcc_tune <- tune.rcc(leaf_data_filtered, leaf_targeted_chem, grid1 = grid1, grid2 = grid2, validation = "loo") 
#leaf_data_tar_rcc_tune <- saveRDS(leaf_data_tar_rcc_tune, file = "leaf_data_tar_rcc_tune_no_mito")
leaf_data_tar_rcc_tune <- readRDS(file = "leaf_data_tar_rcc_tune_no_mito")



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
        cex = c(4, 4), cutoff = 0.3,
        title = '(a) Nutrimouse, rCCA CV comp 1 - 2')

network_leaf_data_untar_CV_rcc_tune <- network(leaf_data_untar_CV_rcc_tune, comp = 1:2, interactive = FALSE, lwd.edge = 1, cutoff = 0.40, symkey = FALSE, block.var.names = FALSE, alpha.node = 0.85)
#write.graph(network_leaf_data_untar_CV_rcc_tune$gR, file = "network_leaf_data_untar_CV_rcc_tune.glm")
#write.graph(network_leaf_data_untar_CV_rcc_tune$gR, fil = "network_leaf_data_untar_CV_rcc_tune.glm")
cim(leaf_data_untar_CV_rcc_tune, comp = 1:2, cutoff = 0.3, margins = c(10, 17))

leaf_untar_network <- network_leaf_data_untar_CV_rcc_tune$gR
class(leaf_untar_network)
leaf_untar_network_matrix <- network_leaf_data_untar_CV_rcc_tune$M
dim(leaf_untar_network_matrix)
saveRDS(leaf_untar_network, file = "leaf_untar_network.RDS")
plot(leaf_untar_network)

Taxonomy[c(which(row.names(Taxonomy) == "SeqVar2")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar10")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar15")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar24")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar25")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar27")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar30")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar31")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar67")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar69")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar112")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar141")), ]
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar164")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar165")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar331")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar670")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar728")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar729")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar1010")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar1586")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar1765")), ] 


leaf_untar_fungal_taxa <- 
  which(row.names(Taxonomy) == "SeqVar2" | row.names(Taxonomy) == "SeqVar10" | row.names(Taxonomy) == "SeqVar15"
        | row.names(Taxonomy) == "SeqVar24" | row.names(Taxonomy) == "SeqVar25" | row.names(Taxonomy) == "SeqVar27"
        | row.names(Taxonomy) == "SeqVar30" | row.names(Taxonomy) == "SeqVar31" | row.names(Taxonomy) == "SeqVar67"
        | row.names(Taxonomy) == "SeqVar69" | row.names(Taxonomy) == "SeqVar112" | row.names(Taxonomy) == "SeqVar141"
        | row.names(Taxonomy) == "SeqVar164" | row.names(Taxonomy) == "SeqVar165" | row.names(Taxonomy) == "SeqVar331"
        | row.names(Taxonomy) == "SeqVar670" | row.names(Taxonomy) == "SeqVar728" | row.names(Taxonomy) == "SeqVar729"
        | row.names(Taxonomy) == "SeqVar1010" | row.names(Taxonomy) == "SeqVar1586" | row.names(Taxonomy) == "SeqVar1765")

write.table (Taxonomy[c(leaf_untar_fungal_taxa) , ], file = "leaf_untar_fungal_taxa.txt")




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
#write.graph(network_leaf_data_tar_CV_rcc_tune$gR, fil = "network_leaf_data_tar_CV_rcc_tune.glm")
#write.graph(network_leaf_data_tar_CV_rcc_tune$gR, fil = "network_leaf_data_tar_CV_rcc_tune.glm")
cim(leaf_data_tar_CV_rcc_tune, comp = 1:2, cutoff = 0.3, margins = c(10, 17))

saveRDS(network_leaf_data_tar_CV_rcc_tune$gR, file = "leaf_tar_network.RDS")


Taxonomy[c(which(row.names(Taxonomy) == "SeqVar4")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar6")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar10")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar12")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar25")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar29")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar30")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar35")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar63")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar89")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar92")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar112")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar153")), ]
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar196")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar211")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar233")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar226")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar291")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar334")), ]
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar336")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar437")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar502")), ]
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar527")), ] 
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar785")), ]
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar1674")), ]
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar1918")), ]
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar2001")), ]
Taxonomy[c(which(row.names(Taxonomy) == "SeqVar2602")), ] 

leaf_tar_fungal_taxa <- 
  which(row.names(Taxonomy) == "SeqVar4" | row.names(Taxonomy) == "SeqVar6" | row.names(Taxonomy) == "SeqVar10" |
          row.names(Taxonomy) =="SeqVar12" | row.names(Taxonomy) == "SeqVar25" | row.names(Taxonomy) == "SeqVar29" |
          row.names(Taxonomy) == "SeqVar30" | row.names(Taxonomy) == "SeqVar35" | row.names(Taxonomy) == "SeqVar63" |
          row.names(Taxonomy) == "SeqVar89" | row.names(Taxonomy) == "SeqVar92" | row.names(Taxonomy) == "SeqVar112" |
          row.names(Taxonomy) == "SeqVar153" | row.names(Taxonomy) == "SeqVar196" | row.names(Taxonomy) == "SeqVar211" |
          row.names(Taxonomy) == "SeqVar233" | row.names(Taxonomy) == "SeqVar226" | row.names(Taxonomy) == "SeqVar291" |
          row.names(Taxonomy) == "SeqVar334" | row.names(Taxonomy) == "SeqVar336" | row.names(Taxonomy) == "SeqVar437" |
          row.names(Taxonomy) == "SeqVar502" | row.names(Taxonomy) == "SeqVar527" | row.names(Taxonomy) == "SeqVar785" |
          row.names(Taxonomy) == "SeqVar1674" | row.names(Taxonomy) == "SeqVar1918" | row.names(Taxonomy) == "SeqVar2001" |
          row.names(Taxonomy) == "SeqVar2602")

write.table (Taxonomy[c(leaf_tar_fungal_taxa) , ], file = "leaf_tar_fungal_taxa.txt")


#plot(leaf_data_filtered[(which(colnames(leaf_data_filtered) == "SeqVar1721")), ])






par(mfrow=c(1,2))

bac_untar_network <- readRDS(file = "leaf_untar_network.RDS")
#plot(bac_untar_network, weighted=TRUE, mode="lower")


V(bac_untar_network)$label
#V(bac_untar_network)$label[[1]] <- "Sphingomonas_SeqVar2"
#V(bac_untar_network)$label[[2]] <- "Sphingomonas_SeqVar10" 
#V(bac_untar_network)$label[[3]] <- "Pseudomonas_SeqVar24"
#V(bac_untar_network)$label[[4]] <- "Microbacteriaceae SeqVar31"
#V(bac_untar_network)$label[[5]] <- "Hymenobacter_SeqVar69"
#V(bac_untar_network)$label[[6]] <- "Aureimonas_SeqVar164"
#V(bac_untar_network)$label[[7]] <- "Paenibacillus_SeqVar1586"
#V(bac_untar_network)$label[[8]] <- "Hymenobacter_SeqVar1765"
#plot(bac_untar_network, weighted = TRUE, mode="lower")

V(bac_untar_network)$color <- c(rep("cyan", times = 8), rep("light green", times = 12))
V(bac_untar_network)$shape <- c(rep("circle", times = 8), rep("circle", times = 12))

V(bac_untar_network)$shape[2] <- c("square")
V(bac_untar_network)$color[2] <- c("blue")
#V(bac_untar_network)$shape <- rep("circle", times = 20)
V(bac_untar_network)$size <- rep(10, times = 20)
#plot(bac_untar_network, layout = layout.circle, vertex.label.dist=1.5)

#plot(bac_untar_network, layout = layout.circle, vertex.label.dist=1.5)
E(bac_untar_network)$color <- ifelse(E(bac_untar_network)$weight  > 0, "black", "red")
#plot(bac_untar_network, layout = layout.circle, vertex.label.dist=1.5)

V(bac_untar_network)$label.cex = 1.2
plot(bac_untar_network, layout = layout.circle, vertex.label.dist=1.5)

bac_untar_network_taxa <- Taxonomy[ which(rownames(Taxonomy) %in% V(bac_untar_network)$label[1:30]) , ]
write.table(bac_untar_network_taxa, file = "bac_untar_network_taxa.txt")



bac_tar_network <- readRDS(file = "leaf_tar_network.RDS")
#plot(bac_tar_network, weighted=TRUE, mode="lower")

V(bac_tar_network)$label
V(bac_tar_network)$label[[29]] <- "Geniposide"
V(bac_tar_network)$label[[30]] <- "Gardoside"
V(bac_tar_network)$label[[31]] <- "Verbascoside"

V(bac_tar_network)$label
#plot(bac_tar_network, weighted = TRUE, mode="lower")

V(bac_tar_network)$color <- c(rep("cyan", times = 28), rep("green3", times = 3))
V(bac_tar_network)$shape <- c(rep("circle", times = 28), rep("circle", times = 3))

V(bac_tar_network)$color[3] <- c("blue")
V(bac_tar_network)$color[5] <- c("blue")
V(bac_tar_network)$color[7] <- c("blue")
V(bac_tar_network)$color[12] <- c("blue")

V(bac_tar_network)$shape[3] <- c("square")
V(bac_tar_network)$shape[5] <- c("square")
V(bac_tar_network)$shape[7] <- c("square")
V(bac_tar_network)$shape[12] <- c("square")

#V(bac_tar_network)$shape <- rep("circle", times = 31)
V(bac_tar_network)$size <- rep(10, times = 31)
#plot(bac_tar_network, layout = layout.circle, vertex.label.dist=1.5)

#plot(bac_tar_network, layout = layout.circle, vertex.label.dist=1.5)
E(bac_tar_network)$color <- ifelse(E(bac_tar_network)$weight  > 0, "black", "red")
#plot(bac_tar_network, layout = layout.circle, vertex.label.dist=1.5)

V(bac_tar_network)$label.cex = 1.2
plot(bac_tar_network, layout = layout.circle, vertex.label.dist=1.5)

bac_tar_network_taxa <- Taxonomy[ which(rownames(Taxonomy) %in% V(bac_tar_network)$label[1:30]) , ]
write.table(bac_tar_network_taxa, file = "bac_tar_network_taxa.txt")









