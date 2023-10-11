
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

#setwd("/Users/jsv736/Dropbox/Danish_Plantago_microbiomes/Bacteria/DADA2_files")
#setwd("c:/Users/drcjb/Dropbox/Danish_Plantago_microbiomes/Bacteria/DADA2_files")
setwd("/Users/christopherbarnes/Dropbox/Danish_Plantago_microbiomes/Bacteria/DADA2_files")

data <- read.table(file = "dada2_seq_table_with_species.txt", sep = "\t", header = TRUE)

Taxonomy <- data[,216]
data <- data[,-c(216)]
rownames(data) <- data[,1]
data <- data[, -c(1)]

rownames(data)
colnames(data)

data <- as.data.frame(t(data))
full_taxonomy <- as.data.frame(str_split_fixed(Taxonomy, ";", 7))

colnames(full_taxonomy) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rownames(full_taxonomy) <- colnames(data[,1])
full_taxonomy[,7] <- str_c(full_taxonomy[,6], "_", full_taxonomy[,7])


Reads <- rowSums(data)
summary(rowSums(data))
a <- which(Reads < 1000)
#metadata <- metadata[-c(a),]
#dim(metadata)
data <- data[-c(a),]
#data <- nthroot(as.matrix(data), 4)
#data <- (make_relative(data)*100)
data <- as.data.frame(data)
rowSums(data)

a <- rownames(data)
#write.table(a, file = "sample_names.txt")
#a <- rownames(data)
#write.table(a, file = "sample_names2.txt")
data <- (data[-c(161, 107, 188, 194, 137, 141, 149, 158, 178, 197, 11, 18, 176, 46, 182, 185,
                 54, 186, 34, 46, 191, 55, 193, 131), ])
dim(data)
#a <- rownames(data)
#write.table(a, file = "sample_names3.txt")

metadata <- read.table(file = "Metadata_new.txt", sep = "\t", header = TRUE)
dim(metadata)

NegData <- data[c(which(metadata$Site == "Neg")),]
dim(NegData)
NegMetadata <- metadata[c(which(metadata$Site == "Neg")),]
dim(NegMetadata)
Contaminants <- colSums(NegData > 0)
data <- data[, -c(which(Contaminants > 1))]
dim(data)
full_taxonomy <- full_taxonomy[-c(which(Contaminants > 1)), ]
dim(full_taxonomy)

Negatives <- which(metadata$Site == "Neg")
data <- data[-c(Negatives),]
dim(data)
metadata <- metadata[-c(Negatives),]
dim(metadata)

Bacteria_Archeaea_taxa <- which(full_taxonomy$Kingdom == "Bacteria" | full_taxonomy$Kingdom =="Archaea")
data <- data[, c(Bacteria_Archeaea_taxa)]
dim(data)
full_taxonomy <- full_taxonomy[c(Bacteria_Archeaea_taxa), ]
dim(full_taxonomy)


Chloroplast_taxa <- which(full_taxonomy$Order == "Chloroplast")
#rowSums(data[,Chloroplast_taxa])
#rowSums(data[,Chloroplast_taxa])/rowSums(data)*100
data <- data[, -c(Chloroplast_taxa)]
dim(data)
full_taxonomy <- full_taxonomy[-c(Chloroplast_taxa), ]
dim(full_taxonomy)

Mitochondrial_taxa  <- which(full_taxonomy$Family == "Mitochondria")
data <- data[, -c(Mitochondrial_taxa)]
dim(data)
full_taxonomy <- full_taxonomy[-c(Mitochondrial_taxa), ]
dim(full_taxonomy)


sort(rowSums(data))
metadata[which(rownames(data) == "S252"),]
metadata[which(rownames(data) == "S248"),]
metadata[which(rownames(data) == "S033"),]
metadata[which(rownames(data) == "S255"),]
metadata[which(rownames(data) == "S253"),]


x <- (colSums(data != 0)) > 0 #OTU as columns
data <- data[, x == "TRUE"] #Removes taxa that are all 0s
dim(data)
full_taxonomy <- full_taxonomy[ x == "TRUE",] #Removes taxa that are all 0s
dim(full_taxonomy)
sort(rowSums(data))

x <- which(rowSums(data) < 200) #Originally 1,000, but lose a lot more samples after chloroplast removal. Therefore 200.
sort(rowSums(data[x,])) 
metadata[c(x),]

data <- data[-c(x),]
dim(data)
metadata <- metadata[-c(x),]
dim(metadata)

min(rowSums(data))
Original_Richness <- rarefy(data, 219) #Originally 1,003, but now 553 after chloroplast removal.
metadata <- data.frame(metadata, "Richness" = Original_Richness)

data <- nthroot(as.matrix(data), 4)
#data <- as.data.frame(make_relative(as.matrix(data))*100)
rel_data <- as.data.frame(make_relative(as.matrix(data))*100)




metadata$pH <- as.numeric(as.character(metadata$pH))
metadata$MeanAnnTemp <- as.numeric(as.character(metadata$MeanAnnTemp))
metadata$ColdestQuartMeanTemp <- as.numeric(as.character(metadata$ColdestQuartMeanTemp))
metadata <- droplevels(metadata)

metadata$Lat <- as.numeric(as.character(metadata$Lat))
metadata$Long <- as.numeric(as.character(metadata$Long))

metadata$Stratum <- factor (metadata$Status, 
                            levels = c("U", "D", "R", "S"), 
                            labels = c("Undamaged leaf", "Damaged leaf", "Root", "Soil"))

metadata$Site <- as.factor(metadata$Site)

spat_data <- as.matrix(cbind(metadata$Long, metadata$Lat))
dist_matrix <-  distm(spat_data) #Creates a distance matrix for how far apart samples are
PCNM <- pcnm(dist_matrix)
PCNM1 <- (PCNM$vectors[,1]*10) #For scaling (keeps them varying at similar levels to other variables)
PCNM2 <- (PCNM$vectors[,2]*10) #For scaling (keeps them varying at similar levels to other variables)
metadata$PCNM1 <- PCNM1
metadata$PCNM2 <- PCNM2





norm_metadata_prep <- as.data.frame((cbind("pH" = metadata$pH, "EC" = metadata$EC, "BulkDensity" = metadata$BulkDensity_g.cm3,
  "MeanAnnTemp" = metadata$MeanAnnTemp, 
  "ColdestQuartMeanTemp" = metadata$ColdestQuartMeanTemp, "AnnPrec" = metadata$AnnPrec, 
  "Richness" = Original_Richness)))
norm_metadata_prep <- as.data.frame(sapply( norm_metadata_prep, as.numeric ))

detach("package:pracma", unload=TRUE)

library(sjmisc)
#library(caret)
#norm_metadata_prep1 <- nearZeroVar(norm_metadata_prep) #Checking to see if there are low variance variables that need to be removed. 
norm_metadata <- as.data.frame(std(norm_metadata_prep, include.fac = TRUE, append = FALSE))
norm_metadata <- center(norm_metadata, append = FALSE)
norm_metadata <- as.data.frame(scale(norm_metadata, scale = TRUE))
colnames(norm_metadata) <- colnames(norm_metadata_prep)

metadata <- data.frame("Site" = metadata$Site, "PlantID" = metadata$PlantID, 
                       "Stratum" = metadata$Stratum, "Habitat" = metadata$Habitat, 
                       "Light" = metadata$Light, "Geography3" = metadata$Geography3,
                       "Island" = metadata$Island, "PCNM1" = metadata$PCNM1, 
                       "PCNM2" = metadata$PCNM2,
                       norm_metadata)

metadata$Habitat[which(metadata$Habitat == "SeasideMeadow")] <- "Meadow"
metadata <- droplevels(metadata)

metadata$Original_Richness <- Original_Richness

metadata$Stratum_simplified <- str_replace(metadata$Stratum, "Damaged leaf", "Leaf")
metadata$Stratum_simplified <- str_replace(metadata$Stratum_simplified, "Undamaged leaf", "Leaf")



library(pracma)
untargeted_chem <- read.table(file = "Untargeted_chem.txt", sep = "\t", header = TRUE)
untargeted_chem <- nthroot(as.matrix(untargeted_chem), 4)
untargeted_chem <- (make_relative(untargeted_chem)*100)
untargeted_chem <- as.data.frame(untargeted_chem)
untargeted_chem <- untargeted_chem[-c(x),]
dim(untargeted_chem)
#untargeted_chem <- as.data.frame(std(untargeted_chem, include.fac = TRUE, append = FALSE))

targeted_chem <- read.table(file = "Targeted_chem.txt", sep = "\t", header = TRUE)
targeted_chem <- nthroot(as.matrix(targeted_chem), 4)
targeted_chem <- (make_relative(targeted_chem)*100)
targeted_chem <- as.data.frame(targeted_chem)
targeted_chem <- targeted_chem[-c(x),]
dim(targeted_chem)
#targeted_chem <- as.data.frame(std(targeted_chem, include.fac = TRUE, append = FALSE))


detach("package:pracma", unload=TRUE)
rowSums(rel_data)
all_bray_data <- vegdist(rel_data, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
all_data_MDS <- monoMDS(all_bray_data) # Performs multidimensional scaling for the data
plot(all_data_MDS)

p <- prcomp(rel_data, scale = FALSE)
s <- summary(p)
unclass(p)
screeplot(p)
biplot(p)

Nmds1 <- all_data_MDS$points[,1] #Extracts dimension 1 as a vector
Nmds2 <- all_data_MDS$points[,2] #Extracts dimension 2 as a vector

####Dimensions are the new variables, which are 'trends' in your data based on shared variation
####Dimension 1 explains the most shared variation of bacteria, dimension 2 the secondmost, then 3rd

Nmds <-cbind(Nmds1, Nmds2) #Makes a new sheet with dimension1 and 2
Nmds <- as.data.frame(Nmds) #Creates a new dataframe with dimensions 1 and 2

#saveRDS(Nmds, file = "/Users/christopherbarnes/Dropbox/Danish_Plantago_microbiomes/Combined_analyses/bacteria_nMDS_no_mito.RDS")
#saveRDS(metadata, file = "/Users/christopherbarnes/Dropbox/Danish_Plantago_microbiomes/Combined_analyses/bacteria_metadata_no_mito.RDS")


Fig1D <- ggplot(data = Nmds, aes(y = Nmds2, x = Nmds1, Type="p", 
                                 colour = metadata$Stratum,
                                 label = rownames(metadata))) + 
  geom_point(size = 4, alpha = 0.5) + labs(x="Dimension 1", y = "Dimension 2", colour = NULL, shape = NULL) + 
  theme(plot.title=element_text(hjust=0)) + 
  ggtitle("D") +
  theme(legend.position = "right") +
  #theme(legend.text=element_text(size=8)) +
  #scale_shape_manual(name = "Skin status", values = c(20, 18, 15)) + #theme(legend.position = "none") +
  scale_colour_brewer(palette = "Set1") #+ facet_grid(metadata$Plant ~ metadata$Site) 
Fig1D

rowSums(rel_data)
rel_data_totals <- aggregate(rel_data, by = list(
  Category = metadata$Stratum_simplified), FUN = sum)

rownames(rel_data_totals) <- rel_data_totals[,1]
rel_data_totals <- rel_data_totals[,-c(1)]
rel_data_totals <- as.data.frame(t(rel_data_totals))
rel_data_totals <- as.data.frame(((rel_data_totals > 0) > 0)*1L)
venn(rel_data_totals)

library("ggVennDiagram")
library(ggvenn)

Soil_OTU <- rownames(rel_data_totals)[which(rel_data_totals$Soil == 1)]
Root_OTU <- rownames(rel_data_totals)[which(rel_data_totals$Root == 1)]
#Damaged_leaf_OTU <- rownames(rel_data_totals)[which(rel_data_totals$`Damaged leaf` == 1)]
#Undamaged_leaf_OTU <- rownames(rel_data_totals)[which(rel_data_totals$`Undamaged leaf` == 1)]
Leaf_OTU <- rownames(rel_data_totals)[which(rel_data_totals$Leaf == 1)]
#x <- list("Soil" = Soil_OTU, "Root" = Root_OTU, "Damaged \n leaf" = Damaged_leaf_OTU, "Undamaged \n leaf" = Undamaged_leaf_OTU)
x <- list("Soil" = Soil_OTU, "Root" = Root_OTU, "Leaf" = Leaf_OTU)

ggVennDiagram(x)
ggvenn(x)

#saveRDS(x, file = "Bacterial_Venn_data.RDS")
#saveRDS(x, file = "Bacterial_Venn_data_simplified.RDS")
#saveRDS(x, file = "Bacterial_Venn_data_simplified_no_mito.RDS")


#Compositional analyses
adonis2(rel_data ~ metadata$Stratum)
adonis2(rel_data ~ metadata$PlantID)

ord <- all_bray_data
groups <- metadata$Stratum
mod <- betadisper(ord, groups)
permutest(mod)
#if the dispersion is different between groups, then examine
plot(mod)
boxplot(mod)
mod.HSD <- TukeyHSD(mod)
mod.HSD
plot(mod.HSD)


kruskal.test(metadata$Richness ~ metadata$Stratum)
dunnTest(metadata$Richness ~ metadata$Stratum, data = metadata)

Fig1E <- ggplot(metadata, aes(x = Stratum, y = Original_Richness, 
  fill = Stratum)) + geom_violin(alpha = 0.5, draw_quantiles = c(0.5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_brewer(palette = "Set1") + ggtitle("E") +
  labs(fill = NULL, x = NULL) +
  theme(legend.position = "none") +
  coord_flip() +
  labs(y = "Richness")
Fig1E


ggplot(metadata, aes(x = Habitat, y = BulkDensity)) + geom_violin(alpha = 0.5, draw_quantiles = c(0.5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  scale_fill_brewer(palette = "Set1") + ggtitle("B") +
  labs(fill = NULL, x = NULL) +
  theme(legend.position = "none")



aggregate(metadata$Original_Richness, by = list(Category = metadata$Stratum), FUN = mean)
aggregate(metadata$Original_Richness, by = list(Category = metadata$Stratum), FUN = sd)

#####################Root data analysis#######################
root_data <- subset(rel_data, metadata$Stratum == "Root")
root_metadata <- subset(metadata, metadata$Stratum == "Root")

root_data <- with(root_data, root_data[order(root_metadata$PlantID),])
root_metadata <- with(root_metadata, root_metadata[order(root_metadata$PlantID),])

#adonis2(root_data ~ root_metadata$Site)
adonis2(root_data ~ root_metadata$Habitat)
#adonis2(root_data ~ root_metadata$Light)

adonis2(root_data ~ root_metadata$AnnPrec)
adonis2(root_data ~ root_metadata$ColdestQuartMeanTemp)
adonis2(root_data ~ root_metadata$MeanAnnTemp)

adonis2(root_data ~ root_metadata$pH)
adonis2(root_data ~ root_metadata$EC)
adonis2(root_data ~ root_metadata$BulkDensity)

#adonis2(root_data ~ root_metadata$PCNM1)
#adonis2(root_data ~ root_metadata$PCNM2)
adonis2(root_data ~ root_metadata$Island)
adonis2(root_data ~ root_metadata$Geography3)

adonis2(root_data ~ root_metadata$Habitat + root_metadata$Island +
  root_metadata$pH + root_metadata$AnnPrec + root_metadata$ColdestQuartMeanTemp + 
  root_metadata$MeanAnnTemp, group = root_metadata$Site)


kruskal.test(root_metadata$Richness ~ root_metadata$Site)
kruskal.test(root_metadata$Richness ~ root_metadata$Habitat)
#kruskal.test(root_metadata$Richness ~ root_metadata$Light)

cor.test(root_metadata$Richness, root_metadata$AnnPrec, method = "spearman")
cor.test(root_metadata$Richness, root_metadata$ColdestQuartMeanTemp, method = "spearman")
cor.test(root_metadata$Richness, root_metadata$MeanAnnTemp, method = "spearman")

cor.test(root_metadata$Richness, root_metadata$pH, method = "spearman")
cor.test(root_metadata$Richness, root_metadata$EC, method = "spearman")
cor.test(root_metadata$Richness, root_metadata$BulkDensity, method = "spearman")

cor.test(root_metadata$Richness, root_metadata$PCNM1, method = "spearman")
cor.test(root_metadata$Richness, root_metadata$PCNM2, method = "spearman")
kruskal.test(root_metadata$Richness ~ root_metadata$Island)
kruskal.test(root_metadata$Richness ~ root_metadata$Geography3)

root_bray_vector <- as.vector(vegdist(root_data, Type = "bray", binary = T))
label = rep("Root", times = length(root_bray_vector))
root_sim_data <- data.frame("Similarity" = root_bray_vector, label)


res2 <- rcorr(as.matrix(cbind(root_metadata[,-c(1:7)])))
corrplot(res2$r, type="upper", 
         p.mat = res2$P, sig.level = 0.01, insig = "blank")

#saveRDS(root_data, file = "root_metadata.RDS")
#saveRDS(root_metadata, file = "root_metadata.RDS")

#####################Soil data analysis#######################
soil_data <- subset(rel_data, metadata$Stratum == "Soil")
soil_metadata <- subset(metadata, metadata$Stratum == "Soil")

soil_data <- with(soil_data, soil_data[order(soil_metadata$PlantID),])
soil_metadata <- with(soil_metadata, soil_metadata[order(soil_metadata$PlantID),])

#adonis2(soil_data ~ soil_metadata$Site)
adonis2(soil_data ~ soil_metadata$Habitat)
#adonis2(soil_data ~ soil_metadata$Light)

adonis2(soil_data ~ soil_metadata$AnnPrec)
adonis2(soil_data ~ soil_metadata$ColdestQuartMeanTemp)
adonis2(soil_data ~ soil_metadata$MeanAnnTemp)

adonis2(soil_data ~ soil_metadata$pH)
adonis2(soil_data ~ soil_metadata$EC)
adonis2(soil_data ~ soil_metadata$BulkDensity)

#adonis2(soil_data ~ soil_metadata$PCNM1)
#adonis2(soil_data ~ soil_metadata$PCNM2)
adonis2(soil_data ~ soil_metadata$Island)
adonis2(soil_data ~ soil_metadata$Geography3)

adonis2(soil_data ~ soil_metadata$Habitat + soil_metadata$Island +
          soil_metadata$pH + soil_metadata$AnnPrec + soil_metadata$MeanAnnTemp, group = soil_metadata$Site)

kruskal.test(soil_metadata$Richness ~ soil_metadata$Site)
kruskal.test(soil_metadata$Richness ~ soil_metadata$Habitat)
kruskal.test(soil_metadata$Richness ~ soil_metadata$Light)

cor.test(soil_metadata$Richness, soil_metadata$AnnPrec, method = "spearman")
cor.test(soil_metadata$Richness, soil_metadata$ColdestQuartMeanTemp, method = "spearman")
cor.test(soil_metadata$Richness, soil_metadata$MeanAnnTemp, method = "spearman")

cor.test(soil_metadata$Richness, soil_metadata$pH, method = "spearman")
cor.test(soil_metadata$Richness, soil_metadata$EC, method = "spearman")
cor.test(soil_metadata$Richness, soil_metadata$BulkDensity, method = "spearman")

cor.test(soil_metadata$Richness, soil_metadata$PCNM1, method = "spearman")
cor.test(soil_metadata$Richness, soil_metadata$PCNM2, method = "spearman")
kruskal.test(soil_metadata$Richness ~ soil_metadata$Island)
kruskal.test(soil_metadata$Richness ~ soil_metadata$Geography3)

soil_bray_vector <- as.vector(vegdist(soil_data, Type = "bray", binary = T))
label = rep("Soil", times = length(soil_bray_vector))
soil_sim_data <- data.frame("Similarity" = soil_bray_vector, label)

#saveRDS(soil_data, file = "root_metadata.RDS")
#saveRDS(soil_metadata, file = "root_metadata.RDS")



#####################Dam data analysis#######################

dam_data <- subset(rel_data, metadata$Stratum == "Damaged leaf")
dam_metadata <- subset(metadata, metadata$Stratum == "Damaged leaf")
dam_untargeted_chem <- subset(untargeted_chem, metadata$Stratum == "Damaged leaf")
dam_targeted_chem <- subset(targeted_chem, metadata$Stratum == "Damaged leaf")

dam_data <- with(dam_data, dam_data[order(dam_metadata$PlantID),])
dam_untargeted_chem <- with(dam_untargeted_chem, dam_untargeted_chem[order(dam_metadata$PlantID),])
dam_targeted_chem <- with(dam_targeted_chem, dam_targeted_chem[order(dam_metadata$PlantID),])
dam_metadata <- with(dam_metadata, dam_metadata[order(dam_metadata$PlantID),])

#adonis2(dam_data ~ dam_metadata$Site)
adonis2(dam_data ~ dam_metadata$Habitat)
#adonis2(dam_data ~ dam_metadata$Light)

adonis2(dam_data ~ dam_metadata$AnnPrec)
adonis2(dam_data ~ dam_metadata$ColdestQuartMeanTemp)
adonis2(dam_data ~ dam_metadata$MeanAnnTemp)

adonis2(dam_data ~ dam_metadata$pH)
adonis2(dam_data ~ dam_metadata$EC)
adonis2(dam_data ~ dam_metadata$BulkDensity)

#adonis2(dam_data ~ dam_metadata$PCNM1)
#adonis2(dam_data ~ dam_metadata$PCNM2)
adonis2(dam_data ~ dam_metadata$Island)
adonis2(dam_data ~ dam_metadata$Geography3)

adonis2(dam_data ~ dam_metadata$Habitat + dam_metadata$Island, group = dam_metadata$Site)


kruskal.test(dam_metadata$Richness ~ dam_metadata$Site)
kruskal.test(dam_metadata$Richness ~ dam_metadata$Habitat)
kruskal.test(dam_metadata$Richness ~ dam_metadata$Light)

cor.test(dam_metadata$Richness, dam_metadata$AnnPrec, method = "spearman")
cor.test(dam_metadata$Richness, dam_metadata$ColdestQuartMeanTemp, method = "spearman")
cor.test(dam_metadata$Richness, dam_metadata$MeanAnnTemp, method = "spearman")

cor.test(dam_metadata$Richness, dam_metadata$pH, method = "spearman")
cor.test(dam_metadata$Richness, dam_metadata$EC, method = "spearman")
cor.test(dam_metadata$Richness, dam_metadata$BulkDensity, method = "spearman")

cor.test(dam_metadata$Richness, dam_metadata$PCNM1, method = "spearman")
cor.test(dam_metadata$Richness, dam_metadata$PCNM2, method = "spearman")
kruskal.test(dam_metadata$Richness ~ dam_metadata$Island)
kruskal.test(dam_metadata$Richness ~ dam_metadata$Geography3)

dam_bray_data <- vegdist(dam_data, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
dam_untargeted_bray_data <- vegdist(dam_untargeted_chem, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
#dam_untargeted_bray_data <- vegdist((dam_untargeted_chem - (min(dam_untargeted_chem))), Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
dam_targeted_bray_data <- vegdist(dam_targeted_chem, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
#dam_targeted_bray_data <- vegdist((dam_targeted_chem - (min(dam_targeted_chem))), Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence

mantel(dam_bray_data, dam_untargeted_bray_data)
mantel(dam_bray_data, dam_targeted_bray_data)


dam_bray_vector <- as.vector(vegdist(dam_data, Type = "bray", binary = T))
label = rep("Damaged leaf", times = length(dam_bray_vector))
dam_sim_data <- data.frame("Similarity" = dam_bray_vector, label)




adonis2(dam_bray_data ~ dam_targeted_chem$GeniposideHeight)
adonis2(dam_bray_data ~ dam_targeted_chem$GardosideHeight)
adonis2(dam_bray_data ~ dam_targeted_chem$MelittosideHeight)
adonis2(dam_bray_data ~ dam_targeted_chem$X10.acetoxymajorHeight)
adonis2(dam_bray_data ~ dam_targeted_chem$VerbascosideHeight) #significant
adonis2(dam_bray_data ~ dam_targeted_chem$AucubinHeight)



dam_MDS <- monoMDS(dam_targeted_bray_data)
Nmds1 <- dam_MDS$points[,1] #Extracts dimension 1 as a vector
Nmds2 <- dam_MDS$points[,2] #Extracts dimension 2 as a vector

####Dimensions are the new variables, which are 'trends' in your data based on shared variation
####Dimension 1 explains the most shared variation of bacteria, dimension 2 the secondmost, then 3rd

dam_Nmds <-cbind(Nmds1, Nmds2) #Makes a new sheet with dimension1 and 2
dam_Nmds <- as.data.frame(dam_Nmds) #Creates a new dataframe with dimensions 1 and 2

#saveRDS(Nmds, file = "/Users/christopherbarnes/Dropbox/Danish_Plantago_microbiomes/Combined_analyses/bacteria_nMDS.RDS")
#saveRDS(metadata, file = "/Users/christopherbarnes/Dropbox/Danish_Plantago_microbiomes/Combined_analyses/bacteria_metadata.RDS")


ggplot(data = dam_Nmds, aes(y = Nmds2, x = Nmds1, Type="p", colour = dam_targeted_chem$VerbascosideHeight)) + 
  geom_point(size = 4, alpha = 0.5) + labs(x="Dimension 1", y = "Dimension 2", colour = NULL, shape = NULL) + 
  theme(plot.title=element_text(hjust=0)) + 
  ggtitle("D") +
  theme(legend.position = "right") 



#saveRDS(dam_data, file = "dam_data.RDS")
#saveRDS(dam_metadata, file = "dam_metadata.RDS")
#saveRDS(dam_untargeted_chem, file = "dam_untargeted_chem.RDS")
#saveRDS(dam_targeted_chem, file = "dam_targeted_chem.RDS")


#####################Undam data analysis#######################

undam_data <- subset(rel_data, metadata$Stratum == "Undamaged leaf")
undam_metadata <- subset(metadata, metadata$Stratum == "Undamaged leaf")
undam_untargeted_chem <- subset(untargeted_chem, metadata$Stratum == "Undamaged leaf")
undam_targeted_chem <- subset(targeted_chem, metadata$Stratum == "Undamaged leaf")

undam_data <- with(undam_data, undam_data[order(undam_metadata$PlantID),])
undam_untargeted_chem <- with(undam_untargeted_chem, undam_untargeted_chem[order(undam_metadata$PlantID),])
undam_targeted_chem <- with(undam_targeted_chem, undam_targeted_chem[order(undam_metadata$PlantID),])
undam_metadata <- with(undam_metadata, undam_metadata[order(undam_metadata$PlantID),])


#adonis2(undam_data ~ undam_metadata$Site)
adonis2(undam_data ~ undam_metadata$Habitat)
#adonis2(undam_data ~ undam_metadata$Light)

adonis2(undam_data ~ undam_metadata$AnnPrec) #1
adonis2(undam_data ~ undam_metadata$ColdestQuartMeanTemp)
adonis2(undam_data ~ undam_metadata$MeanAnnTemp)

adonis2(undam_data ~ undam_metadata$pH)
adonis2(undam_data ~ undam_metadata$EC)
adonis2(undam_data ~ undam_metadata$BulkDensity)

#adonis2(undam_data ~ undam_metadata$PCNM1)
#adonis2(undam_data ~ undam_metadata$PCNM2)
adonis2(undam_data ~ undam_metadata$Island)
adonis2(undam_data ~ undam_metadata$Geography3)

adonis2(undam_data ~ undam_metadata$Habitat + undam_metadata$Island, 
        group = undam_metadata$Site)


kruskal.test(undam_metadata$Richness ~ undam_metadata$Site)
kruskal.test(undam_metadata$Richness ~ undam_metadata$Habitat)
kruskal.test(undam_metadata$Richness ~ undam_metadata$Light)

cor.test(undam_metadata$Richness, undam_metadata$AnnPrec, method = "spearman")
cor.test(undam_metadata$Richness, undam_metadata$ColdestQuartMeanTemp, method = "spearman")
cor.test(undam_metadata$Richness, undam_metadata$MeanAnnTemp, method = "spearman")

cor.test(undam_metadata$Richness, undam_metadata$pH, method = "spearman")
cor.test(undam_metadata$Richness, undam_metadata$EC, method = "spearman")
cor.test(undam_metadata$Richness, undam_metadata$BulkDensity, method = "spearman")

cor.test(undam_metadata$Richness, undam_metadata$PCNM1, method = "spearman")
cor.test(undam_metadata$Richness, undam_metadata$PCNM2, method = "spearman")
kruskal.test(undam_metadata$Richness ~ undam_metadata$Island)
kruskal.test(undam_metadata$Richness ~ undam_metadata$Geography3)

undam_bray_data <- vegdist(undam_data, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
undam_untargeted_bray_data <- vegdist(undam_untargeted_chem, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
#undam_untargeted_bray_data <- vegdist((undam_untargeted_chem - (min(undam_untargeted_chem))), Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
undam_targeted_bray_data <- vegdist(undam_targeted_chem, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
#undam_untargeted_bray_data <- vegdist((undam_untargeted_chem - (min(undam_untargeted_chem))), Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence

mantel(undam_bray_data, undam_untargeted_bray_data)
mantel(undam_bray_data, undam_targeted_bray_data)
mantel(undam_untargeted_bray_data, undam_targeted_bray_data)

adonis2(undam_bray_data ~ undam_targeted_chem$GeniposideHeight) #Significant
adonis2(undam_bray_data ~ undam_targeted_chem$GardosideHeight) #Significant
adonis2(undam_bray_data ~ undam_targeted_chem$MelittosideHeight) 
adonis2(undam_bray_data ~ undam_targeted_chem$X10.acetoxymajorHeight) #Significant
adonis2(undam_bray_data ~ undam_targeted_chem$VerbascosideHeight)
adonis2(undam_bray_data ~ undam_targeted_chem$AucubinHeight)


undam_MDS <- monoMDS(undam_targeted_bray_data)
Nmds1 <- undam_MDS$points[,1] #Extracts dimension 1 as a vector
Nmds2 <- undam_MDS$points[,2] #Extracts dimension 2 as a vector

####Dimensions are the new variables, which are 'trends' in your data based on shared variation
####Dimension 1 explains the most shared variation of bacteria, dimension 2 the secondmost, then 3rd

undam_Nmds <-cbind(Nmds1, Nmds2) #Makes a new sheet with dimension1 and 2
undam_Nmds <- as.data.frame(undam_Nmds) #Creates a new dataframe with dimensions 1 and 2

#saveRDS(Nmds, file = "/Users/christopherbarnes/Dropbox/Danish_Plantago_microbiomes/Combined_analyses/bacteria_nMDS.RDS")
#saveRDS(metadata, file = "/Users/christopherbarnes/Dropbox/Danish_Plantago_microbiomes/Combined_analyses/bacteria_metadata.RDS")


ggplot(data = undam_Nmds, aes(y = Nmds2, x = Nmds1, Type="p", colour = undam_targeted_chem$X10.acetoxymajorHeight)) + 
  geom_point(size = 4, alpha = 0.5) + labs(x="Dimension 1", y = "Dimension 2", colour = NULL, shape = NULL) + 
  theme(plot.title=element_text(hjust=0)) + 
  ggtitle("D") +
  theme(legend.position = "right") 


library(Hmisc)
res <- rcorr(as.matrix(dam_targeted_chem))
round(res$P, 3)


undam_bray_vector <- as.vector(vegdist(undam_data, Type = "bray", binary = T))
label = rep("Undamaged leaf", times = length(undam_bray_vector))
undam_sim_data <- data.frame("Similarity" = undam_bray_vector, label)



#saveRDS(undam_data, file = "undam_data.RDS")
#saveRDS(undam_metadata, file = "undam_metadata.RDS")
#saveRDS(undam_untargeted_chem, file = "undam_untargeted_chem.RDS")
#saveRDS(undam_targeted_chem, file = "undam_targeted_chem.RDS")






long_sim_data <- rbind(soil_sim_data, root_sim_data, dam_sim_data, undam_sim_data)



library(plyr)
mu <- ddply(long_sim_data, "label", summarise, grp.mean = median(Similarity))

long_sim_data$label <- factor (long_sim_data$label, 
                               levels = c("Undamaged leaf", "Damaged leaf", "Root", "Soil"), 
                               labels = c("Undamaged leaf", "Damaged leaf", "Root", "Soil"))

Fig1F <- ggplot(long_sim_data, aes(x = Similarity, color = label)) + geom_density(size = 1) + 
  geom_vline(data = mu, aes(xintercept = grp.mean, color = label),
             linetype="dashed", size = 1) + scale_colour_brewer(palette = "Set1") +
  labs(x="S?rensen dissimilarity", y = "Density", color = "") + ggtitle("F") +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none")
Fig1F


kruskal.test(long_sim_data$Similarity ~ long_sim_data$label)
dunnTest(long_sim_data$Similarity ~ long_sim_data$label, data = long_sim_data)

Fig1.part2 <- Fig1D + Fig1E + Fig1F + plot_layout(guides = "collect")

ggsave(file="Fig1.pdf", width = 297, height = 210, units = "mm")







##################Leaf comparisons


leaves <- which(metadata$Stratum == "Undamaged leaf"| metadata$Stratum == "Damaged leaf")
leaf_data <- rel_data[c(leaves),]
leaf_metadata <- metadata[c(leaves),]


adonis2(leaf_data ~ leaf_metadata$Stratum, group = leaf_metadata$PlantID)
t.test(leaf_metadata$Richness ~ leaf_metadata$Stratum)







#######################################################################
########################Cross strata comparisons#######################
#######################################################################

#################Soil comparisons######################################

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
a <- which(table(c(soil_metadata$PlantID, undam_metadata$PlantID)) > 1)
a <- names(a)

soil_data_undam_soil <- subset(soil_data, subset = soil_metadata$PlantID %in% c(a))
dim(soil_data_undam_soil)
soil_metadata_undam_soil <- subset(soil_metadata, subset = soil_metadata$PlantID %in% c(a))
dim(soil_metadata_undam_soil)

undam_data_undam_soil <- subset(undam_data, subset = undam_metadata$PlantID %in% c(a))
dim(undam_data_undam_soil)
undam_metadata_undam_soil <- subset(undam_metadata, subset = undam_metadata$PlantID %in% c(a))
dim(undam_metadata_undam_soil)

undam_targeted_chem_undam_soil <- subset(undam_targeted_chem, subset = undam_metadata$PlantID %in% c(a))
dim(undam_targeted_chem_undam_soil)
undam_untargeted_chem_undam_soil <- subset(undam_untargeted_chem, subset = undam_metadata$PlantID %in% c(a))
dim(undam_untargeted_chem_undam_soil)

undam_data_undam_soil_bray <- vegdist(undam_data_undam_soil, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
soil_data_undam_soil_bray <- vegdist(soil_data_undam_soil, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
undam_targeted_chem_undam_soil_bary <- vegdist(undam_targeted_chem_undam_soil, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
undam_untargeted_chem_undam_soil_bary <- vegdist(undam_untargeted_chem_undam_soil, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence

mantel(undam_data_undam_soil_bray, soil_data_undam_soil_bray)
cor.test(soil_metadata_undam_soil$Richness, undam_metadata_undam_soil$Richness, method = "spearman")

mantel(soil_data_undam_soil_bray, undam_targeted_chem_undam_soil_bary)
mantel(soil_data_undam_soil_bray, undam_untargeted_chem_undam_soil_bary)

adonis2(soil_data_undam_soil_bray ~ undam_targeted_chem_undam_soil$GeniposideHeight)
adonis2(soil_data_undam_soil_bray ~ undam_targeted_chem_undam_soil$GardosideHeight)
adonis2(soil_data_undam_soil_bray ~ undam_targeted_chem_undam_soil$MelittosideHeight)
adonis2(soil_data_undam_soil_bray ~ undam_targeted_chem_undam_soil$X10.acetoxymajorHeight)
adonis2(soil_data_undam_soil_bray ~ undam_targeted_chem_undam_soil$VerbascosideHeight)
adonis2(soil_data_undam_soil_bray ~ undam_targeted_chem_undam_soil$AucubinHeight)


undam_data_undam_soil_binary <- (undam_data_undam_soil > 0) # logical, or
undam_data_undam_soil_binary <- (undam_data_undam_soil_binary > 0)*1L # integer 01
undam_data_undam_soil_binary <- as.data.frame(undam_data_undam_soil_binary)

soil_data_undam_soil_binary <- (soil_data_undam_soil > 0) # logical, or
soil_data_undam_soil_binary <- (soil_data_undam_soil_binary > 0)*1L # integer 01
soil_data_undam_soil_binary <- as.data.frame(soil_data_undam_soil_binary)

both_undam_soil_binary <- list()
just_soil_data_undam_soil_binary <- list()
just_undam_data_undam_soil_binary <- list()

for (i in 1:nrow(soil_data_undam_soil_binary)){
  both_undam_soil_binary[[i]] <- length(which(undam_data_undam_soil_binary[i,] == 1 & soil_data_undam_soil_binary[i,] == 1))
  just_undam_data_undam_soil_binary[[i]] <- length(which(undam_data_undam_soil_binary[i,] == 1 & soil_data_undam_soil_binary[i,] == 0))
  just_soil_data_undam_soil_binary[[i]] <- length(which(undam_data_undam_soil_binary[i,] == 0 & soil_data_undam_soil_binary[i,] == 1))
}

both_undam_soil_binary <- unlist(both_undam_soil_binary, use.names = FALSE)
just_undam_data_undam_soil_binary <- unlist(just_undam_data_undam_soil_binary, use.names = FALSE)
just_soil_data_undam_soil_binary <- unlist(just_soil_data_undam_soil_binary, use.names = FALSE)
undam_soil_binary_comparisons <- data.frame(both_undam_soil_binary, just_undam_data_undam_soil_binary, just_soil_data_undam_soil_binary)
undam_soil_binary_comparisons

undam_soil_binary_comparisons_means <- colMeans(undam_soil_binary_comparisons)





#soil-damage comparisons
a <- which(table(c(soil_metadata$PlantID, dam_metadata$PlantID)) > 1)
a <- names(a)

soil_data_dam_soil <- subset(soil_data, subset = soil_metadata$PlantID %in% c(a))
dim(soil_data_dam_soil)
soil_metadata_dam_soil <- subset(soil_metadata, subset = soil_metadata$PlantID %in% c(a))
dim(soil_metadata_dam_soil)

dam_data_dam_soil <- subset(dam_data, subset = dam_metadata$PlantID %in% c(a))
dim(dam_data_dam_soil)
dam_metadata_dam_soil <- subset(dam_metadata, subset = dam_metadata$PlantID %in% c(a))
dim(dam_metadata_dam_soil)

dam_targeted_chem_dam_soil <- subset(dam_targeted_chem, subset = dam_metadata$PlantID %in% c(a))
dim(dam_targeted_chem_dam_soil)
dam_untargeted_chem_dam_soil <- subset(dam_untargeted_chem, subset = dam_metadata$PlantID %in% c(a))
dim(dam_untargeted_chem_dam_soil)

dam_data_dam_soil_bray <- vegdist(dam_data_dam_soil, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
soil_data_dam_soil_bray <- vegdist(soil_data_dam_soil, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
dam_targeted_chem_dam_soil_bary <- vegdist(dam_targeted_chem_dam_soil, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
dam_untargeted_chem_dam_soil_bary <- vegdist(dam_untargeted_chem_dam_soil, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence

mantel(dam_data_dam_soil_bray, soil_data_dam_soil_bray)
cor.test(dam_metadata_dam_soil$Richness, soil_metadata_dam_soil$Richness, method = "spearman")

mantel(soil_data_dam_soil_bray, dam_targeted_chem_dam_soil_bary)
mantel(soil_data_dam_soil_bray, dam_untargeted_chem_dam_soil_bary)

adonis2(soil_data_dam_soil_bray ~ dam_targeted_chem_dam_soil$GeniposideHeight)
adonis2(soil_data_dam_soil_bray ~ dam_targeted_chem_dam_soil$GardosideHeight)
adonis2(soil_data_dam_soil_bray ~ dam_targeted_chem_dam_soil$MelittosideHeight)
adonis2(soil_data_dam_soil_bray ~ dam_targeted_chem_dam_soil$X10.acetoxymajorHeight)
adonis2(soil_data_dam_soil_bray ~ dam_targeted_chem_dam_soil$VerbascosideHeight)
adonis2(soil_data_dam_soil_bray ~ dam_targeted_chem_dam_soil$AucubinHeight)



dam_data_dam_soil_binary <- (dam_data_dam_soil > 0) # logical, or
dam_data_dam_soil_binary <- (dam_data_dam_soil_binary > 0)*1L # integer 01
dam_data_dam_soil_binary <- as.data.frame(dam_data_dam_soil_binary)

soil_data_dam_soil_binary <- (soil_data_dam_soil > 0) # logical, or
soil_data_dam_soil_binary <- (soil_data_dam_soil_binary > 0)*1L # integer 01
soil_data_dam_soil_binary <- as.data.frame(soil_data_dam_soil_binary)

both_dam_soil_binary <- list()
just_soil_data_dam_soil_binary <- list()
just_dam_data_dam_soil_binary <- list()

for (i in 1:nrow(soil_data_dam_soil_binary)){
  both_dam_soil_binary[[i]] <- length(which(dam_data_dam_soil_binary[i,] == 1 & soil_data_dam_soil_binary[i,] == 1))
  just_dam_data_dam_soil_binary[[i]] <- length(which(dam_data_dam_soil_binary[i,] == 1 & soil_data_dam_soil_binary[i,] == 0))
  just_soil_data_dam_soil_binary[[i]] <- length(which(dam_data_dam_soil_binary[i,] == 0 & soil_data_dam_soil_binary[i,] == 1))
}

both_dam_soil_binary <- unlist(both_dam_soil_binary, use.names = FALSE)
just_dam_data_dam_soil_binary <- unlist(just_dam_data_dam_soil_binary, use.names = FALSE)
just_soil_data_dam_soil_binary <- unlist(just_soil_data_dam_soil_binary, use.names = FALSE)
dam_soil_binary_comparisons <- data.frame(both_dam_soil_binary, just_dam_data_dam_soil_binary, just_soil_data_dam_soil_binary)
dam_soil_binary_comparisons

dam_soil_binary_comparisons_means <- colMeans(dam_soil_binary_comparisons)




#########Root comparisons############

#Root-damage comparisons
a <- which(table(c(root_metadata$PlantID, dam_metadata$PlantID)) > 1)
a <- names(a)

root_data_dam_root <- subset(root_data, subset = root_metadata$PlantID %in% c(a))
dim(root_data_dam_root)
root_metadata_dam_root <- subset(root_metadata, subset = root_metadata$PlantID %in% c(a))
dim(root_metadata_dam_root)

dam_data_dam_root <- subset(dam_data, subset = dam_metadata$PlantID %in% c(a))
dim(dam_data_dam_root)
dam_metadata_dam_root <- subset(dam_metadata, subset = dam_metadata$PlantID %in% c(a))
dim(dam_metadata_dam_root)

dam_targeted_chem_dam_root <- subset(dam_targeted_chem, subset = dam_metadata$PlantID %in% c(a))
dim(dam_targeted_chem_dam_root)
dam_untargeted_chem_dam_root <- subset(dam_untargeted_chem, subset = dam_metadata$PlantID %in% c(a))
dim(dam_untargeted_chem_dam_root)


dam_data_dam_root_bray <- vegdist(dam_data_dam_root, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
root_data_dam_root_bray <- vegdist(root_data_dam_root, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
dam_targeted_chem_dam_root_bary <- vegdist(dam_targeted_chem_dam_root, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
dam_untargeted_chem_dam_root_bary <- vegdist(dam_untargeted_chem_dam_root, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence

mantel(dam_data_dam_root_bray, root_data_dam_root_bray)
cor.test(dam_metadata_dam_root$Richness, root_metadata_dam_root$Richness, method = "spearman")

mantel(root_data_dam_root_bray, dam_targeted_chem_dam_root_bary)
mantel(root_data_dam_root_bray, dam_untargeted_chem_dam_root_bary)

adonis2(root_data_dam_root_bray ~ dam_targeted_chem_dam_root$GeniposideHeight)
adonis2(root_data_dam_root_bray ~ dam_targeted_chem_dam_root$GardosideHeight)
adonis2(root_data_dam_root_bray ~ dam_targeted_chem_dam_root$MelittosideHeight)
adonis2(root_data_dam_root_bray ~ dam_targeted_chem_dam_root$X10.acetoxymajorHeight)
adonis2(root_data_dam_root_bray ~ dam_targeted_chem_dam_root$VerbascosideHeight)
adonis2(root_data_dam_root_bray ~ dam_targeted_chem_dam_root$AucubinHeight)


root_data_dam_root_binary <- (root_data_dam_root > 0) # logical, or
root_data_dam_root_binary <- (root_data_dam_root_binary > 0)*1L # integer 01
root_data_dam_root_binary <- as.data.frame(root_data_dam_root_binary)

dam_data_dam_root_binary <- (dam_data_dam_root > 0) # logical, or
dam_data_dam_root_binary <- (dam_data_dam_root_binary > 0)*1L # integer 01
dam_data_dam_root_binary <- as.data.frame(dam_data_dam_root_binary)

both_dam_root_binary <- list()
just_dam_data_dam_root_binary <- list()
just_root_data_dam_root_binary <- list()

for (i in 1:nrow(dam_data_dam_root_binary)){
  both_dam_root_binary[[i]] <- length(which(root_data_dam_root_binary[i,] == 1 & dam_data_dam_root_binary[i,] == 1))
  just_root_data_dam_root_binary[[i]] <- length(which(root_data_dam_root_binary[i,] == 1 & dam_data_dam_root_binary[i,] == 0))
  just_dam_data_dam_root_binary[[i]] <- length(which(root_data_dam_root_binary[i,] == 0 & dam_data_dam_root_binary[i,] == 1))
}

both_dam_root_binary <- unlist(both_dam_root_binary, use.names = FALSE)
just_root_data_dam_root_binary <- unlist(just_root_data_dam_root_binary, use.names = FALSE)
just_dam_data_dam_root_binary <- unlist(just_dam_data_dam_root_binary, use.names = FALSE)
dam_root_binary_comparisons <- data.frame(both_dam_root_binary, just_root_data_dam_root_binary, just_dam_data_dam_root_binary)
dam_root_binary_comparisons

dam_root_binary_comparisons_means <- colMeans(dam_root_binary_comparisons)






#Root-undamaged comparisons
a <- which(table(c(root_metadata$PlantID, undam_metadata$PlantID)) > 1)
a <- names(a)

root_data_undam_root <- subset(root_data, subset = root_metadata$PlantID %in% c(a))
dim(root_data_undam_root)
root_metadata_undam_root <- subset(root_metadata, subset = root_metadata$PlantID %in% c(a))
dim(root_metadata_undam_root)

undam_targeted_chem_undam_root <- subset(undam_targeted_chem, subset = undam_metadata$PlantID %in% c(a))
dim(undam_targeted_chem_undam_root)
undam_untargeted_chem_undam_root <- subset(undam_untargeted_chem, subset = undam_metadata$PlantID %in% c(a))
dim(undam_untargeted_chem_undam_root)


undam_data_undam_root <- subset(undam_data, subset = undam_metadata$PlantID %in% c(a))
dim(undam_data_undam_root)
undam_metadata_undam_root <- subset(undam_metadata, subset = undam_metadata$PlantID %in% c(a))
dim(undam_metadata_undam_root)


undam_data_undam_root_bray <- vegdist(undam_data_undam_root, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
root_data_undam_root_bray <- vegdist(root_data_undam_root, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
undam_targeted_chem_undam_root_bary <- vegdist(undam_targeted_chem_undam_root, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
undam_untargeted_chem_undam_root_bary <- vegdist(undam_untargeted_chem_undam_root, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence

mantel(undam_data_undam_root_bray, root_data_undam_root_bray)
cor.test(root_metadata_undam_root$Richness, undam_metadata_undam_root$Richness, method = "spearman")

mantel(root_data_undam_root_bray, undam_targeted_chem_undam_root_bary)
mantel(root_data_undam_root_bray, undam_untargeted_chem_undam_root_bary)

adonis2(root_data_undam_root_bray ~ undam_targeted_chem_undam_root$GeniposideHeight)
adonis2(root_data_undam_root_bray ~ undam_targeted_chem_undam_root$GardosideHeight)
adonis2(root_data_undam_root_bray ~ undam_targeted_chem_undam_root$MelittosideHeight)
adonis2(root_data_undam_root_bray ~ undam_targeted_chem_undam_root$X10.acetoxymajorHeight)
adonis2(root_data_undam_root_bray ~ undam_targeted_chem_undam_root$VerbascosideHeight)
adonis2(root_data_undam_root_bray ~ undam_targeted_chem_undam_root$AucubinHeight)


undam_data_undam_root_binary <- (undam_data_undam_root > 0) # logical, or
undam_data_undam_root_binary <- (undam_data_undam_root_binary > 0)*1L # integer 01
undam_data_undam_root_binary <- as.data.frame(undam_data_undam_root_binary)

root_data_undam_root_binary <- (root_data_undam_root > 0) # logical, or
root_data_undam_root_binary <- (root_data_undam_root_binary > 0)*1L # integer 01
root_data_undam_root_binary <- as.data.frame(root_data_undam_root_binary)

both_undam_root_binary <- list()
just_root_data_undam_root_binary <- list()
just_undam_data_undam_root_binary <- list()

for (i in 1:nrow(root_data_undam_root_binary)){
  both_undam_root_binary[[i]] <- length(which(undam_data_undam_root_binary[i,] == 1 & root_data_undam_root_binary[i,] == 1))
  just_undam_data_undam_root_binary[[i]] <- length(which(undam_data_undam_root_binary[i,] == 1 & root_data_undam_root_binary[i,] == 0))
  just_root_data_undam_root_binary[[i]] <- length(which(undam_data_undam_root_binary[i,] == 0 & root_data_undam_root_binary[i,] == 1))
}

both_undam_root_binary <- unlist(both_undam_root_binary, use.names = FALSE)
just_undam_data_undam_root_binary <- unlist(just_undam_data_undam_root_binary, use.names = FALSE)
just_root_data_undam_root_binary <- unlist(just_root_data_undam_root_binary, use.names = FALSE)
dam_root_binary_comparisons <- data.frame(both_undam_root_binary, just_undam_data_undam_root_binary, just_root_data_undam_root_binary)
dam_root_binary_comparisons

dam_root_binary_comparisons_means <- colMeans(dam_root_binary_comparisons)






############Leaf comparisons##############



#Undamaged-damaged comparison
a <- which(table(c(undam_metadata$PlantID, dam_metadata$PlantID)) > 1)
a <- names(a)

dam_data_undam_dam <- subset(dam_data, subset = dam_metadata$PlantID %in% c(a))
dim(dam_data_undam_dam)
dam_targeted_chem_undam_dam <- subset(dam_targeted_chem, subset = dam_metadata$PlantID %in% c(a))
dim(dam_targeted_chem_undam_dam)
dam_untargeted_chem_undam_dam <- subset(dam_untargeted_chem, subset = dam_metadata$PlantID %in% c(a))
dim(dam_untargeted_chem_undam_dam)
dam_metadata_undam_dam <- subset(dam_metadata, subset = dam_metadata$PlantID %in% c(a))
dim(dam_metadata_undam_dam)

dam_metadata_undam_dam <- subset(dam_metadata, subset = dam_metadata$PlantID %in% c(a))
dim(dam_metadata_undam_dam)

undam_data_undam_dam <- subset(undam_data, subset = undam_metadata$PlantID %in% c(a))
dim(undam_data_undam_dam)
undam_targeted_chem_undam_dam <- subset(undam_targeted_chem, subset = undam_metadata$PlantID %in% c(a))
dim(undam_targeted_chem_undam_dam)
undam_untargeted_chem_undam_dam <- subset(undam_untargeted_chem, subset = undam_metadata$PlantID %in% c(a))
dim(undam_untargeted_chem_undam_dam)
undam_metadata_undam_dam <- subset(undam_metadata, subset = undam_metadata$PlantID %in% c(a))
dim(undam_metadata_undam_dam)

undam_data_undam_dam_bray <- vegdist(undam_data_undam_dam, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
dam_data_undam_dam_bray <- vegdist(dam_data_undam_dam, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence


dam_targeted_chem_undam_dam_bray <- vegdist(dam_targeted_chem_undam_dam, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
dam_untargeted_chem_undam_dam_bray <- vegdist(dam_untargeted_chem_undam_dam, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence

undam_targeted_chem_undam_dam_bray <- vegdist(undam_targeted_chem_undam_dam, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence
undam_untargeted_chem_undam_dam_bray <- vegdist(undam_untargeted_chem_undam_dam, Type = "bray", binary = FALSE) # Creates a similarity matrix (Jaccard is the type), based on presence/absence


mantel(undam_data_undam_dam_bray, dam_data_undam_dam_bray)
cor.test(dam_metadata_undam_dam$Richness, undam_metadata_undam_dam$Richness, method = "spearman")

mantel(undam_data_undam_dam_bray, dam_targeted_chem_undam_dam_bray)
mantel(undam_data_undam_dam_bray, dam_untargeted_chem_undam_dam_bray)

mantel(dam_data_undam_dam_bray, undam_targeted_chem_undam_dam_bray)
mantel(dam_data_undam_dam_bray, undam_untargeted_chem_undam_dam_bray)

adonis2(dam_data_undam_dam_bray ~ undam_targeted_chem_undam_dam$GeniposideHeight)
adonis2(dam_data_undam_dam_bray ~ undam_targeted_chem_undam_dam$GardosideHeight)
adonis2(dam_data_undam_dam_bray ~ undam_targeted_chem_undam_dam$MelittosideHeight)
adonis2(dam_data_undam_dam_bray ~ undam_targeted_chem_undam_dam$X10.acetoxymajorHeight)
adonis2(dam_data_undam_dam_bray ~ undam_targeted_chem_undam_dam$VerbascosideHeight)
adonis2(dam_data_undam_dam_bray ~ undam_targeted_chem_undam_dam$AucubinHeight)

adonis2(undam_data_undam_dam_bray ~ dam_targeted_chem_undam_dam$GeniposideHeight)
adonis2(undam_data_undam_dam_bray ~ dam_targeted_chem_undam_dam$GardosideHeight)
adonis2(undam_data_undam_dam_bray ~ dam_targeted_chem_undam_dam$MelittosideHeight)
adonis2(undam_data_undam_dam_bray ~ dam_targeted_chem_undam_dam$X10.acetoxymajorHeight)
adonis2(undam_data_undam_dam_bray ~ dam_targeted_chem_undam_dam$VerbascosideHeight)
adonis2(undam_data_undam_dam_bray ~ dam_targeted_chem_undam_dam$AucubinHeight)


plot(dam_targeted_chem_undam_dam$VerbascosideHeight, undam_targeted_chem_undam_dam$VerbascosideHeight)

undam_data_undam_dam_binary <- (undam_data_undam_dam > 0) # logical, or
undam_data_undam_dam_binary <- (undam_data_undam_dam_binary > 0)*1L # integer 01
undam_data_undam_dam_binary <- as.data.frame(undam_data_undam_dam_binary)

dam_data_undam_dam_binary <- (dam_data_undam_dam > 0) # logical, or
dam_data_undam_dam_binary <- (dam_data_undam_dam_binary > 0)*1L # integer 01
dam_data_undam_dam_binary <- as.data.frame(dam_data_undam_dam_binary)

both_undam_dam_binary <- list()
just_dam_data_undam_dam_binary <- list()
just_undam_data_undam_dam_binary <- list()

for (i in 1:nrow(dam_data_undam_dam_binary)){
  both_undam_dam_binary[[i]] <- length(which(undam_data_undam_dam_binary[i,] == 1 & dam_data_undam_dam_binary[i,] == 1))
  just_undam_data_undam_dam_binary[[i]] <- length(which(undam_data_undam_dam_binary[i,] == 1 & dam_data_undam_dam_binary[i,] == 0))
  just_dam_data_undam_dam_binary[[i]] <- length(which(undam_data_undam_dam_binary[i,] == 0 & dam_data_undam_dam_binary[i,] == 1))
}

both_undam_dam_binary <- unlist(both_undam_dam_binary, use.names = FALSE)
just_undam_data_undam_dam_binary <- unlist(just_undam_data_undam_dam_binary, use.names = FALSE)
just_dam_data_undam_dam_binary <- unlist(just_dam_data_undam_dam_binary, use.names = FALSE)
dam_root_binary_comparisons <- data.frame(both_undam_dam_binary, just_undam_data_undam_dam_binary, just_dam_data_undam_dam_binary)
dam_root_binary_comparisons

dam_root_binary_comparisons_means <- colMeans(dam_root_binary_comparisons)
dam_root_binary_comparisons_means













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
long_jaccards <- stack(jaccards)$values

long_data <- data.frame("Stratum" = long_status, "PlantID" = long_PlantID, "jaccard" = long_jaccards)
long_data$Stratum_simplified <- str_replace(long_data$Stratum, "Damaged leaf", "Leaf")
long_data$Stratum_simplified <- str_replace(long_data$Stratum_simplified, "Undamaged leaf", "Leaf")

jaccard_Stratums <- dist_groups(jaccards, metadata$Stratum)
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


ggplot(Big_data, aes(x = Stratum, y = Jaccards, 
                     fill = Group)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_grid(~ Barrier, scales = 'free', space = 'free')

ggplot(Big_data, aes(x = Stratum, y = Similarity, 
                     fill = Group)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_grid(~ Barrier, scales = 'free', space = 'free')

ggplot(Big_data, aes(x = Stratum, y = pH, 
                     fill = Group)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_grid(~ Group2, scales = 'free', space = 'free')


Big_data$Group <- factor(Big_data$Group, levels = c("Within plants", "Within sites", "Between sites"))

ggplot(subset(Big_data, Big_data$Barrier == "Above-belowground transfer"), 
  aes(x = Group, y = Jaccards, fill = Stratum_simplified))  +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_boxplot(alpha = 0.8) + labs(x = NULL, y = "Shared ASVs", fill = NULL) + 
  scale_fill_brewer(palette = "Accent") #+ facet_grid(metadata$Plant ~ metadata$Site) 

saveRDS(Big_data, file = "Bacterial_big_data_no_mito.RDS")

class(Big_data$EC)
Big_data$Status

ggplot(Big_data, aes(x = Status, y = Jaccards, 
                     fill = Group)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(Big_data, aes(x = Group, y = EC, 
                     fill = Group)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggplot(Big_data, aes(x = Group, y = pH, 
                     fill = Group)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))




Big_data_within_plant <- Big_data [which(Big_data$Group == "Within plants"), ]
Big_data_within_plant_interactions <- Big_data_within_plant[which (Big_data_within_plant$Stratum_simplified != "Within Leaf"),]

ggplot(Big_data_within_plant_interactions, aes(x = Stratum_simplified, y = Jaccards, 
                     fill = Group)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


metadata$Stratum <- factor (Big_data_within_plant_interactions$Stratum_simplified, 
                            levels = c("U", "D", "R", "S"), 
                            labels = c("Undamaged leaf", "Damaged leaf", "Root", "Soil"))





sum(data)
mean(rowSums(data))

count(metadata$Status)

summary(rowSums(subset(data, metadata$Status == "R")))

x <- readRDS(file = "metadata_data_no_mito.RDS")
count(x$Stratum)
