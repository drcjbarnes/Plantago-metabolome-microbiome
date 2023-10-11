
library(tidyr)
library(tidyverse)
library(reshape2)
library(dplyr)
library(plyr)
library(tibble)
library(taRifx)
library(stringr)
library(vegan)
library(lme4)
library(dplyr)
library(MASS)
library(RColorBrewer)
library(ggpubr)
library(ggplot2); theme_set(theme_bw(base_size = 16))
library(RColorBrewer)
library(corrplot)
library(patchwork)
library(ape)
library(phytools)
library(caper)
library(picante)
library(funrar) #Converting into relative abundances 
library(stringr)
library(usedist)
library("pracma")
library("FUNGuildR")
library(FSA)

library(mixOmics)
library(mixKernel)
library(igraph)
library(qgraph)

setwd("/Users/christopherbarnes/Dropbox/Danish_Plantago_microbiomes/Fine_root_endophytes_and_AMF")

FRE_data <- read.table("/Users/christopherbarnes/Dropbox/Danish_Plantago_microbiomes/Fungi/FRE_data.txt", sep = "\t", header = TRUE, row.names = 1)
FRE_data <- FRE_data[-c(68),]
dim(FRE_data)
FRE_data <- as.matrix(FRE_data)
FRE_data <- nthroot(FRE_data, 4)
rel_FRE_data <- (make_relative(FRE_data)*100)
rel_FRE_data <- FRE_data

FRE_metadata <- read.table("/Users/christopherbarnes/Dropbox/Danish_Plantago_microbiomes/Fungi/FRE_metadata_initial.txt", sep = "\t", header = TRUE, row.names = 1)
FRE_metadata <- FRE_metadata[-c(68),]

rel_FRE_data <- rel_FRE_data[which(FRE_metadata$Status == "Root"),]
dim(rel_FRE_data)
FRE_metadata <- FRE_metadata[which(FRE_metadata$Status == "Root"),]
dim(FRE_metadata)

FRE_metadata$merging_var <- paste(FRE_metadata$PlantID, FRE_metadata$Status, sep = "_")
rownames(rel_FRE_data) <- FRE_metadata$merging_var
rel_FRE_data <- rel_FRE_data[order(rownames(rel_FRE_data)),]
FRE_metadata <- FRE_metadata[order(FRE_metadata$merging_var),]

targeted_chem_mean_leaves <- readRDS(file = "/Users/christopherbarnes/Dropbox/Danish_Plantago_microbiomes/Fungi/targeted_chem_mean_leaves.RDS")
untargeted_chem_mean_leaves <- readRDS(file = "/Users/christopherbarnes/Dropbox/Danish_Plantago_microbiomes/Fungi/untargeted_chem_mean_leaves.RDS")
metadata_mean_leaves <- readRDS(file = "/Users/christopherbarnes/Dropbox/Danish_Plantago_microbiomes/Fungi/metadata_mean_leaves.RDS")

root_metadata <- subset(metadata_mean_leaves, metadata_mean_leaves$Stratum_simplified == "Root")
root_untargeted_chem <- subset(untargeted_chem_mean_leaves, metadata_mean_leaves$Stratum_simplified == "Leaf")
root_targeted_chem <- subset(targeted_chem_mean_leaves, metadata_mean_leaves$Stratum_simplified == "Leaf")

#FRE_Taxonomy <- read.table("/Users/christopherbarnes/Dropbox/Danish_Plantago_microbiomes/Fungi/FRE_Taxonomy_curated.txt", sep = '\t', header = TRUE, row.names = 1)
#FRE_Guilds <- funguild_assign(FRE_Taxonomy)
#colnames(FRE_Taxonomy) <- "taxa"
#FRE_Taxonomy <- str_split_fixed(FRE_Taxonomy$taxa, ";", 11)
#rownames(FRE_Taxonomy) <- colnames(FRE_data)
#FRE_Taxonomy <- as.data.frame(FRE_Taxonomy)
#dim(FRE_Taxonomy)
#FRE_Taxonomy <- data.frame(FRE_Taxonomy, FRE_Guilds[,-c(1)])
#saveRDS(FRE_Taxonomy, file = "FRE_Taxonomy.RDS")
FRE_Taxonomy <- readRDS(file = "FRE_Taxonomy.RDS")
FRE_Taxonomy <- FRE_Taxonomy[,-c(1:4, 7)]
colnames(FRE_Taxonomy)[1:6] <- c("Phylum", "Subphylum", "Order", "Family", "Genus", "Species")

empty_ASVs <- (colSums(rel_FRE_data != 0)) > 0 #OTU as columns
rel_FRE_data <- rel_FRE_data[, empty_ASVs == "TRUE"] #Removes taxa that are all 0s
dim(rel_FRE_data)
FRE_Taxonomy <- FRE_Taxonomy[ empty_ASVs == "TRUE",] #Removes taxa that are all 0s
dim(FRE_Taxonomy)
rowSums(rel_FRE_data)

x <- which(FRE_Taxonomy$Subphylum == " Glomeromycotina" | FRE_Taxonomy$Subphylum == " Mucoromycotina")
rel_FRE_data <- rel_FRE_data[, x] #Removes taxa that are all 0s
dim(rel_FRE_data)
FRE_Taxonomy <- FRE_Taxonomy[ x ,] #Removes taxa that are all 0s
dim(FRE_Taxonomy)
rowSums(rel_FRE_data)

rel_FRE_data <- as.matrix(rel_FRE_data)
rel_FRE_data <- (make_relative(rel_FRE_data)*100)
rel_FRE_data <- as.data.frame(rel_FRE_data)



Glom.taxa <- which(FRE_Taxonomy$Subphylum == " Glomeromycotina")
Glom.data <- rel_FRE_data[, c(Glom.taxa)]
Glom.taxonomy <- FRE_Taxonomy[Glom.taxa , ]
dim(Glom.taxonomy)
dim(Glom.data)

Glom_totals <- rowSums(Glom.data)
Glom_Richness <- rowSums(Glom.data > 0)
mean(Glom_totals)
sd(Glom_totals)
min(rowSums(Glom.data))
sort(Glom_totals)
mean(Glom_totals)
sd(Glom_totals)
mean(Glom_Richness)
sd(Glom_Richness)

Muc.taxa <- which(FRE_Taxonomy$Subphylum == " Mucoromycotina")
Muc_data <- rel_FRE_data[, c(Muc.taxa)]
Muc_taxonomy <- FRE_Taxonomy[Muc.taxa , ]
dim(Muc_data)
dim(Muc_taxonomy)
Muc_totals <- rowSums(Muc_data)
Muc_richness <- rowSums(Muc_data > 0)
mean(Muc_totals)
sd(Muc_totals)
mean(Muc_richness)
sd(Muc_richness)

Wide.data <- data.frame("Site" = FRE_metadata$Site, 
                        "Habitat" = FRE_metadata$Habitat, "pH" = FRE_metadata$pH, 
                        "EC" = FRE_metadata$EC, "BulkDensity" = FRE_metadata$BulkDensity, "Geography3" = FRE_metadata$Geography3,
                        "PlantID" = FRE_metadata$PlantID,
                        "Muc_totals" = Muc_totals,
                        "Glom_totals" = Glom_totals)

long.data <- melt(Wide.data, id.vars = c("Site", "Habitat", "pH", "EC", "BulkDensity", "Geography3", "PlantID"), 
                  variable.name = "Guild", value.name = "Abundance")


FigureS1 <- ggplot(long.data, aes(x = PlantID, y = Abundance, fill = Guild)) +
  geom_bar(stat = "identity", colour="black") + facet_grid(~ Habitat, scales = "free", space = "free_x") +
  theme(strip.text.x = element_text(angle = 0)) + theme(text = element_text(size=12)) +
  theme(axis.text.x = element_text(angle = 90, size = 8))

colourCount = length(unique(long.order.data$Order))
#getPalette = colorRampPalette(brewer.pal(4, "Accent"))

FigureS1 <- FigureS1 + scale_fill_manual(values = getPalette(colourCount))
FigureS1 #Shows you our graph


Glom.bray <- vegdist(Glom.data, Type = 'Bray')

root_untargeted_chem_bray <- vegdist(root_untargeted_chem, Type = 'Bray')
root_targeted_chem_bray <- vegdist(root_targeted_chem, Type = 'Bray')
mantel(root_untargeted_chem_bray, root_targeted_chem_bray)




#######Final AMF analyses############
#adonis2(Glom.bray ~ root_metadata$Site)
adonis2(Glom.bray ~ root_metadata$Habitat)
#adonis2(Glom.bray ~ root_metadata$Light)

adonis2(Glom.bray ~ root_metadata$AnnPrec)
adonis2(Glom.bray ~ root_metadata$ColdestQuartMeanTemp)
adonis2(Glom.bray ~ root_metadata$MeanAnnTemp)

adonis2(Glom.bray ~ root_metadata$pH)
adonis2(Glom.bray ~ root_metadata$EC)
adonis2(Glom.bray ~ root_metadata$BulkDensity)

#adonis2(Glom.bray ~ root_metadata$PCNM1)
#adonis2(Glom.bray ~ root_metadata$PCNM2)
adonis2(Glom.bray ~ root_metadata$Island)
adonis2(Glom.bray ~ root_metadata$Geography3)

adonis2(Glom.bray ~ root_metadata$Habitat + root_metadata$pH +
          root_metadata$BulkDensity, group = root_metadata$Site)


kruskal.test(Glom_Richness ~ root_metadata$Site)
kruskal.test(Glom_Richness ~ root_metadata$Habitat)
#kruskal.test(FRE_Richness ~ root_metadata$Light)

cor.test(Glom_Richness, root_metadata$AnnPrec, method = "spearman")
cor.test(Glom_Richness, root_metadata$ColdestQuartMeanTemp, method = "spearman")
cor.test(Glom_Richness, root_metadata$MeanAnnTemp, method = "spearman")

cor.test(Glom_Richness, root_metadata$pH, method = "spearman")
cor.test(Glom_Richness, root_metadata$EC, method = "spearman")
cor.test(Glom_Richness, root_metadata$BulkDensity, method = "spearman")

cor.test(Glom_Richness, root_metadata$PCNM1, method = "spearman")
cor.test(Glom_Richness, root_metadata$PCNM2, method = "spearman")
kruskal.test(Glom_Richness ~ root_metadata$Island)
kruskal.test(Glom_Richness ~ root_metadata$Geography3)



mantel(Glom.bray, root_untargeted_chem_bray)
mantel(Glom.bray, root_targeted_chem_bray)

adonis2(Glom.bray ~ root_targeted_chem$GeniposideHeight)
adonis2(Glom.bray ~ root_targeted_chem$GardosideHeight)
adonis2(Glom.bray ~ root_targeted_chem$MelittosideHeight)
adonis2(Glom.bray ~ root_targeted_chem$X10.acetoxymajorHeight)
adonis2(Glom.bray ~ root_targeted_chem$VerbascosideHeight)
adonis2(Glom.bray ~ root_targeted_chem$AucubinHeight)



###########Final FRE Analyses################
length(c(which(rowSums(Muc_data) == 0)))
filtered_muc_targeted_chem <- root_targeted_chem[ -c(which(rowSums(Muc_data) == 0)), ]
filtered_muc_untargeted_chem <- root_untargeted_chem[ -c(which(rowSums(Muc_data) == 0)), ]
filtered_muc_data <- Muc_data[ -c(which(rowSums(Muc_data) == 0)), ]
filtered_muc_bray  <- vegdist(filtered_muc_data, Type = 'Bray')
filtered_muc_richness <- Muc_richness[-c(which(rowSums(Muc_data) == 0))]
filtered_muc_metadata <- root_metadata[ -c(which(rowSums(Muc_data) == 0)), ]

#adonis2(filtered_muc_bray ~ filtered_muc_metadata$Site)
adonis2(filtered_muc_bray ~ filtered_muc_metadata$Habitat)
#adonis2(filtered_muc_bray ~ filtered_muc_metadata$Light)

adonis2(filtered_muc_bray ~ filtered_muc_metadata$AnnPrec)
adonis2(filtered_muc_bray ~ filtered_muc_metadata$ColdestQuartMeanTemp)
adonis2(filtered_muc_bray ~ filtered_muc_metadata$MeanAnnTemp)

adonis2(filtered_muc_bray ~ filtered_muc_metadata$pH)
adonis2(filtered_muc_bray ~ filtered_muc_metadata$EC)
adonis2(filtered_muc_bray ~ filtered_muc_metadata$BulkDensity)

#adonis2(filtered_muc_bray ~ filtered_muc_metadata$PCNM1)
#adonis2(filtered_muc_bray ~ filtered_muc_metadata$PCNM2)
adonis2(filtered_muc_bray ~ filtered_muc_metadata$Island)
adonis2(filtered_muc_bray ~ filtered_muc_metadata$Geography3)

adonis2(filtered_muc_bray ~ filtered_muc_metadata$Habitat + filtered_muc_metadata$Island +
  filtered_muc_metadata$pH, group = filtered_muc_metadata$Site)


kruskal.test(filtered_muc_richness ~ filtered_muc_metadata$Site)
kruskal.test(filtered_muc_richness ~ filtered_muc_metadata$Habitat)
#kruskal.test(Muc_Richness ~ filtered_muc_metadata$Light)

cor.test(filtered_muc_richness, filtered_muc_metadata$AnnPrec, method = "spearman")
cor.test(filtered_muc_richness, filtered_muc_metadata$ColdestQuartMeanTemp, method = "spearman")
cor.test(filtered_muc_richness, filtered_muc_metadata$MeanAnnTemp, method = "spearman")

cor.test(filtered_muc_richness, filtered_muc_metadata$pH, method = "spearman")
cor.test(filtered_muc_richness, filtered_muc_metadata$EC, method = "spearman")
cor.test(filtered_muc_richness, filtered_muc_metadata$BulkDensity, method = "spearman")

cor.test(filtered_muc_richness, filtered_muc_metadata$PCNM1, method = "spearman")
cor.test(filtered_muc_richness, filtered_muc_metadata$PCNM2, method = "spearman")
kruskal.test(filtered_muc_richness ~ filtered_muc_metadata$Island)
kruskal.test(filtered_muc_richness ~ filtered_muc_metadata$Geography3)


filtered_muc_untargeted_chem_bray  <- vegdist(filtered_muc_untargeted_chem, Type = 'Bray')
mantel(filtered_muc_bray, filtered_muc_untargeted_chem_bray, na.rm = TRUE)

filtered_muc_targeted_chem_bray  <- vegdist(filtered_muc_targeted_chem, Type = 'Bray')
mantel(filtered_muc_bray, filtered_muc_targeted_chem_bray, na.rm = TRUE)

adonis2(filtered_muc_bray ~ filtered_muc_targeted_chem$GeniposideHeight)
adonis2(filtered_muc_bray ~ filtered_muc_targeted_chem$GardosideHeight)
adonis2(filtered_muc_bray ~ filtered_muc_targeted_chem$MelittosideHeight)
adonis2(filtered_muc_bray ~ filtered_muc_targeted_chem$X10.acetoxymajorHeight)
adonis2(filtered_muc_bray ~ filtered_muc_targeted_chem$VerbascosideHeight)
adonis2(filtered_muc_bray ~ filtered_muc_targeted_chem$AucubinHeight)


















rownames(Glom.taxonomy)
colnames(Glom.data)

###########################Network analyses#################################
colnames(rel_data_mean_leaves) <- gsub("SeqVar", "ASV", colnames(rel_data_mean_leaves))
genus_descriptors <- paste(Taxonomy$Genus, colnames(rel_data_mean_leaves), sep = "_")
colnames(leaf_data) <- genus_descriptors
rownames(Taxonomy) <- colnames(leaf_data)


Glom.taxonomy <- FRE_Taxonomy[Glom.taxa , ]
genus_descriptors <- paste(Glom.taxonomy$Genus, rownames(Glom.taxonomy), sep = "_")
genus_descriptors <- gsub("FRE_AMF", "aOTU", genus_descriptors)

colnames(Glom.data) <- genus_descriptors
rownames(Glom.taxonomy) <- colnames(Glom.data)

Glom.data.filtered <- Glom.data[,c(which(colSums(Glom.data > 0) > 4))]

x <- nearZeroVar(root_untargeted_chem)
root_untargeted_chem_filtered <- root_untargeted_chem[,-c(x$Position)]
dim(root_untargeted_chem_filtered)

rownames(root_untargeted_chem_filtered) <- rownames(Glom.data)
rownames(root_targeted_chem) <- rownames(Glom.data)

imgCor(Glom.data.filtered, root_untargeted_chem_filtered, sideColors = c("purple", "green")) 
imgCor(Glom.data.filtered, root_targeted_chem, sideColors = c("purple", "green")) 

# set grid search values for each regularisation parameter
grid1 <- seq(0.001, 0.2, length = 10) 
grid2 <- seq(0.001, 0.2, length = 10)


#root_data_untar_rcc_tune <- tune.rcc(Glom.data.filtered, root_untargeted_chem_filtered, grid1 = grid1, grid2 = grid2, validation = "loo") 
#root_data_untar_rcc_tune <- saveRDS(root_data_untar_rcc_tune, file = "root_data_untar_rcc_tune")
#root_data_untar_rcc_tune <- rcc(root_data_filtered, root_untargeted_chem_filtered, method = 'shrinkage') 
root_data_untar_rcc_tune <- readRDS(file = "root_data_untar_rcc_tune")

#root_data_tar_rcc_tune <- tune.rcc(Glom.data.filtered, root_targeted_chem, grid1 = grid1, grid2 = grid2, validation = "loo") 
#root_data_tar_rcc_tune <- saveRDS(root_data_tar_rcc_tune, file = "root_data_tar_rcc_tune")
root_data_tar_rcc_tune <- readRDS(file = "root_data_tar_rcc_tune")



# formed optimised CV rCCA. Need to adjust for this for all four comparisons. 
root_data_untar_rcc_tune_opt.l1 <- root_data_untar_rcc_tune$opt.lambda1 # extract the optimal lambda values
root_data_untar_rcc_tune_opt.l2 <- root_data_untar_rcc_tune$opt.lambda2

root_data_untar_CV_rcc_tune <- rcc(Glom.data.filtered, root_untargeted_chem_filtered, method = "ridge", lambda1 = root_data_untar_rcc_tune_opt.l1, lambda2 = root_data_untar_rcc_tune_opt.l2) 

colnames(root_data_untar_CV_rcc_tune$X) <- colnames(Glom.data.filtered)
rownames(root_data_untar_CV_rcc_tune$loadings$X) <- colnames(Glom.data.filtered)
root_data_untar_CV_rcc_tune$names$colnames$X <- colnames(Glom.data.filtered)

colnames(root_data_untar_CV_rcc_tune$Y) <- paste("Untar_metab", seq(1: length(colnames(root_data_untar_CV_rcc_tune$Y))), sep = "_")
rownames(root_data_untar_CV_rcc_tune$loadings$Y) <- paste("Untar_metab", seq(1: length(colnames(root_data_untar_CV_rcc_tune$Y))), sep = "_")
root_data_untar_CV_rcc_tune$names$colnames$Y <- paste("Untar_metab", seq(1: length(colnames(root_data_untar_CV_rcc_tune$Y))), sep = "_")

plot(root_data_untar_CV_rcc_tune, type = "barplot", main = "Cross Validation") 

plotIndiv(root_data_untar_CV_rcc_tune, comp = 1:2, 
          ind.names = root_metadata$Habitat,
          group = root_metadata$Habitat, rep.space = "XY-variate", 
          legend = TRUE, title = '(a) Nutrimouse, rCCA CV XY-space')


plotArrow(root_data_untar_CV_rcc_tune, group = root_metadata$Habitat, 
          col.per.group = color.mixo(1:5),
          title = '(a) Nutrimouse, CV method')


plotVar(root_data_untar_CV_rcc_tune, var.names = c(TRUE, TRUE),
        cex = c(4, 4), cutoff = 0.3,
        title = '(a) Nutrimouse, rCCA CV comp 1 - 2')

network_root_data_untar_CV_rcc_tune <- network(root_data_untar_CV_rcc_tune, comp = 1:2, interactive = FALSE, lwd.edge = 1, cutoff = 0.30, symkey = FALSE, block.var.names = FALSE, alpha.node = 0.85)
#write.graph(network_root_data_untar_CV_rcc_tune$gR, file = "network_root_data_untar_CV_rcc_tune.glm")
#write.graph(network_root_data_untar_CV_rcc_tune$gR, fil = "network_root_data_untar_CV_rcc_tune.glm")
cim(root_data_untar_CV_rcc_tune, comp = 1:2, cutoff = 0.3, margins = c(10, 17)) 


root_untar_network <- network_root_data_untar_CV_rcc_tune$gR
class(root_untar_network)
plot(root_untar_network)

saveRDS(root_untar_network, file = "root_untar_network.RDS")


root_untar_network_matrix <- network_root_data_untar_CV_rcc_tune$M
dim(root_untar_network_matrix)

Glom.taxonomy[c(which(row.names(Glom.taxonomy) == "FRE_AMF11")), ] 
Glom.taxonomy[c(which(row.names(Glom.taxonomy) == "FRE_AMF18")), ] 
Glom.taxonomy[c(which(row.names(Glom.taxonomy) == "FRE_AMF22")), ] 
Glom.taxonomy[c(which(row.names(Glom.taxonomy) == "FRE_AMF53")), ] 
Glom.taxonomy[c(which(row.names(Glom.taxonomy) == "FRE_AMF94")), ] 
Glom.taxonomy[c(which(row.names(Glom.taxonomy) == "FRE_AMF126")), ] 



root_untar_fungal_taxa <- 
  which(row.names(Glom.taxonomy) == "FRE_AMF11" | row.names(Glom.taxonomy) == "FRE_AMF18" | row.names(Glom.taxonomy) == "FRE_AMF22"
        | row.names(Glom.taxonomy) == "FRE_AMF53" | row.names(Glom.taxonomy) == "FRE_AMF94" | row.names(Glom.taxonomy) == "FRE_AMF126")

#write.table (Glom.taxonomy[c(root_untar_fungal_taxa) , ], file = "root_endophytes_untar_fungal_taxa.txt")


key_taxa <- c("FRE_AMF11", "FRE_AMF18", "FRE_AMF22", "FRE_AMF53", "FRE_AMF94", "FRE_AMF126")


# formed optimised CV rCCA. Need to adjust for this for all four comparisons. 
root_data_tar_rcc_tune_opt.l1 <- root_data_tar_rcc_tune$opt.lambda1 # extract the optimal lambda values
root_data_tar_rcc_tune_opt.l2 <- root_data_tar_rcc_tune$opt.lambda2

root_data_tar_CV_rcc_tune <- rcc(Glom.data, root_targeted_chem, method = "ridge", lambda1 = root_data_tar_rcc_tune_opt.l1, lambda2 = root_data_tar_rcc_tune_opt.l2) 
plot(root_data_tar_CV_rcc_tune, type = "barplot", main = "Cross Validation") 

plotIndiv(root_data_tar_CV_rcc_tune, comp = 1:2, 
          ind.names = root_metadata$Habitat,
          group = root_metadata$Habitat, rep.space = "XY-variate", 
          legend = TRUE, title = '(a) Nutrimouse, rCCA CV XY-space')

plotArrow(root_data_tar_CV_rcc_tune, group = root_metadata$Habitat, 
          col.per.group = color.mixo(1:5),
          title = '(a) Nutrimouse, CV method')

plotVar(root_data_tar_CV_rcc_tune, var.names = c(TRUE, TRUE),
        cex = c(4, 4), cutoff = 0.35,
        title = '(a) Nutrimouse, rCCA CV comp 1 - 2')

network_root_data_tar_CV_rcc_tune <- network(root_data_tar_CV_rcc_tune, comp = 1:2, interactive = FALSE, lwd.edge = 1, cutoff = 0.30, symkey = FALSE, block.var.names = FALSE, alpha.node = 0.85)
plot(igraph)

saveRDS(network_root_data_tar_CV_rcc_tune$gR, file = "root_tar_network.RDS")

write.graph(network_root_data_tar_CV_rcc_tune$gR, "Targeted_fungi_undamaged.graphml", format = "graphml")
#write.graph(network_root_data_tar_CV_rcc_tune$gR, fil = "network_root_data_tar_CV_rcc_tune.glm")
#write.graph(network_root_data_tar_CV_rcc_tune$gR, fil = "network_root_data_tar_CV_rcc_tune.glm")
a <- cim(root_data_tar_CV_rcc_tune, comp = 1:2, cutoff = 0.3, margins = c(10, 17))

a

class(a)
max(a$mat)



a$mat[a$mat < 0.3] <- 0
a$mat[a$mat < -0.3] <- 0

b <- which(a$mat > 0.3 | a$mat < -0.3)
a$mat[-b] <- 0



#Taxonomy[c(which(row.names(Taxonomy) == "SeqVar138")), ] #Too poorly assigned


root_tar_fungal_taxa <- 
  which(row.names(Taxonomy) == "SeqVar138" | row.names(Taxonomy) == "SeqVar151" | row.names(Taxonomy) == "SeqVar549"
        | row.names(Taxonomy) == "SeqVar910" | row.names(Taxonomy) == "SeqVar2342" | row.names(Taxonomy) == "SeqVar2979"
        | row.names(Taxonomy) == "SeqVar3327" | row.names(Taxonomy) == "SeqVar3325" | row.names(Taxonomy) == "SeqVar3331"
        | row.names(Taxonomy) == "SeqVar3337" | row.names(Taxonomy) == "SeqVar3660" | row.names(Taxonomy) == "SeqVar4568")

#write.table (Taxonomy[c(root_tar_fungal_taxa) , ], file = "root_tar_fungal_taxa.txt")


#plot(root_data_filtered[(which(colnames(root_data_filtered) == "SeqVar1721")), ])

key_taxa_untar <- c("FRE_AMF11", "FRE_AMF18", "FRE_AMF22", "FRE_AMF53", "FRE_AMF94", "FRE_AMF126")

key_taxa_tar <- c("FRE_AMF60", "FRE_AMF317", "FRE_AMF389", "FRE_AMF434", "FRE_AMF464", "FRE_AMF745", "FRE_AMF900", 
                  "FRE_AMF986", "FRE_AMF1063", "FRE_AMF1115")



root_untar_fungal_taxa <- 
  which(row.names(Glom.taxonomy) == "FRE_AMF11" | row.names(Glom.taxonomy) == "FRE_AMF18" | row.names(Glom.taxonomy) == "FRE_AMF22"
        | row.names(Glom.taxonomy) == "FRE_AMF53" | row.names(Glom.taxonomy) == "FRE_AMF94" | row.names(Glom.taxonomy) == "FRE_AMF126")

#write.table (Glom.taxonomy[c(key_taxa_tar) , ], file = "root_endophytes_tar_fungal_taxa.txt")



























par(mfrow=c(1,2)) 

fungi_untar_network <- readRDS(file = "root_untar_network.RDS")
#plot(fungi_untar_network, weighted=TRUE, mode="lower")

#coords <- layout_(fungi_untar_network, as_star())
#plot(fungi_untar_network, layout = coords)
#coords <- layout_(fungi_untar_network, as_star())
#plot(fungi_untar_network, layout = coords)

V(fungi_untar_network)$label
V(fungi_untar_network)$label[[1]] <- "Glomeromycotina_aOTU11"
#plot(fungi_untar_network, weighted = TRUE, mode="lower")

V(fungi_untar_network)$color <- c(rep("cyan", times = 6), rep("light green", times = 19))
V(fungi_untar_network)$shape <- c(rep("circle", times = 6), rep("circle", times = 19))
V(fungi_untar_network)$shape <- rep("circle", times = 25)
V(fungi_untar_network)$size <- rep(10, times = 25)
#plot(fungi_untar_network, layout = layout.circle, vertex.label.dist=1.5)

#plot(fungi_untar_network, layout = layout.circle, vertex.label.dist=1.5)
E(fungi_untar_network)$color <- ifelse(E(fungi_untar_network)$weight  > 0, "black", "red")
#plot(fungi_untar_network, layout = layout.circle, vertex.label.dist=1.5)

V(fungi_untar_network)$label.cex = 1.2

plot(fungi_untar_network, layout = layout.circle, vertex.label.dist=1.5)










fungi_tar_network <- readRDS(file = "root_tar_network.RDS")
#plot(fungi_tar_network, weighted=TRUE, mode="lower")



V(fungi_tar_network)$label
V(fungi_tar_network)$label[[1]] <- "Glomeromycotina_aOTU6"
V(fungi_tar_network)$label[[4]] <- "Glomerales_aOTU434"
V(fungi_tar_network)$label[[5]] <- "Glomerales_aOTU464"
V(fungi_tar_network)$label[[6]] <- "Glomerales_aOTU745"
V(fungi_tar_network)$label[[7]] <- "Glomerales_aOTU900"
V(fungi_tar_network)$label[[11]] <- "Geniposide"
V(fungi_tar_network)$label[[12]] <- "Gardoside"
V(fungi_tar_network)$label[[13]] <- "Melittoside"
V(fungi_tar_network)$label[[14]] <- "Aucubin"
V(fungi_tar_network)$label
#plot(fungi_tar_network, weighted = TRUE, mode="lower")

V(fungi_tar_network)$color <- c(rep("cyan", times = 10), rep("green3", times = 4))
V(fungi_tar_network)$shape <- c(rep("circle", times = 10), rep("circle", times = 4))
V(fungi_tar_network)$shape <- rep("circle", times = 14)
V(fungi_tar_network)$size <- rep(10, times = 14)
#plot(fungi_tar_network, layout = layout.circle, vertex.label.dist=1.5)

#plot(fungi_tar_network, layout = layout.circle, vertex.label.dist=1.5)
E(fungi_tar_network)$color <- ifelse(E(fungi_tar_network)$weight  > 0, "black", "red")
#plot(fungi_tar_network, layout = layout.circle, vertex.label.dist=1.5)

V(fungi_tar_network)$label.cex = 1.2
plot(fungi_tar_network, layout = layout.circle, vertex.label.dist=1.5)

