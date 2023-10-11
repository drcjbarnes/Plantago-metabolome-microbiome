

library(msa)
library(ape)
library(ips)
library(phangorn)
library(DECIPHER)
library("seqinr")
library(phytools)



big_bacteria_fasta <- read.fasta(file = "dada_seqs_editted.txt") #FASTA file with all seqs in.
#big_bacteria_fasta <- read.fasta(file = "/Users/christopherbarnes/Dropbox/Danish_Plantago_microbiomes/bacteria/DADA2_files/dada_seqs_editted.txt") #FASTA file with all seqs in.

rel_data <- readRDS(file = "data_data_no_mito.RDS")
taxonomy <- readRDS(file = "final_taxonomy_no_mito.RDS")
metadata <- readRDS(file = "final_metadata.RDS")

Occurrences <- colSums(rel_data > 0)
a <- sort(Occurrences, decreasing = TRUE)
names(a)
b <- rel_data[,names(a)]
rel_data_reduced <- b[,1:1000]
dim(rel_data_reduced)

#Interesting_bacterial_OTU_fasta <- big_bacteria_fasta[names(big_bacteria_fasta) %in% colnames(rel_data_reduced)] #Extract interesting OTU from the big fasta.
Interesting_bacterial_OTU_fasta <- big_bacteria_fasta[names(big_bacteria_fasta) %in% colnames(rel_data)] #Extract interesting OTU from the big fasta.
length(Interesting_bacterial_OTU_fasta)

Interesting_bacterial_OTU_fasta <- Interesting_bacterial_OTU_fasta[order(names(Interesting_bacterial_OTU_fasta))] #Ensuring they are in the same order
names(Interesting_bacterial_OTU_fasta)

rel_data_reduced <- rel_data_reduced[,order(colnames(rel_data_reduced))] #Ensuring they are in the same order
colnames(rel_data_reduced)

#write.fasta(sequences = Interesting_bacterial_OTU_fasta, names = names(Interesting_bacterial_OTU_fasta), file.out = "Interesting_bacterial_OTU_fasta.fasta")
#saveRDS(rel_data_reduced, file = "rel_data_reduced.RDS")


#Interesting_bacterial_OTU_fasta <- readDNAStringSet(file = "Interesting_bacterial_OTU_fasta.fasta")
#Reading as fasta file corrects the format for now. Need to convert in script, probsbly with MSA. 

Interesting_bacterial_OTU_fasta_aligned <- msa(Interesting_bacterial_OTU_fasta, method = "ClustalW", type = "dna", order = "input")

alignment2Fasta <- function(alignment, filename) {
  sink(filename)
  
  n <- length(rownames(alignment))
  for(i in seq(1, n)) {
    cat(paste0('>', rownames(alignment)[i]))
    cat('\n')
    the.sequence <- toString(unmasked(alignment)[[i]])
    cat(the.sequence)
    cat('\n')  
  }
  
  sink(NULL)
}

alignment2Fasta(Interesting_bacterial_OTU_fasta_aligned, 'bacterial_alignment.fasta')
#Interesting_OTU_alignment_DNAbin <- msaConvert(Interesting_bacterial_OTU_fasta_aligned, "ape::DNAbin")
sink() 

tree <- read.tree(file = "bacterial_tree.tree")
plot(tree)

t_rel_data_reduced <- as.data.frame(t(rel_data_reduced))


phylocom$traits
phylocom$sample

phy <- phylocom$phylo
phy <- tree
class(phy)
phy$tip.label

phylocom
data(phylocom)
names(phylocom)

class(phylocom$phylo)
plot(phylocom$phylo)


pd_result <- pd(rel_data_reduced, tree, include.root = FALSE)
saveRDS(pd_result, file = "pd_result.RDS")
aggregate(pd_result$PD, by = list(metadata$Stratum), FUN = mean)

phydist <- cophenetic(tree)
comdist.result <- comdist(rel_data_reduced, phydist)
saveRDS(comdist.result, file = "comdist.result.RDS")

comdistnt.data <- comdistnt(rel_data_reduced, phydist, 
          abundance.weighted = FALSE, exclude.conspecifics = FALSE)


#ses.mpd.result <- 
# ses.mpd(rel_data_reduced, phydist, null.model="taxa.labels", abundance.weighted=FALSE, runs=99)
#saveRDS(ses.mpd.result, file = "ses.mpd.result.RDS")
ses.mpd.result <- readRDS(file = "ses.mpd.result.RDS")

aggregate(ses.mpd.result$mpd.obs.z, by = list(metadata$Stratum), FUN = mean)
#mpd.obs.z is equivalent to -NRI

ses.mntd.result <- 
  ses.mntd(rel_data_reduced, phydist, null.model="taxa.labels", abundance.weighted=FALSE, runs=99)
saveRDS(ses.mntd.result, file = "ses.mntd.result.RDS")

aggregate(ses.mntd.result$mntd.obs.z, by = list(metadata$Stratum), FUN = mean)
aggregate(ses.mntd.result$mntd.obs.p, by = list(metadata$Stratum), FUN = mean)

mntd_results <- mntd(rel_data_reduced, phydist, abundance.weighted=FALSE)    
aggregate(mntd_results, by = list(metadata$Stratum), FUN = mean)
saveRDS(mntd_results, file = "mntd_results.RDS")




ses.mpd.result.weighted <- 
  ses.mpd(rel_data_reduced, phydist, null.model="taxa.labels", abundance.weighted = TRUE, runs=99)
saveRDS(ses.mpd.result.weighted, file = "ses.mpd.result.weighted")

aggregate(ses.mpd.result.weighted$mpd.obs.z, by = list(metadata$Stratum), FUN = mean)
#mpd.obs.z is equivalent to -NRI

ses.mntd.result.weighted <- 
  ses.mntd(rel_data_reduced, phydist, null.model="taxa.labels", abundance.weighted = TRUE, runs=99)
saveRDS(ses.mntd.result.weighted, file = "ses.mntd.result.weighted.RDS")

aggregate(ses.mntd.result.weighted$mntd.obs.z, by = list(metadata$Stratum), FUN = mean)
aggregate(ses.mntd.result.weighted$mntd.obs.p, by = list(metadata$Stratum), FUN = mean)

mntd_results.weighted <- mntd(rel_data_reduced, phydist, abundance.weighted = TRUE)    
aggregate(mntd_results.weighted, by = list(metadata$Stratum), FUN = mean)
saveRDS(mntd_results.weighted, file = "mntd_results.weighted.RDS")

mntd_results.weighted <- readRDS(file = "mntd_results.weighted.RDS")
