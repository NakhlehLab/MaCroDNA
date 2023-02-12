library(phylosignal)
library(adephylo)
library(ape)
library(phylobase)
library(phangorn)

setwd("/Users/edrisi/Documents/MaCroDNA/")
biopsy.name <- "PAT14_HGD"
dir_ <- "./data/Busslinger_data/macrodna_paired/"

dat.dna <- read.csv(paste(dir_, biopsy.name, "_dna.csv", sep = ""), header = TRUE)
# sets the names of cells as the index of the data frame
rownames(dat.dna) <- dat.dna$X
# print the number of rows in the data frame
print(nrow(dat.dna))
# print the number of columns in the data frame
print(ncol(dat.dna))
# drop the first column 
dat.dna <- dat.dna[-c(1)]
dat.dna <- as.data.frame(dat.dna)
dat.dna.t <- t(dat.dna)
# print the number of rows in the data frame
print(nrow(dat.dna.t))
# print the number of columns in the data frame
print(ncol(dat.dna.t))

# compute the Euclidean pairwise distance matrix of the scDNA-seq cells
dm <- dist(dat.dna.t)
# infer the UPGMA tree 
tre <- upgma(dm)

# read the rna data with rpm normalization
dat <- read.csv(paste(dir_, biopsy.name, "_rna.csv", sep = ""), header = TRUE)
# sets the names of cells as the index of the data frame
rownames(dat) <- dat$X
# drop the first column 
dat <- dat[-c(1)]
dat <- as.data.frame(dat)

# create the phylo4d object using the count data and the tree 
p4d <- phylo4d(tre, dat)
any(is.na(dat))

# compute the correlation values for all the genes using all the methods
phylo.signal <- phyloSignal(p4d = p4d, method = "K.star")
# save the correlation values and their p-values into CSV files
write.csv(phylo.signal$stat, paste(biopsy.name, "_stat", ".csv", sep = ""))
write.csv(phylo.signal$pvalue, paste(biopsy.name, "_pval", ".csv", sep = ""))
# save the phylogenetic tree of the current biopsy into a newick file
ape::write.tree(tre, file=paste(biopsy.name, "_tre", ".nwk", sep = ""))