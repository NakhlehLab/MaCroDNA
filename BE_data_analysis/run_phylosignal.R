library(phylosignal)
library(adephylo)
library(ape)
library(phylobase)
library(phangorn)

setwd("/Users/edrisi/Documents/ongoing_projects/MaCroDNA/Dec_2022/")
biopsy.name <- "PAT20_CARD"

dat.dna <- read.csv(paste("./data/Busslinger_data/macrodna_paired/",biopsy.name, "_dna.csv", sep = ""), header = TRUE)
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

dm <- dist(dat.dna.t)
#tre <- NJ(dm)
#tre$edge.length[which(tre$edge.length <= 0)] <- 0.0000001
tre <- upgma(dm)
#plot(tre, main = "NJ")
#plot(tre, main = "UPGMA")

# read the rna data with rpm normalization
dat <- read.csv(paste("./data/Busslinger_data//macrodna_paired/",biopsy.name, "_rna.csv", sep = ""), header = TRUE)
# sets the names of cells as the index of the data frame
rownames(dat) <- dat$X
# print the number of rows in the data frame
print(nrow(dat))
# print the number of columns in the data frame
print(ncol(dat))
# drop the first column 
dat <- dat[-c(1)]
dat <- as.data.frame(dat)
# drop the columns whose number of zeros is above 95 percent
# dat <- dat[, which((colSums(dat==0)/nrow(dat))<0.95)]
# dat <- dat[, which(colSums(dat==0)>1)]
# print the number of rows in the data frame
print(nrow(dat))
# print the number of columns in the data frame
print(ncol(dat))
# put all the cosmic gene names into a list
# cosmic.genes <- colnames(dat)
# create the phylo4d object
# tre <- read.tree(paste("./trees/",biopsy.name,"_upgma.newick", sep = ""))
# z score normalization for the gene expression data 
# dat <- scale(dat, center = TRUE, scale = TRUE)
p4d <- phylo4d(tre, dat)
any(is.na(dat))
 
# # create a new directory to store the Moran's indices for each gene and sample there
# dir.create(paste(biopsy.name, "_Moran_plot_info", sep = ""))
# # iterate over all cosmic genes and compute Moran's index for each
# cnt = 1
# for (gene in cosmic.genes){
#   crlg <- phyloCorrelogram(p4d, trait = gene, ci.bs = 2, n.points = 50)
#   write.csv(crlg$res, file = paste(biopsy.name, "_Moran_plot_info/", gene, ".csv", sep = ""))
#   print(paste("processing is done for gene ", gene, sep = ":"))
#   print(cnt)
#   cnt = cnt + 1
# }

# compute the correlation values for all the genes using all the methods
phylo.signal <- phyloSignal(p4d = p4d, method = "K.star")
write.csv(phylo.signal$stat, paste(biopsy.name, "_stat", ".csv", sep = ""))
write.csv(phylo.signal$pvalue, paste(biopsy.name, "_pval", ".csv", sep = ""))
ape::write.tree(tre, file=paste(biopsy.name, "_tre", ".nwk", sep = ""))














######################################
# add pseudocounts to the data frame
#dat <- dat + 1

# phylo.signal <- phyloSignal(p4d = p4d, method = "all")
# print(phylo.signal)
# name_ <- "PAT9_NDBE"
# write.csv(phylo.signal$stat, paste(name_, "_stat.csv"))
# write.csv(phylo.signal$pvalue, paste(name_, "_pvalue.csv"))
#
# crlg <- phyloCorrelogram(p4d, trait = c(2))
# # pdf('PAT9_NDBE_nj_l1.pdf')
# # print(crlg)
# # typeof(crlg)
# # typeof(crlg$res)
# plot(crlg)
# # dev.off()
# #
# # #capture.output(crlg$res, file = "crlg_res.txt")
# # write.csv(crlg$res, file = "crlg_res.csv")
#
# for (i in colnames(dat)){
#   print(i)
# }
