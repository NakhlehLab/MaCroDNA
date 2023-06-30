rm(list=ls())
setwd("/home/xiru/Dropbox/Research/2022/MaCroDNA/Preprocess/")
df_cnv <- read.csv("../Data/CRC/all_cnv.csv",sep = ',',row.names=1)
# df_cnv = df_cnv[,-1]
df_cnv = tibble::tibble(df_cnv)

#load("dat.RData")
#df_cnv = cnv_data[,-5]
#df_cnv = tibble::tibble(df_cnv)
#names = unique(df_cnv$clone)
#df_cnv$clone[df_cnv$clone == names[1]] = "A"
#df_cnv$clone[df_cnv$clone == names[2]] = "B"
# library(clonealign)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(dplyr)
library(org.Hs.eg.db)
library(tidyr)
library(matrixStats)
# data(df_cnv)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
g <- genes(txdb, single.strand.genes.only=FALSE)

#df_cnv <- dplyr::mutate(df_cnv, chr = paste0("chr", chr))

gene_list = read.csv("gene_ordering_file.txt", 
                     sep = "\t",
                     col.names=c("names","chr","start","end"))
gene_gr <- makeGRangesFromDataFrame(gene_list, keep.extra.columns = TRUE)
gene_gr = setNames(gene_gr,gene_list$names)

cnv_gr <- makeGRangesFromDataFrame(df_cnv, keep.extra.columns = TRUE)

gene_olaps<-findOverlaps(gene_gr, cnv_gr)

v2_df_gene <- data_frame(chr = as.vector(seqnames(cnv_gr)[subjectHits(gene_olaps)]),
                         start = start(ranges(cnv_gr))[subjectHits(gene_olaps)],
                         end = end(ranges(cnv_gr))[subjectHits(gene_olaps)],
                         width = width(ranges(cnv_gr))[subjectHits(gene_olaps)],
                         gene = names(gene_gr)[queryHits(gene_olaps)], 
                         copy_number = mcols(cnv_gr)$copy_number[subjectHits(gene_olaps)],
                         clone = mcols(cnv_gr)$clone[subjectHits(gene_olaps)])

v2_df_gene <- dplyr::count(v2_df_gene,gene) %>%
  dplyr::filter(n == length(unique(v2_df_gene$clone))) %>%
  inner_join(v2_df_gene) %>%
  dplyr::select(-n)

V2_df_gene_expanded <- spread(v2_df_gene, clone, copy_number)

V2_cnv_mat <- dplyr::select(V2_df_gene_expanded, -chr,-start,-end,-width, -gene) %>% 
  as.matrix()
rownames(V2_cnv_mat) <- V2_df_gene_expanded$gene

write.csv(V2_cnv_mat,"../Data/CRC/all_cnv_processed_gene_name.csv")
