rm(list=ls())
setwd("/home/xiru/Dropbox/Research/comparison/Preprocess/")
df_cnv <- read.csv("../Data/all_cnv_no_X.csv",sep = ',',row.names=1)
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
cnv_gr <- makeGRangesFromDataFrame(df_cnv, keep.extra.columns = TRUE)
olaps <- findOverlaps(g, cnv_gr)

# add chr start end width 2022.1.27
df_gene <- data_frame(chr = as.vector(seqnames(cnv_gr)[subjectHits(olaps)]),
                      start = start(ranges(cnv_gr))[subjectHits(olaps)],
                      end = end(ranges(cnv_gr))[subjectHits(olaps)],
                      width = width(ranges(cnv_gr))[subjectHits(olaps)],
                      entrezgene = names(g)[queryHits(olaps)], 
                      copy_number = mcols(cnv_gr)$copy_number[subjectHits(olaps)],
                      clone = mcols(cnv_gr)$clone[subjectHits(olaps)])


entrezgene_ensembl_map <- as.list(org.Hs.egENSEMBL)
entrezgene_ensembl_map <- lapply(entrezgene_ensembl_map, `[`, 1)

df_gene <- dplyr::filter(df_gene, entrezgene %in% names(entrezgene_ensembl_map)) %>% 
  dplyr::mutate(ensembl_gene_id = unlist(entrezgene_ensembl_map[entrezgene])) %>% 
  dplyr::select(chr, start, end, width, ensembl_gene_id, entrezgene, copy_number, clone) %>% 
  drop_na()



df_gene <- dplyr::count(df_gene,ensembl_gene_id) %>%
  dplyr::filter(n == length(unique(df_gene$clone))) %>%
  inner_join(df_gene) %>%
  dplyr::select(-n)

df_gene_expanded <- spread(df_gene, clone, copy_number)
cnv_mat_region <- dplyr::select(df_gene_expanded, -ensembl_gene_id, -entrezgene) %>% 
  as.matrix()

cnv_mat <- dplyr::select(df_gene_expanded, -ensembl_gene_id, -entrezgene, -chr,-start,-end,-width) %>% 
  as.matrix()

rownames(cnv_mat) <- df_gene_expanded$ensembl_gene_id
rownames(cnv_mat_region) <- df_gene_expanded$ensembl_gene_id
# keep_gene <- rowMins(cnv_mat) <= 6 & rowVars(cnv_mat) > 0
# 
# cnv_mat <- cnv_mat[keep_gene,]

# write.csv(cnv_mat_region,"../Data/cnv_processed_region.csv")
write.csv(cnv_mat,"../Data/all_cnv_no_X_processed.csv")

cnv_save <- as.data.frame(cnv_mat)
