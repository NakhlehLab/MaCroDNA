rm(list=ls())
# setwd("/home/xiru/Dropbox/Research/Datasets/scTrio_seq2/CNV_Bian/")
df_rna <- read.csv("../Data/all_rna.csv",sep = ',',row.names=1)

library("AnnotationDbi")
library("org.Hs.eg.db")

gene_name = rownames(df_rna)

df_rna$name=rownames(df_rna)

df_rna$gene_id <- mapIds(org.Hs.eg.db,
                           keys=gene_name, 
                           column="ENSEMBL",
                           keytype="SYMBOL",
                           fuzzy = TRUE,
                           multiVals="first")
write.csv(df_rna,"../Data/all_gene_id_name.csv",row.names = TRUE,sep="," )
