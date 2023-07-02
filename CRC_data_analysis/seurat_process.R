rm(list=ls())
library(Seurat)
setwd("/home/xiru/Downloads/scripts/scripts/CRC_data_analysis/")
set.seed(333)

# get the # of genes need to be kept, either "all" or "2000"
myargs = commandArgs(trailingOnly=TRUE)
gene_amount = myargs[1]

# read dna and rna data
DataDir = "Data/"
dna_data = read.csv(paste0(DataDir,"dna_all_genes_raw.csv"),row.names=1)
rna_data = read.csv(paste0(DataDir,"rna_all_genes_raw.csv"),row.names=1)

# create Seurat object
rna_seurat = CreateSeuratObject(counts=rna_data, assay="RNA", 
                             names.field=1, names.delim="_")
dna_seurat = CreateSeuratObject(counts=dna_data, assay="DNA", 
                             names.field=1, names.delim="_")

# Seurat processing
rna_seurat <- NormalizeData(rna_seurat)
rna_seurat <- ScaleData(rna_seurat)
rna_seurat <- FindVariableFeatures(rna_seurat)

dna_seurat <- NormalizeData(dna_seurat)
dna_seurat <- ScaleData(dna_seurat)
dna_seurat <- FindVariableFeatures(dna_seurat)


rna_res_all = as.data.frame(GetAssayData(rna_seurat))
dna_res_all =as.data.frame(GetAssayData(dna_seurat))

if (gene_amount != "all"){
  gene_amount = as.integer(gene_amount)
  IntgFeature <- SelectIntegrationFeatures(c(rna_seurat, dna_seurat),
                                           nfeatures=gene_amount)
  rna_res_part = rna_res_all[IntgFeature,]
  dna_res_part = dna_res_all[IntgFeature,]
  
}

# save data
if(gene_amount=="all"){
  write.csv(rna_res_all, paste0(DataDir,"rna_all_genes_log.csv"))
  write.csv(dna_res_all, paste0(DataDir,"dna_all_genes_log.csv"))
} else {
  write.csv(rna_res_part, paste0(DataDir,"rna_",gene_amount,"_genes_log.csv"))
  write.csv(dna_res_part, paste0(DataDir,"dna_",gene_amount,"_genes_log.csv"))
} 
  
