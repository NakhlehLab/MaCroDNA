rm(list=ls())
library(Seurat)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(dplyr)
library(caret)
setwd("/home/xiru/Dropbox/Research/comparison/Preprocess/")


patient_list = c("crc04","crc10","crc11")


tmp_patient = patient_list[[1]]

  
mainDir <- "../Data/same_gene_with_X_no_0/"
resultDir <- "../Data/Seurat_preprocess/"
rna_name <- "_rna"
dna_name <- "_dna"

dir.create(file.path(resultDir), showWarnings = FALSE)

for (tmp_patient in patient_list) {
  rna_path = paste0(mainDir,tmp_patient,"_rna_raw.csv")
  dna_path = paste0(mainDir,tmp_patient,"_dna_raw.csv")
  rna_data = read.csv(rna_path,row.names = 1)
  cnv_data = read.csv(dna_path,row.names = 1)
  cnv_data = data.matrix(cnv_data)
  rna_data = data.matrix(rna_data)
  
  crc.rna = CreateSeuratObject(counts=rna_data, assay="RNA", 
                               names.field=3, names.delim="_")
  crc.cnv = CreateSeuratObject(counts=cnv_data, assay="DNA", 
                               names.field=3, names.delim="_")
  
  crc.rna <- NormalizeData(crc.rna)
  crc.rna <- ScaleData(crc.rna)
  crc.rna <- FindVariableFeatures(crc.rna, nfeatures = 4000)
  
  
  crc.cnv <- NormalizeData(crc.cnv)
  crc.cnv <- ScaleData(crc.cnv)
  crc.cnv <- FindVariableFeatures(crc.cnv, nfeatures = 4000)
  
  
  
  
  IntgFeature <- SelectIntegrationFeatures(c(crc.cnv, crc.rna))
  
  
  rna_pre_all = as.data.frame(GetAssayData(crc.rna))
  dna_pre_all =as.data.frame(GetAssayData(crc.cnv))
  
  rna_pre_2000 = rna_pre_all[IntgFeature,]
  dna_pre_2000 = dna_pre_all[IntgFeature,]
  
  resultDir <- "../Data/Seurat_preprocess/"
  all.resultDir = paste0(resultDir,"all_genes/")
  dir.create(file.path(all.resultDir),showWarnings = FALSE)
  write.csv(rna_pre_all, paste0(all.resultDir, tmp_patient,"_rna_raw.csv"))
  write.csv(dna_pre_all, paste0(all.resultDir, tmp_patient,"_dna_raw.csv"))
  
  two.resultDir = paste0(resultDir,"2000_genes/")
  dir.create(file.path(two.resultDir),showWarnings = FALSE)
  write.csv(rna_pre_2000, paste0(two.resultDir, tmp_patient,"_rna_raw.csv"))
  write.csv(dna_pre_2000, paste0(two.resultDir, tmp_patient,"_dna_raw.csv"))
  saveRDS(crc.cnv, file=paste0(two.resultDir,tmp_patient, "_dna.RDS"))
  saveRDS(crc.rna, file=paste0(two.resultDir,tmp_patient,  "_rna.RDS"))
  saveRDS(IntgFeature, file=paste0(two.resultDir,tmp_patient,  "_2000genes.RDS"))
}
  
  