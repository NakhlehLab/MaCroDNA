rm(list=ls())
set.seed(1)
library(Seurat)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(dplyr)
library(caret)
# setwd("/home/xiru/Dropbox/Research/GitHub/MaCroDNA/analysis/")

patient_list = c("crc04","crc10","crc11")
cluster_method = c("intNMF","Agglomerative")
norm_method = c("log","untransformed")
intnmf_cluster_num <- data.frame(patient = c("crc04","crc10","crc11","crc04","crc10","crc11"),
                                 norm = c("untransformed","untransformed","untransformed","log","log","log"),
                                 cluster_num=c(2,2,2,2,3,3))
  
dataDir <- "../data/2000_genes_log/"
clusterDir <- "../data/clusters/"
resultDir <- "Result/"
dir.create(file.path(resultDir), showWarnings = FALSE)

for (tmp_patient in patient_list) {
  for (tmp_cluster_method in cluster_method) {
    for (tmp_norm_method in norm_method) {
      
      rna_path = paste0(dataDir,tmp_patient,"_rna.csv")
      dna_path = paste0(dataDir,tmp_patient,"_dna.csv")
      
      dna_df = read.csv(dna_path, row.names=1)
      rna_df = read.csv(rna_path, row.names=1)
      
      cnv_data = data.matrix(dna_df)
      rna_data = data.matrix(rna_df)
      
      crc.rna = CreateSeuratObject(counts=rna_data, assay="RNA", 
                                   names.field=3, names.delim="_")
      crc.cnv = CreateSeuratObject(counts=cnv_data, assay="DNA", 
                                   names.field=3, names.delim="_")
      

      crc.anchors <- FindTransferAnchors(reference = crc.cnv, query = crc.rna,
                                         features = rownames(dna_df), k.filter = NA)
      
      tmp_clusterDir = paste0(clusterDir,tmp_cluster_method,"/",tmp_norm_method,"/")
      
      if (tmp_cluster_method == "Agglomerative"){
        cluster_num = 4
      } 
      
      if(tmp_cluster_method == "intNMF"){
        cluster_num = intnmf_cluster_num[((intnmf_cluster_num$patient==tmp_patient)&
                                            (intnmf_cluster_num$norm==tmp_norm_method)), ]$cluster_num
      }
      
      dna_label_path = paste0(tmp_clusterDir,tmp_norm_method,"_",
                              tmp_patient,"_",cluster_num,".csv")
      ref_dna = read.csv(dna_label_path, header=TRUE,row.names = 1)
      rownames(ref_dna) <- ref_dna$cell
      ref_dna = ref_dna[Cells(crc.cnv),]
      name <- ref_dna$cell
      ref_dna <- as.vector(ref_dna$clone)
      names(ref_dna) <- name
      ref_dna <- as.factor(ref_dna)
      crc.prediction <- TransferData(anchorset = crc.anchors, refdata = ref_dna, 
                                     k.weight=15)
      
      clones = crc.prediction$predicted.id
      result_dna <- data.frame("predict"=clones)
      result_dna$cell <- Cells(crc.rna)
      result_path = paste0(resultDir,"Seurat_",tmp_cluster_method,"_", 
                           tmp_norm_method,"_",tmp_patient, "_",cluster_num,"_result.csv")
      write.csv(result_dna,result_path,row.names = TRUE)
      }
    }
  }
