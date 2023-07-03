rm(list=ls())
set.seed(1)
setwd("/home/xiru/Dropbox/Research/2022/MaCroDNA/Seurat/")
library(Seurat)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(dplyr)
library(caret)


patient_list = c("crc04","crc10","crc11")

cluster_method = c("intNMF","agglomerative")
norm_method = c("log","raw")
intnmf_cluster_num <- data.frame(patient = c("crc04","crc10","crc11","crc04","crc10","crc11"),
                                 norm = c("raw","raw","raw","log","log","log"),
                                 cluster_num=c(2,2,2,2,3,3))


  
dataDir <- "../Data/CRC/Seurat_preprocess/2000_genes/"
clusterDir <- "../Data/CRC/same_gene_with_X_no_0/Clusters/"
resultDir <- "Result/"
rna_name <- "_rna_"
dna_name <- "_dna_"

dir.create(file.path(resultDir), showWarnings = FALSE)

for (tmp_patient in patient_list) {
  for (tmp_cluster_method in cluster_method) {
    for (tmp_norm_method in norm_method) {
      
      rna_path = paste0(dataDir,tmp_patient,"_rna.RDS")
      dna_path = paste0(dataDir,tmp_patient,"_dna.RDS")
      gene.2000_path = paste0(dataDir,tmp_patient,"_2000genes.RDS")
      crc.rna = readRDS(rna_path)
      crc.cnv = readRDS(dna_path)
      gene.2000 = readRDS(gene.2000_path)
      
      crc.anchors <- FindTransferAnchors(reference = crc.cnv, query = crc.rna,features=gene.2000)
      
      tmp_clusterDir = paste0(clusterDir,tmp_cluster_method,"_",tmp_norm_method,"/")
      
      if (tmp_cluster_method == "agglomerative"){
        cluster_num = 4
      } 
      
      if(tmp_cluster_method == "intNMF"){
        cluster_num = intnmf_cluster_num[((intnmf_cluster_num$patient==tmp_patient)&
                                            (intnmf_cluster_num$norm==tmp_norm_method)), ]$cluster_num
      }
      
      dna_label_path = paste0(tmp_clusterDir,"dna_cluster_",tmp_patient,"_",cluster_num,".csv")
      
      
      c = cluster_num
      ref_dna = read.csv(dna_label_path, header=TRUE,row.names = 1)
      # ref_dna <- ref_dna[match(colnames(rna_data), ref_dna$cell), ] # ysb
      rownames(ref_dna) <- ref_dna$cell
      ref_dna = ref_dna[Cells(crc.cnv),]
      # result_dna = ref_dna
      name <- ref_dna$cell
      ref_dna <- as.vector(ref_dna$clone)
      names(ref_dna) <- name
      ref_dna <- as.factor(ref_dna)
      crc.prediction <- TransferData(anchorset = crc.anchors, refdata = ref_dna, 
                                     k.weight=5)
      
      clones = crc.prediction$predicted.id
      result_dna <- data.frame("predict"=clones)
      result_dna$cell <- Cells(crc.rna)
      result_path = paste0(resultDir,tmp_patient, "_",tmp_cluster_method,"_",
                           tmp_norm_method,"_cluster_",c,"_result.csv")
      write.csv(result_dna,result_path,row.names = TRUE, col.names = TRUE)
      }
    }
  }
