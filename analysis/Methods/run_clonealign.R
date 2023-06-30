rm(list=ls())
set.seed(1)
setwd("/home/xiru/Dropbox/Research/2022/MaCroDNA/clonealign/")

library(tensorflow)
use_condaenv(condaenv = "py3")
library(clonealign)
library(SingleCellExperiment)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(caret)
library(stringr)



patient_list = c("crc04","crc10","crc11")
cluster_method = c("agglomerative", "intNMF")
# cluster_method = c("intNMF")
norm_method = c("log","raw")
intnmf_cluster_num <- data.frame(patient = c("crc04","crc10","crc11","crc04","crc10","crc11"),
                                 norm = c("raw","raw","raw","log","log","log"),
                                 cluster_num=c(2,2,2,2,3,3))


clusterDir <- "../Data/CRC/same_gene_with_X_no_0/Clusters/"
rna_name <- "_rna_"
dna_name <- "_dna_"

add1_dataDir <- "../Data/CRC/clonealign_preprocess/add1/"
rm0_dataDir <- "../Data/CRC/clonealign_preprocess/rm0/"

dir.create("add1/", showWarnings = FALSE)
dir.create("rm0/", showWarnings = FALSE)

add1_resultDir <- "add1/Result/"
rm0_resultDir <- "rm0/Result/"

dir.create(add1_resultDir, showWarnings = FALSE)
dir.create(rm0_resultDir, showWarnings = FALSE)

run_method <- function(dataDir,resultDir,patient_list,cluster_method,norm_method){
  
  for (tmp_patient in patient_list) {
    for (tmp_cluster_method in cluster_method) {
      for (tmp_norm_method in norm_method) {
        
   
        rna_path = paste0(dataDir,tmp_patient,rna_name,"raw.csv")
        rna = read.csv(rna_path,row.names=1)
        rna_matrix = data.matrix(rna)
        rna_sce = SingleCellExperiment(assays=list(counts=rna_matrix))
        
        dna_path = paste0(dataDir,tmp_patient,dna_name,"raw.csv")
        dna_data = read.csv(dna_path, row.names = 1)
        
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
        
        dna_label = read.csv(dna_label_path,header=TRUE, sep=",",row.names=1)
        dna_cnv <- data.frame(matrix(ncol = 0, nrow=dim(dna_data)[1]))
        rownames(dna_cnv) <- rownames(dna_data)
        
        if (min(dna_label$clone)!=0){
          dna_label$clone = dna_label$clone-min(dna_label$clone)
        } 
        
        for (c_i in seq(0,c-1)) {
          c_i_cell = dna_label[dna_label$clone==c_i,]$cell
          c_i_cnv_all = dna_data[,c_i_cell]
          if (length(c_i_cell) == 1){
            dna_cnv[paste0("clone",c_i)] = c_i_cnv_all
          } else{
            c_i_cnv_median = apply(c_i_cnv_all,1,median)
            c_i_cnv_median = round(c_i_cnv_medium, digits = 0)
            dna_cnv[paste0("clone",c_i)] = c_i_cnv_median
          }
        }
        
        
        colnames(dna_cnv)=seq(0,dim(dna_cnv)[2]-1)
        
        rownames(dna_label) <- dna_label$cell
        same_cell= intersect(dna_label$cell, colnames(rna_sce))
        
        dna_label = dna_label[same_cell,]
        
        dna_label = dna_label$clone
        
        ca_data <- preprocess_for_clonealign(rna_sce, dna_cnv)
        
        cal <- run_clonealign(ca_data$gene_expression_data, ca_data$copy_number_data)
        result_dna <- data.frame("predict"=cal$clone)
        rownames(result_dna) = colnames(rna)
        
        result_dna$cell = colnames(rna)
        result_path = paste0(resultDir,tmp_patient, "_",tmp_cluster_method,"_",tmp_norm_method,"_cluster_",c,"_result.csv")
        write.csv(result_dna,result_path,row.names = TRUE, col.names = TRUE)
        
      }
    }
  }
}

run_method(add1_dataDir,add1_resultDir,patient_list,cluster_method,norm_method)

run_method(rm0_dataDir,rm0_resultDir, patient_list,cluster_method,norm_method)


