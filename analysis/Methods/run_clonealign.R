rm(list=ls())
set.seed(1)


library(tensorflow)
# you should change the conda environment to the one where tensorflow is
use_condaenv(condaenv = "py3")
library(clonealign)
library(SingleCellExperiment)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(caret)
library(stringr)



patient_list = c("crc04","crc10","crc11")
cluster_method = c("intNMF","Agglomerative")
norm_method = c("log","untransformed")
intnmf_cluster_num <- data.frame(patient = c("crc04","crc10","crc11","crc04","crc10","crc11"),
                                 norm = c("untransformed","untransformed","untransformed","log","log","log"),
                                 cluster_num=c(2,2,2,2,3,3))

# tmp_patient = patient_list[[2]]
# tmp_cluster_method = cluster_method[[2]]
# tmp_norm_method = norm_method[[1]]


dataDir <- "../data/noX_genes_raw/"
clusterDir <- "../data/clusters/"
resultDir <- "Result/"
dir.create(file.path(resultDir), showWarnings = FALSE)

for (tmp_patient in patient_list) {
  for (tmp_cluster_method in cluster_method) {
    for (tmp_norm_method in norm_method) {
      
      rna_path = paste0(dataDir,tmp_patient,"_rna.csv")
      rna = read.csv(rna_path,row.names=1)
      rna_matrix = data.matrix(rna)
      rna_sce = SingleCellExperiment(assays=list(counts=rna_matrix))
      
      dna_path = paste0(dataDir,tmp_patient,"_dna.csv")
      dna_data = read.csv(dna_path, row.names = 1)
      
      # if don't use this, it will cause "Initial elbo is NA"
      # https://github.com/kieranrcampbell/clonealign/issues/4
      dna_data[dna_data == 0] <- 1
      
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
      

      
      dna_label = read.csv(dna_label_path,header=TRUE, sep=",",row.names=1)
      dna_cnv <- data.frame(matrix(ncol = 0, nrow=dim(dna_data)[1]))
      rownames(dna_cnv) <- rownames(dna_data)
      
      if (min(dna_label$clone)!=0){
        dna_label$clone = dna_label$clone-min(dna_label$clone)
      } 
      
      for (c_i in seq(0,cluster_num-1)) {
        c_i_cell = dna_label[dna_label$clone==c_i,]$cell
        c_i_cnv_all = dna_data[,c_i_cell]
        if (length(c_i_cell) == 1){
          dna_cnv[paste0("clone",c_i)] = c_i_cnv_all
        } else{
          c_i_cnv_mean = apply(c_i_cnv_all,1,median)
          # c_i_cnv_mean = rowMeans(c_i_cnv_all)
          c_i_cnv_mean = round(c_i_cnv_mean, digits = 0)
          dna_cnv[paste0("clone",c_i)] = c_i_cnv_mean
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
      result_path = paste0(resultDir,"clonealign_",tmp_cluster_method,"_", 
                           tmp_norm_method,"_",tmp_patient, "_",cluster_num,"_result.csv")
      write.csv(result_dna,result_path,row.names = TRUE, col.names = TRUE)
      
    }
  }
}
