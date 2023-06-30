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


norm_method = "log"

cluster_method = "agglomerative"
cluster_num = 4


get_data <- function(p,c_method,c_norm,c_num){
  clusterDir <- "../Data/CRC/same_gene_with_X_no_0/Clusters/"
  rna_name <- "_rna_"
  dna_name <- "_dna_"
  
  rm0_dataDir <- "../Data/CRC/clonealign_preprocess/rm0/"
  dataDir = rm0_dataDir
  
  rna_path = paste0(dataDir,p,rna_name,"raw.csv")
  rna = read.csv(rna_path,row.names=1)
  rna_matrix = data.matrix(rna)
  rna_sce = SingleCellExperiment(assays=list(counts=rna_matrix))
  
  dna_path = paste0(dataDir,p,dna_name,"raw.csv")
  dna_data = read.csv(dna_path, row.names = 1)
  
  tmp_clusterDir = paste0(clusterDir,c_method,"_",c_norm,"/")
  dna_label_path = paste0(tmp_clusterDir,"dna_cluster_",p,"_",c_num,".csv")
  
  dna_label = read.csv(dna_label_path,header=TRUE, sep=",",row.names=1)
  dna_cnv <- data.frame(matrix(ncol = 0, nrow=dim(dna_data)[1]))
  rownames(dna_cnv) <- rownames(dna_data)
  
  if (min(dna_label$clone)!=0){
    dna_label$clone = dna_label$clone-min(dna_label$clone)
  } 
  
  for (c_i in seq(0,c_num-1)) {
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
  
  return(list(rna_sce=rna_sce,dna=dna_cnv))
}

analyze_clonealign <- function(dna,rna_sce){
  ca_data <- preprocess_for_clonealign(rna_sce, dna)
  cal<- run_clonealign(ca_data$gene_expression_data, 
                       ca_data$copy_number_data)
  qplot(seq_along(cal$convergence_info$elbo), 
        cal$convergence_info$elbo, geom = c("point", "line")) +
    labs(x = "Iteration", y = "ELBO")
  # get assignment prob
  clone_prob = cal$ml_params$clone_probs
  clone_prob = colnames(dna)[max.col(clone_prob)]
  clone_prob = factor(clone_prob, levels=colnames(dna))
  clone_table=table(clone_prob)
  
  # highest proportions
  high_pp = apply(dna, 2, function(x) max(table(x))/sum(table(x)))
  high_item = apply(dna, 2, function(x) names(table(x))[which.max(table(x))])
  
  # combine info
  high_pp_df = t(data.frame(high_pp))
  res_table <- rbind(high_pp_df,high_item,
                     clone_table,apply(cal$ml_params$clone_probs, 2, mean))

  rownames(res_table) <- c("highest_pp", "highest_item","clonealign_result","result_probability")
  
  return(list(res_table=res_table,clone_prob=cal$ml_params$clone_probs))
}

cnv_rm_exit_0 <- function(cnv_data){
  col_var <- apply(cnv_data,2,var)
  cnv_data <- cnv_data[,col_var != 0]
  return(cnv_data)
}

cnv_add_var0 <- function(cnv_data,cnv_num){
  # add new column
  cnv_data[[paste0("var0_cnv",cnv_num)]] <- cnv_num
  return(cnv_data)
}

# get the data
crc04_data = get_data("crc04",cluster_method,norm_method,cluster_num)
crc10_data = get_data("crc10",cluster_method,norm_method,cluster_num)
crc11_data = get_data("crc11",cluster_method,norm_method,cluster_num)

############################
# original
############################
exp0 <- function(p_data,p_name){
  new_dna = p_data$dna
  res = analyze_clonealign(new_dna,p_data$rna_sce)
  print("########################")
  print(paste0("Result for patient: ", p_name))
  print(res$res_table)
  print("########################")
  return(res)
}

crc04_exp0 = exp0(crc04_data,"crc04")
crc10_exp0 = exp0(crc10_data,"crc10")
crc11_exp0 = exp0(crc11_data,"crc11")


############################
# original: rm 
############################
exp0_rm <- function(p_data,p_name){
  new_dna = cnv_rm_exit_0(p_data$dna)
  res = analyze_clonealign(new_dna,p_data$rna_sce)
  print("########################")
  print(paste0("Result for patient: ", p_name))
  print(res$res_table)
  print("########################")
  return(res)
}

crc04_exp0_rm = exp0_rm(crc04_data,"crc04")
crc10_exp0_rm = exp0_rm(crc10_data,"crc10")
crc11_exp0_rm = exp0_rm(crc11_data,"crc11")



############################
# exp1: add clone with all 1
############################
exp_add_cnv1 <- function(p_data,p_name){
  new_dna = cnv_rm_exit_0(p_data$dna)
  new_dna = cnv_add_var0(new_dna,1)
  res = analyze_clonealign(new_dna,p_data$rna_sce)
  print("########################")
  print(paste0("Result for patient: ", p_name))
  print(res$res_table)
  print("########################")
  return(res)
}

crc04_add_cnv1 = exp_add_cnv1(crc04_data,"crc04")
crc10_add_cnv1 = exp_add_cnv1(crc10_data,"crc10")
crc11_add_cnv1 = exp_add_cnv1(crc11_data,"crc11")


############################
# exp1: add clone with all 5
############################
exp_add_cnv5 <- function(p_data,p_name){
  new_dna = cnv_rm_exit_0(p_data$dna)
  new_dna = cnv_add_var0(new_dna,5)
  res = analyze_clonealign(new_dna,p_data$rna_sce)
  print("########################")
  print(paste0("Result for patient: ", p_name))
  print(res$res_table)
  print("########################")
  return(res)
}

crc04_add_cnv5 = exp_add_cnv5(crc04_data,"crc04")
crc10_add_cnv5 = exp_add_cnv5(crc10_data,"crc10")
crc11_add_cnv5 = exp_add_cnv5(crc11_data,"crc11")


############################
# exp2: change proportion
############################
exp_high_pp <- function(p_data,p_name,portion){
  
  new_dna = cnv_rm_exit_0(p_data$dna)
  
  # calculate the number of values that should be equal to 2
  n <- nrow(new_dna)
  n_2 <- ceiling(n*portion)
  
  
  # generate new column
  values <- c(rep(2,n_2), 
              sample(c(1,3,4,5), n-n_2, replace = TRUE))
  values <- sample(values)
  new_dna[[paste0("pp_",portion)]] = values
  
  res = analyze_clonealign(new_dna,p_data$rna_sce)
  print("########################")
  print(paste0("Result for patient: ", p_name))
  print(res$res_table)
  print("########################")
  return(res)
}

crc04_high_pp_1 = exp_high_pp(crc04_data,"crc04",0.9)
crc10_high_pp_1 = exp_high_pp(crc10_data,"crc10",0.9)
crc11_high_pp_1 = exp_high_pp(crc11_data,"crc11",0.9)

crc04_high_pp_2 = exp_high_pp(crc04_data,"crc04",0.8)
crc10_high_pp_2 = exp_high_pp(crc10_data,"crc10",0.8)
crc11_high_pp_2 = exp_high_pp(crc11_data,"crc11",0.8)

crc04_high_pp_3 = exp_high_pp(crc04_data,"crc04",0.7)
crc10_high_pp_3 = exp_high_pp(crc10_data,"crc10",0.7)
crc11_high_pp_3 = exp_high_pp(crc11_data,"crc11",0.7)

############################
# exp3: 2 uniform clones
############################

exp_2uniform <- function(p_data,p_name,cnv_num,s){
  new_dna = cnv_add_var0(p_data$dna,cnv_num)
  set.seed(s)
  res = analyze_clonealign(new_dna,p_data$rna_sce)
  print("########################")
  print(paste0("Result for patient: ", p_name))
  print(res$res_table)
  print("########################")
  return(res)
}

crc04_2uniform_2 = exp_2uniform(crc04_data,"crc04",5,222)
crc10_2uniform_2 = exp_2uniform(crc10_data,"crc10",5,222)
crc11_2uniform_2 = exp_2uniform(crc11_data,"crc11",5,222)

crc04_2uniform_1 = exp_2uniform(crc04_data,"crc04",5,333)
crc10_2uniform_1 = exp_2uniform(crc10_data,"crc10",5,333)
crc11_2uniform_1 = exp_2uniform(crc11_data,"crc11",5,333)

crc04_2uniform_3 = exp_2uniform(crc04_data,"crc04",5,5)
crc10_2uniform_3 = exp_2uniform(crc10_data,"crc10",5,5)
crc11_2uniform_3 = exp_2uniform(crc11_data,"crc11",5,5)
