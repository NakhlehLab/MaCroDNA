library(IntNMF)
setwd("/home/xiru/Downloads/scripts/scripts/CRC_data_analysis/")



norm_method = c("untransformed", "log")

dataDir <- "Cluster/cluster_data/"
clusterDir <- "Cluster/"


for (tmp_norm_method in norm_method){
  
  resultDir = paste0(clusterDir,"intNMF_cluster_",tmp_norm_method,"/")
  dir.create(file.path(resultDir), showWarnings = FALSE)
  dna_path = paste0(dataDir,"dna_",tmp_norm_method,".csv")
  dna_data = read.csv(dna_path,row.names = 1)
  dna_matrix = data.matrix(dna_data)
  dna_matrix = t(dna_matrix)
  zero_var_columns <- which(apply(dna_matrix, 2, var) == 0)
  dna_matrix <- dna_matrix[, -zero_var_columns]
  opt_k = nmf.opt.k(dat = dna_matrix)
  opt_k_1 = as.integer(substring(names(which.max(rowMeans(opt_k))),2))
  fit <- nmf.mnnals(dat=dna_matrix, k=opt_k_1, seed = FALSE)
  cluster_id = fit$clusters
  cell_ls = colnames(dna_data)
  
  cluster_result = list(clone=cluster_id-min(cluster_id), cell=cell_ls)
  cluster_result = as.data.frame(do.call(cbind, cluster_result))
  rownames(cluster_result) = NULL
  filename = paste0(opt_k_1,"_clusters.csv")
  write.csv(cluster_result,paste0(resultDir,filename), row.names = seq(0,length(cell_ls)-1))
}
