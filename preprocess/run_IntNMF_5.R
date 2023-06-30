library(IntNMF)
setwd("/home/xiru/Dropbox/Research/comparison/Preprocess/")

# patient_list = c("crc04")
# patient_list = c("crc04","crc10","crc11")
# norm_method = c("log","raw")
patient_list = c("crc11")
norm_method = c("raw")
intnmf_cluster_num <- data.frame(patient = c("crc04","crc10","crc11","crc04","crc10","crc11"),
                                 norm = c("raw","raw","raw","log","log","log"),
                                 cluster_num=c(2,2,2,2,3,3))

# tmp_patient = patient_list[[1]]
# tmp_norm_method = norm_method[[1]]
# ncluster = cluster_num[[1]]
# NOT_SAME = TRUE

  

mainDir <- "../Data/same_gene_with_X_no_0/"
rna_name <- "_rna_"
dna_name <- "_dna_"

# dir.create(file.path(resultDir), showWarnings = FALSE)

for (tmp_norm_method in norm_method){
  
  resultDir = paste0(mainDir,"Clusters/intNMF_",tmp_norm_method,"/")
  dir.create(file.path(resultDir), showWarnings = FALSE)
  
  for (tmp_patient in patient_list) {
    dna_path = paste0(mainDir,tmp_patient,dna_name,tmp_norm_method,".csv")
    dna_data = read.csv(dna_path,row.names = 1)
    dna_matrix = data.matrix(dna_data)
    dna_matrix = t(dna_matrix)
    # opt_k = nmf.opt.k(dat = dna_matrix)
    # opt_k_1 = as.integer(substring(names(which.max(rowMeans(opt_k))),2))
    opt_k_1 = intnmf_cluster_num[((intnmf_cluster_num$patient==tmp_patient)&
                                       (intnmf_cluster_num$norm==tmp_norm_method)), ]$cluster_num
    
    set.seed(330)
    fit <- nmf.mnnals(dat=dna_matrix, k=opt_k_1, seed = FALSE)
    cluster_id = fit$clusters
    cell_ls = colnames(dna_data)
    
    cluster_result = list(clone=cluster_id-min(cluster_id), cell=cell_ls)
    cluster_result = as.data.frame(do.call(cbind, cluster_result))
    rownames(cluster_result) = NULL
    filename = paste0("dna_cluster_",tmp_patient,"_",opt_k_1,".csv")
    write.csv(cluster_result,paste0(resultDir,filename), row.names = seq(0,length(cell_ls)-1))
    
  }
}

  

