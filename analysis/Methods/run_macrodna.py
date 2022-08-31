import pandas as pd
import sys
import os

sys.path.append("../../src/MaCroDNA/")
from MaCroDNA import MaCroDNA

if __name__ == "__main__":
    # gene_preprocess should be in
    # ["2000_genes_log", "all_genes_log", "all-genes_raw", "noX_genes_raw"]

    print(os.getcwd())
    gene_preprocess = sys.argv[1]
    assert gene_preprocess in ["2000_genes_log", "all_genes_log", "all-genes_raw", "noX_genes_raw"]
    patient_list = ["crc04", "crc10", "crc11"]
    cluster_method = ["intNMF", "Agglomerative"]
    norm_method = ["log", "untransformed"]
    intnmf_cluster_num = {("crc04", "untransformed"): 2,
                          ("crc10", "untransformed"): 2,
                          ("crc11", "untransformed"): 2,
                          ("crc04", "log"): 2,
                          ("crc10", "log"): 3,
                          ("crc11", "log"): 3}

    dataDir = "../data/" + gene_preprocess + "/"
    clusterDir = "../data/clusters/"
    resultDir = "Result/"
    if not os.path.exists(resultDir):
        os.mkdir(resultDir)

    for tmp_patient in patient_list:
        for tmp_cluster_method in cluster_method:
            for tmp_norm_method in norm_method:

                rna_path = dataDir + tmp_patient + "_rna.csv"
                dna_path = dataDir + tmp_patient + "_dna.csv"

                rna_data = pd.read_csv(rna_path, index_col=0)
                dna_data = pd.read_csv(dna_path, index_col=0)

                if tmp_cluster_method == "Agglomerative":
                    cluster_num = 4
                if tmp_cluster_method == "intNMF":
                    cluster_num = intnmf_cluster_num[(tmp_patient, tmp_norm_method)]

                tmp_clusterDir = clusterDir + tmp_cluster_method + "/" + tmp_norm_method + "/"
                dna_label_path = tmp_clusterDir + tmp_norm_method + "_" + tmp_patient + \
                                 "_" + str(cluster_num) + ".csv"
                dna_label = pd.read_csv(dna_label_path, index_col=0)

                run_model = MaCroDNA(rna_data, dna_data, dna_label)
                cell2clone = run_model.cell2clone_assignment()
                cell2clone["cell"] = cell2clone.index
                cell2clone.index.name = None
                cell2clone = cell2clone.drop(columns=["predict_cell"])
                cell2clone = cell2clone.rename(columns={"predict_clone": "predict"})
                result_filename = "_".join(["MaCroDNA", tmp_cluster_method, tmp_norm_method,
                                            tmp_patient, str(cluster_num), "result.csv"])
                cell2clone.to_csv(resultDir + result_filename)
