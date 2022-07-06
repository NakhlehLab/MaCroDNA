import pandas as pd
import sys


def get_rna_cells(rna_folder, patient_list=("crc04", "crc10", "crc11")):
    """
    get the rnba cell ids for each patient
    @param rna_folder: string, the directory where the rna files are
    @param patient_list:
    @return: rna_cell: dict, key is the patient id,
    value is the list of cell ids for this patient
    """
    rna_cell = {}
    for tmp_patient in patient_list:
        file_name = tmp_patient + "_rna.csv"
        rna_data = pd.read_csv(rna_folder + file_name, index_col=0)
        tmp_rna_cell = list(rna_data.columns)
        rna_cell[tmp_patient] = tmp_rna_cell
    return rna_cell

def get_Agg_clusters(cluster_tran, cluster_info,
                     cluster_method="Agglomerative",
                     patient_list=("crc04", "crc10", "crc11")):
    """
    Get the agglomerative clusters
    @param cluster_tran: string, normalisation method for the data being clustered,
        either "untransformed" or "log"
    @param cluster_info: ict, extra info for each clustering method
    @return: dna_label, dict, key is the tuple (cluster_method, cluster_tran, patient_id, cluster_key)
    value is the agglomerative clustering results and the clone id list
    """
    cluster_key = cluster_info[cluster_method]["cluster_key"][cluster_tran]
    cluster_dir = cluster_info[cluster_method]["cluster_dir"] + cluster_tran + "/"

    dna_label = {}
    for c_k in cluster_key:
        dna_file_name = c_k + ".csv"
        dna_clone_list_file_name = c_k + "_label.csv"
        tmp_dna_label = pd.read_csv(cluster_dir + dna_file_name, index_col=0)
        tmp_clone_list = pd.read_csv(cluster_dir + dna_clone_list_file_name, index_col=0)
        tmp_clone_list = list(tmp_clone_list["cluster_label"])
        tmp_dna_label = tmp_dna_label.set_index("cell")
        dna_cells = list(tmp_dna_label.index)

        for tmp_patient in patient_list:
            # select cells by patient
            tmp_patient_cell = [x for x in dna_cells if x.split("_")[0] == tmp_patient.upper()]
            tmp_patient_dna_label = tmp_dna_label.loc[tmp_patient_cell, :]
            tmp_patient_clone_list = [x for x in tmp_clone_list if x.split("_")[0] == tmp_patient]
            tmp_key = (cluster_method, cluster_tran, tmp_patient, c_k)
            tmp_patient_dna_label["clone"] = tmp_patient_dna_label["agg_clone"]
            tmp_patient_dna_label = tmp_patient_dna_label.drop(columns=["agg_clone", "original_clone"])
            dna_label[tmp_key] = {"dna_cluster":tmp_patient_dna_label, "clone_list": tmp_patient_clone_list}
    return dna_label


def get_intNMF_clusters(cluster_tran, cluster_info,
                        cluster_method="intNMF",
                        patient_list=("crc04", "crc10", "crc11")):
    """
    Get the agglomerative clusters
    @param cluster_tran: string, normalisation method for the data being clustered,
        either "untransformed" or "log"
    @param cluster_info: dict, extra info for each clustering method
    @return: dna_label, dict, key is the tuple (cluster_method, cluster_tran, patient_id, "None")
    value is the intNMF clustering results and the clone id list
    """
    cluster_key = cluster_info[cluster_method]["cluster_key"][cluster_tran]
    cluster_dir = cluster_info[cluster_method]["cluster_dir"] + cluster_tran + "/"

    dna_label = {}
    for c_k in cluster_key:
        dna_file_name = c_k + ".csv"
        dna_clone_list_file_name = c_k + "_label.csv"
        tmp_dna_label = pd.read_csv(cluster_dir + dna_file_name, index_col=0)
        tmp_dna_label = tmp_dna_label.set_index("cell")
        tmp_clone_list = pd.read_csv(cluster_dir + dna_clone_list_file_name, index_col=0)
        tmp_clone_list = list(tmp_clone_list["cluster_label"])

        tmp_patient = c_k.split("_")[1]
        tmp_dna_label["clone"] = [tmp_patient+"_clone"+str(x) for x in tmp_dna_label["clone"]]
        tmp_clone_list = [tmp_patient+"_clone"+str(x) for x in tmp_clone_list]
        tmp_key = (cluster_method, cluster_tran, tmp_patient, "None")
        dna_label[tmp_key] = {"dna_cluster":tmp_dna_label, "clone_list": tmp_clone_list}
    return dna_label

def random_predict_num(n_dna_ci, n_rna_cj, n_dna):
    """
    expected number of RNA cells mapped to the some clone when randomly select the DNA clone id to the RNA cells
    @param n_dna_ci: num of dna cells in clone ci
    @param n_rna_cj: num of rna cells in clone cj
    @param n_dna: num of all dna cells
    @return: expected number of RNA cells being mapped to clone ci, and those RNA cells' original clone id is cj
    """
    return round(n_dna_ci * n_rna_cj / n_dna, 1)


def random_baseline(cluster_tran, cluster_info, rna_folder, cluster_method,
                    patient_list=("crc04", "crc10", "crc11")):
    """
    Random baseline
    @param cluster_tran: string, normalisation method for the data being clustered,
        either "untransformed" or "log"
    @param cluster_info: dict, extra info for each clustering method
    @param rna_folder: string, the directory where the rna files are
    @param cluster_method: string, method used for clustering, either "Agglomerative" or "intNMF"
    @return:
    """
    rna_cell = get_rna_cells(rna_folder)
    if cluster_method == "Agglomerative":
        dna_cluster_dict = get_Agg_clusters(cluster_tran, cluster_info)
    elif cluster_method == "intNMF":
        dna_cluster_dict = get_intNMF_clusters(cluster_tran, cluster_info)
    else:
        sys.exit('Wrong cluster_method. It should be Agglomerative or intNMF.')

    baseline = {"rna_clone": [],
                "rna_number": [],
                "dna_clone": [],
                "dna_number": [],
                "dna_all_number": [],
                "expected_random_number": [],
                "setting":[]}

    for dna_key in dna_cluster_dict.keys():
        tmp_dna_cluster_df = dna_cluster_dict[dna_key]["dna_cluster"]
        tmp_clone_list = dna_cluster_dict[dna_key]["clone_list"]
        tmp_patient = dna_key[2]
        tmp_rna_cell = rna_cell[tmp_patient]
        tmp_dna_cell = list(tmp_dna_cluster_df.index)
        tmp_same_cell = set(tmp_rna_cell).intersection(tmp_dna_cell)
        tmp_rna_cluster_df = tmp_dna_cluster_df.loc[tmp_same_cell, :]
        n_dna = len(tmp_dna_cluster_df.index)
        tmp_setting = "_".join([dna_key[0], dna_key[1],dna_key[3]])
        for dna_clone in tmp_clone_list:
            for rna_clone in tmp_clone_list:
                n_dna_ci = len([x for x in list(tmp_dna_cluster_df["clone"]) if x == dna_clone])
                n_rna_cj = len([x for x in list(tmp_rna_cluster_df["clone"]) if x == rna_clone])
                n_expected = random_predict_num(n_dna_ci, n_rna_cj, n_dna)
                baseline["rna_clone"].append(rna_clone)
                baseline["dna_clone"].append(dna_clone)
                baseline["rna_number"].append(n_rna_cj)
                baseline["dna_number"].append(n_dna_ci)
                baseline["dna_all_number"].append(n_dna)
                baseline["expected_random_number"].append(n_expected)
                baseline["setting"].append(tmp_setting)
    baseline = pd.DataFrame.from_dict(baseline)
    return baseline



if __name__ == '__main__':
    patient_list = ["crc04", "crc10", "crc11"]
    cluster_tran = ["untransformed", "log"]
    cluster_method = ["Agglomerative", "intNMF"]
    rna_dir = "../data/2000_genes_log/"
    cluster_info = {"Agglomerative": {"cluster_key":
                                          {"untransformed": [str(x) + "_clusters" for x in range(6, 13)],
                                           "log": [str(x) + "_clusters" for x in range(6, 13)]},
                                      "cluster_dir": "../data/clusters/Agglomerative/",
                                      },
                    "intNMF": {"cluster_key":
                                   {"untransformed": ["untransformed_crc04_2", "untransformed_crc10_2",
                                                      "untransformed_crc11_2"],
                                    "log": ["log_crc04_2", "log_crc10_3", "log_crc11_3"]},
                               "cluster_dir": "../data/clusters/intNMF/",
                               }}
    result = []
    for c_method in cluster_method:
        for c_tran in cluster_tran:
            tmp_result = random_baseline(c_tran, cluster_info, rna_dir, c_method)
            result.append(tmp_result)
    result = pd.concat(result)
    resultDir = "Result/"
    result_filename = "random_baseline_result.csv"
    result.to_csv(resultDir + result_filename)
