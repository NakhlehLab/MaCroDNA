import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix,accuracy_score


def read_intNMF_result(cluster_key, cluster_dir, cluster_tran, intgr_dir,
                       intgr_methods=("clonealign","Seurat","MaCroDNA"),
                       patient_list=("crc04", "crc10", "crc11")):
    """
    Prepare the intNMF integration results
    @param cluster_key: dict, information for each cluster's name
    @param cluster_dir: string, specify where the cluster files are
    @param cluster_tran: string, normalisation method for the data being clustered,
        either "untransformed" or "log"
    @param intgr_dir: string, specify where the integration results are
    @param intgr_methods: list, integration methods
    @return: intgr_result: dict, {(intgr_method,cluster_method,cluster_tran,patient, cluster_key): result_df},
    result_df is a df with 2 cols: "label", "predict"
    """


    # read the integration results from the integration methods
    intgr_result = {}
    for c_k in cluster_key:
        # read dna clusters
        dna_file_name = c_k + ".csv"
        tmp_dna_label = pd.read_csv(cluster_dir+dna_file_name, index_col=0)
        tmp_dna_label = tmp_dna_label.set_index("cell")


        # read rna results
        for tmp_intgr_method in intgr_methods:
            rna_file_name = "_".join([tmp_intgr_method, "intNMF", c_k, "result.csv"])
            tmp_rna_result = pd.read_csv(intgr_dir + rna_file_name, index_col=0)
            tmp_rna_result = tmp_rna_result.set_index("cell")

            # transfer true label from dna label to rna result
            tmp_same_cell = list(set(list(tmp_rna_result.index)).intersection(list(tmp_dna_label.index)))
            tmp_rna_result_same = tmp_rna_result.loc[tmp_same_cell, :]
            tmp_dna_label_same = tmp_dna_label.loc[tmp_same_cell, :]
            tmp_result = pd.DataFrame(list(zip(list(tmp_dna_label_same["clone"]),
                                               list(tmp_rna_result_same["predict"]))),
                                      columns=["label", "predict"])
            assert cluster_tran == c_k.split("_")[0]
            intgr_result[(tmp_intgr_method, "intNMF", cluster_tran,c_k.split("_")[1],c_k)] = tmp_result

    return intgr_result



def read_Agg_result(cluster_key, cluster_dir, cluster_tran, intgr_dir,
                    intgr_methods=("clonealign","Seurat","MaCroDNA"),
                    patient_list=("crc04", "crc10", "crc11")):
    """
    Prepare the agglomerative integration results
    @param cluster_key: dict, information for each cluster's name
    @param cluster_dir: string, specify where the cluster files are
    @param cluster_tran: string, normalisation method for the data being clustered,
        either "untransformed" or "log"
    @param intgr_dir: string, specify where the integration results are
    @param intgr_methods: list, integration methods
    @return: intgr_result: dict, {(intgr_method,cluster_method,cluster_tran,patient, cluster_key): result_df},
    result_df is a df with 2 cols: "label", "predict"
    """

    # read the integration results from the integration methods
    intgr_result = {}
    for c_k in cluster_key:
        # read dna clusters
        dna_file_name = c_k + ".csv"
        tmp_dna_label = pd.read_csv(cluster_dir + dna_file_name, index_col=0)
        tmp_dna_label = tmp_dna_label.set_index("cell")

        # read rna results
        for tmp_patient in patient_list:
            # obtain the mapping from the original_clone to agg_clone
            # original_clone is the clone id when inferring 4 clusters in each patient_list
            # agg_clone is the clone id after assembling all clusters throughout all the patients
            tmp_clone_map = {}
            for tmp_cell in list(tmp_dna_label.index):
                tmp_original_clone = tmp_dna_label.loc[tmp_cell, "original_clone"]
                tmp_agg_clone = tmp_dna_label.loc[tmp_cell, "agg_clone"]
                if tmp_agg_clone.split("_")[0] ==tmp_patient:
                    if tmp_original_clone not in list(tmp_clone_map.keys()):
                        tmp_clone_map[tmp_original_clone] = tmp_agg_clone
                    else:
                        assert tmp_clone_map[tmp_original_clone] == tmp_agg_clone

            for tmp_intgr_method in intgr_methods:
                rna_file_name = "_".join([tmp_intgr_method, "Agglomerative", cluster_tran, tmp_patient, "4_result.csv"])
                tmp_rna_result = pd.read_csv(intgr_dir + rna_file_name, index_col=0)
                tmp_rna_result = tmp_rna_result.set_index("cell")

                # transfer true label from dna label to rna result
                tmp_same_cell = list(set(list(tmp_rna_result.index)).intersection(list(tmp_dna_label.index)))
                tmp_rna_result_same = tmp_rna_result.loc[tmp_same_cell, :]
                tmp_dna_label_same = tmp_dna_label.loc[tmp_same_cell, :]

                tmp_rna_result_same = tmp_rna_result_same.replace({"predict": tmp_clone_map})
                tmp_result = pd.DataFrame(list(zip(list(tmp_dna_label_same["agg_clone"]),
                                                   list(tmp_rna_result_same["predict"]))),
                                          columns=["label", "predict"])
                intgr_result[(tmp_intgr_method, "Agglomerative", cluster_tran, tmp_patient, c_k)] = tmp_result
    return intgr_result


def read_integration_result(cluster_method, cluster_tran, cluster_info, intgr_dir):
    """
    Prepare the integration results for the further analysis
    @param cluster_method: string, method used for clustering, either "Agglomerative" or "intNMF"
    @param cluster_tran: string, normalisation method for the data being clustered,
        either "untransformed" or "log"
    @param cluster_info: dict, extra info for each clustering method
    @param intgr_dir: string, specify where the integration results are
    @return: intgr_result: dict, {(intgr_method,cluster_tran,patient, cluster_key): result_df},
    result_df is a df with 2 cols: "label", "predict"
    """
    cluster_key = cluster_info[cluster_method]["cluster_key"][cluster_tran]
    cluster_dir = cluster_info[cluster_method]["cluster_dir"] + cluster_tran + "/"

    if cluster_method == "Agglomerative":
        intgr_result = read_Agg_result(cluster_key, cluster_dir, cluster_tran, intgr_dir)

    elif cluster_method == "intNMF":
        intgr_result = read_intNMF_result(cluster_key, cluster_dir, cluster_tran, intgr_dir)

    else:
        sys.exit('Wrong cluster_method. It should be Agglomerative or intNMF.')

    return intgr_result


def plot_results_barplot(intgr_result,one_intgr_method, cluster_method, cluster_tran,cluster_info,
             patient_list=("crc04", "crc10", "crc11"), MaCroDNA_gene=None):
    """
    Barplots for the integration results
    @param intgr_result: intgr_result: dict, {(intgr_method,cluster_method,cluster_tran,patient, cluster_key): result_df},
    result_df is a df with 2 cols: "label", "predict"
    @param one_intgr_method: string, integration method, should be one from clonealign, Seurat, MaCroDNA
    @param cluster_method: string, method used for clustering, either "Agglomerative" or "intNMF"
    @param cluster_tran: string, normalisation method for the data being clustered,
        either "untransformed" or "log"
    @param cluster_info: dict, extra info for each clustering method
    @param MaCroDNA_gene: string, should be one from 2000_genes_log, all_genes_log, all_genes_raw, noX_genes_raw
    @return: N/A
    """

    # prepare data
    cluster_key = cluster_info[cluster_method]["cluster_key"][cluster_tran]
    plot_index = []
    plot_right = []
    plot_wrong = []

    if cluster_method == "Agglomerative":
        plot_index = cluster_key
        for c_k in cluster_key:
            tmp_right = 0
            tmp_wrong = 0
            for tmp_patient in patient_list:
                tmp_key = (one_intgr_method, "Agglomerative", cluster_tran, tmp_patient, c_k)
                tmp_result = intgr_result[tmp_key]
                r = sum(tmp_result["predict"] == tmp_result["label"])
                w = len(tmp_result.index) - r
                tmp_right = tmp_right + r
                tmp_wrong = tmp_wrong + w
            plot_right.append(tmp_right)
            plot_wrong.append(tmp_wrong)

    elif cluster_method == "intNMF":
        plot_index = [tmp_p.upper() for tmp_p in patient_list]
        for c_k in cluster_key:
            tmp_patient = c_k.split("_")[1]
            tmp_key = (one_intgr_method, "intNMF", cluster_tran, tmp_patient, c_k)
            tmp_result = intgr_result[tmp_key]
            tmp_right = sum(tmp_result["predict"] == tmp_result["label"])
            tmp_wrong = len(tmp_result.index) - tmp_right
            plot_right.append(tmp_right)
            plot_wrong.append(tmp_wrong)
    else:
        sys.exit('Wrong cluster_method. It should be Agglomerative or intNMF.')

    font_color = '#525252'
    hfont = {'fontname': 'DejaVu Sans'}
    facecolor = "white"
    color_red = "steelblue"
    color_blue = "sandybrown"
    index = plot_index

    column0 = plot_wrong
    column1 = plot_right
    title0 = 'incorrect predictions'
    title1 = 'correct predictions'

    col_font = 18
    title_font = 18
    axlabel_font = 15

    fig, axes = plt.subplots(figsize=(8, 5), facecolor=facecolor, ncols=2, sharey=True)

    axes[0].barh(index, column0, align='center', color=color_red, zorder=10)
    axes[0].set_title(title0, fontsize=col_font, pad=15, color=color_red, **hfont)
    axes[1].barh(index, column1, align='center', color=color_blue, zorder=10)
    axes[1].set_title(title1, fontsize=col_font, pad=15, color=color_blue, **hfont)

    axes[0].set(yticks=index, yticklabels=index)
    axes[0].yaxis.tick_left()
    axes[0].tick_params(axis='y', colors='white')  # tick color

    axes[1].set_xticks(list(np.arange(0, 360, 50)))
    axes[1].set_xticklabels(list(np.arange(0, 360, 50)))

    axes[0].set_xticks(list(np.arange(0, 360, 50)))
    axes[0].set_xticklabels(list(np.arange(0, 360, 50)))

    axes[0].invert_xaxis()

    for label in (axes[0].get_xticklabels() + axes[0].get_yticklabels()):
        label.set(fontsize=axlabel_font, color=font_color, **hfont)
    for label in (axes[1].get_xticklabels() + axes[1].get_yticklabels()):
        label.set(fontsize=axlabel_font, color=font_color, **hfont)

    for i1, v1 in enumerate(column0):
        axes[0].text(v1 + 1, i1 - 0.1, str(v1), color=color_red, fontweight='bold', ha='right')
    for i2, v2 in enumerate(column1):
        axes[1].text(v2 + 1, i2 - 0.1, str(v2), color=color_blue, fontweight='bold', ha='left')

    plt.subplots_adjust(wspace=0, top=0.80, bottom=0.1, left=0.18, right=0.95)

    if one_intgr_method != "MaCroDNA":
        fig_title = "[Method] " + one_intgr_method + "\n[Clone method] " + cluster_method + \
                    "\n[Clone transformation] " + cluster_tran
    else:
        fig_title = "[Method] " + one_intgr_method + "\n[Clone method] " + cluster_method + \
                    "\n[Clone transformation] " + cluster_tran + "\n[Gene] " + MaCroDNA_gene
    plt.suptitle(fig_title,
                 fontsize=title_font, color="black", horizontalalignment="left",
                 fontweight="bold", x=0.1)
    fig.tight_layout(w_pad=-1)
    plt.show()


def prepare_random_results(random_df,cluster_method, cluster_tran,
                           patient_list=("crc04", "crc10", "crc11")):
    """
    Return the accuracy and the clonal prevelance of the random baseline
    @param random_df: df, 7 cols, rna_clone, rna_number, dna_clone, dna_number, dna_all_number,
    expected_random_number, setting
    @return: acc, dict, key is the patient id, value is the accuracy
    clone_prevalence, df, 3 cols, "clone_id","true_percentage","predict_percentage"
    """

    # select data based on cluster_method and cluster_tran
    select_random_df = random_df[random_df["setting"].str.contains(cluster_method+"_"+cluster_tran)]
    if cluster_method == "Agglomerative":
        select_random_df = random_df[random_df["setting"].str.contains("_".join([cluster_method, cluster_tran,
                                                                                 "12_clusters"]))]
    # generate confusion matrix
    clone_label = list(set(select_random_df["dna_clone"]))
    cf = np.empty([len(clone_label), len(clone_label)])
    for i_dna in range(len(clone_label)):
        for i_rna in range(len(clone_label)):
            dna_clone = clone_label[i_dna]
            rna_clone = clone_label[i_rna]
            expected_num = select_random_df[(select_random_df["rna_clone"] == rna_clone) &
                                            (select_random_df["dna_clone"] == dna_clone)]["expected_random_number"]
            expected_num = list(expected_num)
            if rna_clone[:5] != dna_clone[:5]:
                expected_num = [0]
            assert len(expected_num) == 1
            expected_num = expected_num[0]
            expected_num = round(expected_num, 1)
            cf[i_rna, i_dna] = expected_num


    # seperate by patient
    acc = {}
    clone_prevalence = pd.DataFrame(columns=["clone_id","true_percentage","predict_percentage"])
    for tmp_patient in patient_list:
        tmp_clone_label = [x for x in clone_label if x[:5] == tmp_patient]
        tmp_c_label_index = [clone_label.index(x) for x in tmp_clone_label]
        tmp_cf = cf[tmp_c_label_index, :]
        tmp_cf = tmp_cf[:, tmp_c_label_index]
        tmp_accuracy = tmp_cf.trace()/tmp_cf.sum()
        acc[tmp_patient] = tmp_accuracy

        true_cluster_num = np.sum(tmp_cf, axis=1)
        predict_cluster_num = np.sum(tmp_cf, axis=0)
        true_cluster_per = true_cluster_num / sum(true_cluster_num)
        predict_cluster_per = predict_cluster_num / sum(predict_cluster_num)

        for clone_i in range(len(tmp_clone_label)):
            clone_prevalence = clone_prevalence.append({"clone_id": tmp_clone_label[clone_i],
                                                        "true_percentage": true_cluster_per[clone_i],
                                                        "predict_percentage": predict_cluster_per[clone_i]},
                                                       ignore_index=True)
    return acc, clone_prevalence


def prepare_intgr_plot_info(intgr_result, cluster_method, cluster_tran, cluster_info,
                            intgr_methods=("clonealign","Seurat","MaCroDNA"),
                            patient_list=("crc04", "crc10", "crc11")):
    cluster_dir = cluster_info[cluster_method]["cluster_dir"] + cluster_tran + "/"
    result_key = intgr_result.keys()
    acc = {}
    clone_prevalence = {}
    for tmp_method in intgr_methods:
        acc[tmp_method] = {}
        clone_prevalence[tmp_method] = pd.DataFrame(columns=["clone_id","true_percentage","predict_percentage"])
        for tmp_patient in patient_list:
            acc[tmp_method][tmp_patient] = 0

    for result_k in result_key:
        c_k = result_k[4]
        tmp_patient = result_k[3]
        tmp_method = result_k[0]
        if cluster_method == "intNMF" or c_k == "12_clusters":
            dna_clone_list_file_name = c_k + "_label.csv"
            tmp_clone_list = pd.read_csv(cluster_dir + dna_clone_list_file_name, index_col=0)
            tmp_clone_list = list(tmp_clone_list["cluster_label"])
            tmp_patient_clone_list = [x for x in tmp_clone_list if x[:5] == tmp_patient]
            tmp_data = intgr_result[result_k]
            data_true = list(tmp_data["label"])
            data_predict = list(tmp_data["predict"])
            tmp_cf = confusion_matrix(data_true, data_predict, labels=tmp_patient_clone_list)
            tmp_acc = accuracy_score(data_true, data_predict)
            acc[tmp_method][tmp_patient] = tmp_acc

            true_cluster_num = np.sum(tmp_cf, axis=1)
            predict_cluster_num = np.sum(tmp_cf, axis=0)
            true_cluster_per = true_cluster_num / sum(true_cluster_num)
            predict_cluster_per = predict_cluster_num / sum(predict_cluster_num)
            for clone_i in range(len(tmp_patient_clone_list)):
                clone_prevalence[tmp_method] = clone_prevalence[tmp_method].append({"clone_id":
                                                                                        tmp_patient_clone_list[clone_i],
                                                                                    "true_percentage":
                                                                                        true_cluster_per[clone_i],
                                                                                    "predict_percentage":
                                                                                        predict_cluster_per[clone_i]},
                                                                                   ignore_index=True)
    return acc, clone_prevalence


def plot_results_accuracy(intgr_result, cluster_method, cluster_tran, cluster_info,
                          random_path,
                          intgr_methods=("Baseline","Seurat","clonealign","MaCroDNA"),
                          patient_list=("crc04", "crc10", "crc11"), MaCroDNA_gene=None,
                          trans=0.4, scale=200, axis_label_font=16,axis_tick_font = 14,
                          pad=5,panel_font = 17):
    """
    Accuracy plots
    @param intgr_result: intgr_result: dict, {(intgr_method,cluster_method,cluster_tran,patient, cluster_key): result_df},
    result_df is a df with 2 cols: "label", "predict"
    @param cluster_method: string, method used for clustering, either "Agglomerative" or "intNMF"
    @param cluster_tran: string, normalisation method for the data being clustered,
        either "untransformed" or "log"
    @param cluster_info: dict, extra info for each clustering method
    @param random_path: string, directory+file_name of the random baseline
    @return: N/A
    """

    model_marker = {"Baseline": "^", "Seurat": "X", "clonealign": "P",
                    "MaCroDNA": "o"}
    model_color = {"Baseline": "darkcyan", "Seurat": "yellowgreen", "clonealign": "gray",
                   "MaCroDNA": "darkviolet"}

    intgr_acc,_ = prepare_intgr_plot_info(intgr_result, cluster_method, cluster_tran, cluster_info)

    fig, ax = plt.subplots(nrows=1, ncols=1)
    for model_i in range(len(intgr_methods)):
        tmp_model = intgr_methods[model_i]
        if tmp_model == "Baseline":
            random_model = pd.read_csv(random_path, index_col=0)
            acc, _ = prepare_random_results(random_model, cluster_method, cluster_tran)
        else:
            acc = intgr_acc[tmp_model]

        acc_array = list(acc.values())
        ax.scatter([model_i] * 3, acc_array, color=model_color[tmp_model],
                   label=tmp_model,
                   marker=model_marker[tmp_model],
                   alpha=trans, s=scale)

    ax.set_ylabel("Accuracy",fontsize=axis_label_font)
    ax.set(ylim=(-0.05, 1))
    ax.set_xticks([0, 1, 2, 3])
    ax.set_xticklabels(intgr_methods, rotation=30)

    ax.xaxis.set_tick_params(labelsize=axis_tick_font)
    ax.yaxis.set_tick_params(labelsize=axis_tick_font)
    ax.annotate(cluster_tran, xy=(0.5, 1.02), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=panel_font, ha='center', va='baseline',
                weight="bold")
    ax.annotate(cluster_method, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                xycoords=ax.yaxis.label, textcoords='offset points',
                size=panel_font, ha='right', va='center',
                weight="bold", rotation="vertical")
    fig.subplots_adjust(left=0.2, top=0.9, bottom=0.2)
    plt.show()


def plot_results_clone_prevalence(intgr_result, cluster_method, cluster_tran, cluster_info,
                                  random_path,
                                  intgr_methods=("Baseline","Seurat","clonealign","MaCroDNA"),
                                  patient_list=("crc04", "crc10", "crc11"), MaCroDNA_gene=None,
                                  trans=0.4, scale=200, axis_label_font=16,axis_tick_font = 14,
                                  pad=5,panel_font = 17, legend_font = 14):
    """
    Clonal prevalence plots
    @param intgr_result: intgr_result: dict, {(intgr_method,cluster_method,cluster_tran,patient, cluster_key): result_df},
    result_df is a df with 2 cols: "label", "predict"
    @param cluster_method: string, method used for clustering, either "Agglomerative" or "intNMF"
    @param cluster_tran: string, normalisation method for the data being clustered,
        either "untransformed" or "log"
    @param cluster_info: dict, extra info for each clustering method
    @param random_path: string, directory+file_name of the random baseline
    @return: N/A
    """

    model_marker = {"Baseline": "^", "Seurat": "X", "clonealign": "P",
                    "MaCroDNA": "o"}
    model_color = {"Baseline": "darkcyan", "Seurat": "yellowgreen", "clonealign": "gray",
                   "MaCroDNA": "darkviolet"}

    _, intgr_clone = prepare_intgr_plot_info(intgr_result, cluster_method, cluster_tran, cluster_info)

    fig, ax = plt.subplots(nrows=1, ncols=1)
    for model_i in range(len(intgr_methods)):
        tmp_model = intgr_methods[model_i]
        if tmp_model == "Baseline":
            random_model = pd.read_csv(random_path, index_col=0)
            _, clone_pre = prepare_random_results(random_model, cluster_method, cluster_tran)
        else:
            clone_pre = intgr_clone[tmp_model]

        ax.scatter(clone_pre["true_percentage"],
                   clone_pre["predict_percentage"],
                   c=model_color[tmp_model], label=tmp_model,
                   marker=model_marker[tmp_model],
                   alpha=trans, s=scale)
    ax.plot([0, 1], [0, 1], alpha=0.1, c="gray")
    ax.set_ylabel("Prediction", fontsize=axis_label_font)
    ax.set_xlabel("Ground truth", fontsize=axis_label_font)
    ax.set(ylim=(-0.05, 1.05))
    ax.set(xlim=(-0.05, 1.05))


    ax.xaxis.set_tick_params(labelsize=axis_tick_font)
    ax.yaxis.set_tick_params(labelsize=axis_tick_font)
    ax.annotate(cluster_tran, xy=(0.5, 1.02), xytext=(0, pad),
                xycoords='axes fraction', textcoords='offset points',
                size=panel_font, ha='center', va='baseline',
                weight="bold")
    ax.annotate(cluster_method, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                xycoords=ax.yaxis.label, textcoords='offset points',
                size=panel_font, ha='right', va='center',
                weight="bold", rotation="vertical")

    axLine, axLabel = ax.get_legend_handles_labels()
    fig.legend(loc='lower center',fancybox=True, shadow=True, ncol=4)
    # fig.legend(axLine, axLabel, loc='lower center', bbox_to_anchor=(0.5, 0.01),
    #            fancybox=True, shadow=True, ncol=4, fontsize=legend_font)
    fig.subplots_adjust(left=0.2, top=0.9, bottom=0.25, right=0.9)
    plt.show()





if __name__ == '__main__':

    patient_list = ["CRC04", "CRC10", "CRC11"]
    method_list = ["clonealign", "Seurat", "MaCroDNA"]
    agglomerative_cluster_key = [str(x)+"_clusters" for x in range(6, 13)]

    cluster_info = {"Agglomerative": {"cluster_key":
                                      {"untransformed": [str(x)+"_clusters" for x in range(6, 13)],
                                      "log": [str(x)+"_clusters" for x in range(6, 13)]},
                                     "cluster_dir": "../data/clusters/Agglomerative/",
                                     },
                   "intNMF": {"cluster_key":
                              {"untransformed": ["untransformed_crc04_2", "untransformed_crc10_2","untransformed_crc11_2"],
                              "log": ["log_crc04_2", "log_crc10_3","log_crc11_3"]},
                             "cluster_dir": "../data/clusters/intNMF/",
                             }}

    """
    Specify the cluster method and cluster transformation method
    for infering clones from the CNV (scDNA-seq) data
    
    ----------
    Variables
    ----------
    CLUSTER_METHOD: string, method used for getting clusters, either "Agglomerative" or "intNMF"
    CLUSTER_TRANS: string, the transformation used for the data before clustering, either "untransformed" or "log"
    INTGR_DIR: string, where the integration results are
    
    """

    CLUSTER_METHOD = "Agglomerative"
    CLUSTER_TRANS = "untransformed"
    INTGR_DIR = "Result/"
    MaCroDNA_GENE = "2000_genes_log"
    RANDOM_PATH = "Result/random_baseline_result.csv"

    assert CLUSTER_METHOD == "Agglomerative" or CLUSTER_METHOD == "intNMF"
    assert CLUSTER_TRANS == "untransformed" or CLUSTER_TRANS == "log"

    intgr_results = read_integration_result(CLUSTER_METHOD, CLUSTER_TRANS, cluster_info, INTGR_DIR)

    plot_results_barplot(intgr_results, "MaCroDNA", CLUSTER_METHOD, CLUSTER_TRANS, cluster_info,
                               MaCroDNA_gene=MaCroDNA_GENE)

    plot_results_accuracy(intgr_results, CLUSTER_METHOD, CLUSTER_TRANS, cluster_info, RANDOM_PATH)
    plot_results_clone_prevalence(intgr_results, CLUSTER_METHOD, CLUSTER_TRANS, cluster_info, RANDOM_PATH)