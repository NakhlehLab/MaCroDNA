import numpy as np
import pandas as pd
import os
import copy
from sklearn.cluster import AgglomerativeClustering
from scipy.spatial.distance import cdist


def non_zero_var_gene(cell_df):
    cell_var = list(cell_df.var(axis=1))
    gene_all = list(cell_df.index)
    gene = [gene_all[i] for i in range(len(cell_var)) if cell_var[i] != 0]
    return gene


def non_zero_gene(cell_df, percentage=0.01):
    cell_amount = len(cell_df.columns)
    sum_thresh = cell_amount * percentage
    df_mask = copy.deepcopy(cell_df)
    df_mask = df_mask.mask(df_mask > 0, 1)
    df_mask = df_mask[(df_mask.sum(axis=1) > sum_thresh)]
    keep_gene = list(df_mask.index)
    return keep_gene


def find_genes(dna_df, rna_df, remove_var0=False,remove_0_percent=0.01):
    rna_gene = list(rna_df.index)
    dna_gene = list(dna_df.index)
    same_gene = set(rna_gene).intersection(dna_gene)

    #     remove 0
    rna_gene_filter_0 = non_zero_gene(rna_df, percentage=remove_0_percent)
    dna_gene_filter_0 = non_zero_gene(dna_df,percentage=remove_0_percent)
    result_gene = same_gene.intersection(rna_gene_filter_0)
    result_gene = result_gene.intersection(dna_gene_filter_0)

    if remove_var0:
        rna_gene_nonzero = non_zero_var_gene(rna_df)
        dna_gene_nonzero = non_zero_var_gene(dna_df)
        result_gene = result_gene.intersection(rna_gene_nonzero)
        result_gene = result_gene.intersection(dna_gene_nonzero)
    else:
        result_gene = result_gene
    #     return result_gene
    result_gene = list(result_gene)
    return dna_df.loc[result_gene, :], rna_df.loc[result_gene, :]


# add pseudo 1 to all data
def pseudo_1(df):
    return df.applymap(lambda x:x+1)


# remove data with any 0 copy number
def rm_0_gene(df):
    return df.loc[(df != 0).all(axis=1)]


# agglomerative clustering
class Agg_CLuster_Combine:
    def __init__(self, data_df, initial_cluster_num, clusterdir):
        self.data_df = data_df
        self.cluster_num = initial_cluster_num
        self.clusterDir = clusterdir

    def patient_cluster_cnv(self,dna_cluster, dna_df):
        cluster_id = list(set(list(dna_cluster["clone"])))
        for i in range(len(cluster_id)):
            iid = cluster_id[i]
            tmp_cluster_id = list(dna_cluster[dna_cluster["clone"] == iid]["cell"])
            tmp_cluster_df = dna_df[tmp_cluster_id]
            tmp_cluster_mean = tmp_cluster_df.mean(axis=1)

            if i == 0:
                cluster_cnv_df = tmp_cluster_mean.to_frame(name="clone" + str(iid))

            else:
                cluster_cnv_df["clone" + str(iid)] = tmp_cluster_mean

        return cluster_cnv_df

    def read_patient_cluster(self):
        # get initial clusters
        # Return: org_cluster_cnv = {cnv: cnv df row id are genes, col id are clone id
        #                                            cell_num: {clone0: num, clone1:num}
        #                                             cell2clone: df,2cols: clone and cell}
        org_cluster_cnv = {}
        cnv = self.data_df
        cell_id = list(cnv.columns)

        initial_result = AgglomerativeClustering(n_clusters=self.cluster_num,affinity="euclidean",linkage="ward").\
            fit_predict(cnv.to_numpy().T)
        initial_cluster_df = pd.DataFrame(initial_result, columns=["clone"])
        initial_cluster_df["cell"] = cell_id
        initial_cluster_cnv_df = self.patient_cluster_cnv(initial_cluster_df, cnv)

        cell_num = initial_cluster_df.groupby("clone")["cell"].count().to_dict()
        cell_num = dict(("clone" + str(key), value) for (key, value) in cell_num.items())

        org_cluster_cnv["cnv"] = initial_cluster_cnv_df
        org_cluster_cnv["cell2clone"] = initial_cluster_df
        org_cluster_cnv["cell_num"] = cell_num
        return org_cluster_cnv

    def least_dis_one_patient(self, patient_df):
        # least distance among clones in one patient_df, use euclidean distance
        # Input: patient_df -- dataframe, col ids are clone ids, row ids are genes
        # Output: min_c1, min_c2 -- two clones with minimum distance
        # tmp_dist / len(patient_df.index) -- minimum distance
        clone_ls = list(patient_df.columns)
        min_dis = 9999999
        min_c1 = "test"
        min_c2 = "test"
        for c1 in clone_ls:
            for c2 in clone_ls:
                if c1 != c2:
                    v1 = patient_df[c1].to_numpy()
                    v2 = patient_df[c2].to_numpy()
                    v1 = v1.reshape(1, -1)
                    v2 = v2.reshape(1, -1)
                    tmp_dist = np.sum(np.square(v1 - v2))
                    if tmp_dist < min_dis:
                        min_c1 = c1
                        min_c2 = c2
                        min_dis = tmp_dist

        return min_c1, min_c2, min_dis / len(patient_df.index)

    def patient_agg_one_turn(self, cluster_cnv, cluster_cell_num,agg_cluster_num):
        # one iteration for updating the clone based on the minimum dis from least_dis_one_patient

        # INPUT:
        # cluster_cnv -- dict, keys are patient id, values are cnv df with cols are clones, rows are genes
        # cluster_cell_num -- dict, keys are patient id, values are dict, keys are clones, cell number in that clone
        # patient_cluster_num -- dict, keys are patient id, values are clone number in that patient

        # OUTPUT:
        # new_cluster_cnv -- dict,updated cluster_cnv,keys are patient id, values are cnv df with cols are clones, rows are genes
        # new_cluster_cell_num -- dict, updated cell number, keys are patient id, values are cell number in that clone
        # clone_map -- dict, keys are clones of output, values are updated clones

        assert len(cluster_cnv.keys()) == len(cluster_cell_num.keys())
        patient_list = list(cluster_cnv.keys())
        min_dis = 99999
        min_c1 = "test"
        min_c2 = "test"
        clone_map = {}

        #     find least distance
        tmp_cluster_cnv = cluster_cnv
        min_c1, min_c2, min_dis = self.least_dis_one_patient(tmp_cluster_cnv)


        #     update the df in least distance
        min_cluster_cnv = cluster_cnv
        min_cluster_cell_num = cluster_cell_num
        #     get new cnv and cell num
        new_cell_cnv = min_cluster_cell_num[min_c1] * min_cluster_cnv[min_c1] + \
                       min_cluster_cell_num[min_c2] * min_cluster_cnv[min_c2]
        new_cell_num = min_cluster_cell_num[min_c1] + min_cluster_cell_num[min_c2]
        new_cell_cnv = new_cell_cnv / new_cell_num

        patient_new_cluster_cnv = new_cell_cnv.to_frame("clone0")

        patient_cluster_cell_num = {"clone0": new_cell_num}
        clone_map[min_c1] = "clone0"
        clone_map[min_c2] = "clone0"

        #     update cnv df and cell num dict
        i = 1
        for tmp_clone in min_cluster_cell_num.keys():
            if tmp_clone != min_c1 and tmp_clone != min_c2:
                patient_new_cluster_cnv["clone" + str(i)] = min_cluster_cnv[tmp_clone]
                patient_cluster_cell_num["clone" + str(i)] = min_cluster_cell_num[tmp_clone]
                clone_map[tmp_clone] = "clone" + str(i)
                i = i + 1

        new_patient_cluster_num = copy.deepcopy(patient_cluster_num)
        new_patient_cluster_num[min_patient] = new_patient_cluster_num[min_patient] - 1

        new_cluster_cnv = copy.deepcopy(cluster_cnv)
        new_cluster_cnv[min_patient] = patient_new_cluster_cnv

        new_cluster_cell_num = copy.deepcopy(cluster_cell_num)
        new_cluster_cell_num[min_patient] = patient_cluster_cell_num

        return new_cluster_cnv, new_cluster_cell_num, min_patient, clone_map, new_patient_cluster_num

    def save_cluster_data(self, original_cluster_cell2clone,
                          cluster_map, folder_path, agg_cluster_num):
        # save data for each iteration
        # INPUT:
        # original_cluster_cell2clone -- df, original clone for each cell, 2 cols, clone and cell
        # patient_cluster_num -- dict, keys are patient ids, values are cluster num in each patient_cluster_num
        # cluster_map -- dict, keys are patient ids, values are dict, keys are original clone id, values are current clone id
        # folder_path -- str, directory to save the data
        # agg_cluster_num -- int, current cluster numbers across all patients

        # OUPUT:
        # None, save the data

        save_cluster_cell2clone = {"cell": [], "original_clone": [], "agg_clone": []}
        save_patient_cluster = {"cluster_label": []}
        #     save_cluster_cell_num ={"original_clone":[],"agg_clone":[],"cell_amount":[]}
        tmp_cell2clone = original_cluster_cell2clone
        tmp_cell = list(tmp_cell2clone["cell"])
        tmp_original_clone = list(tmp_cell2clone["clone"])
        tmp_agg_clone = [str(cluster_map[x]) for x in tmp_original_clone]
        save_cluster_cell2clone["cell"].extend(tmp_cell)
        save_cluster_cell2clone["original_clone"].extend(tmp_original_clone)
        save_cluster_cell2clone["agg_clone"].extend(tmp_agg_clone)
        save_cluster_cell2clone = pd.DataFrame.from_dict(save_cluster_cell2clone)
        save_cluster_cell2clone.to_csv(folder_path + str(agg_cluster_num) + "_clusters.csv", header=True)


    def patient_agg_all_turn(self,cluster_cnv, cluster_cell_num, original_cluster_cell2clone,
                             tmp_norm_method):
        # for one norm method, get the cluster resolution from self.cluster_num to 2 clusters in each patient
        # INPUT:
        # cluster_cnv -- dict, keys are patient id, values are cnv df with cols are clones, rows are genes
        # cluster_cell_num -- dict, keys are patient id, values are dict, keys are clones, cell number in that clone
        # original_cluster_cell2clone -- dict, keys are patient id, values are  df, original clone for each cell, 2 cols, clone and cell

        # OUTPUT:
        # None, save the agg clusters

        num_cluster = self.cluster_num
        start_clone_num = num_cluster
        end_clone_num = 2
        agg_cluster_num = copy.deepcopy(start_clone_num)
        new_cluster_cnv = copy.deepcopy(cluster_cnv)
        new_cluster_cell_num = copy.deepcopy(cluster_cell_num)

        #     save 4 clusters
        folder_path = self.clusterDir+"agglomerative_"+tmp_norm_method+"/"

        #     make result folder
        if not os.path.isdir(folder_path):
            os.mkdir(folder_path)

        #     original cluster_map
        cluster_map = {}
        for n_c in range(num_cluster):
            cluster_map[n_c] = "clone" + str(n_c)


        #     save 12 clusters
        self.save_cluster_data(original_cluster_cell2clone,
                               cluster_map, folder_path, agg_cluster_num)

        while agg_cluster_num > end_clone_num:
            agg_cluster_num = agg_cluster_num - 1
            new_cluster_cnv, new_cluster_cell_num, min_patient, min_patient_cluster_map, patient_cluster_num = \
                self.patient_agg_one_turn(new_cluster_cnv, new_cluster_cell_num, patient_cluster_num)
            # print("agg_cluster_num", agg_cluster_num, " ,patient_cluster_num", patient_cluster_num)

            #         update cluster map
            tmp_patient_cluster_map = cluster_map[min_patient]
            new_patient_cluster_map = {}

            for p_keys in tmp_patient_cluster_map.keys():
                new_patient_cluster_map[p_keys] = min_patient_cluster_map[tmp_patient_cluster_map[p_keys]]
            cluster_map[min_patient] = new_patient_cluster_map
            # print("Clusters: ", agg_cluster_num)
            # print("patient: ", min_patient)
            # print("map: ", cluster_map)

            #         save data
            self.save_cluster_data(original_cluster_cell2clone, patient_cluster_num,
                              cluster_map, folder_path, agg_cluster_num)



    def agg_cluster_combine(self):
        org_cluster_cnv = self.read_patient_cluster()
        self.patient_agg_all_turn(org_cluster_cnv["cnv"],
                                  org_cluster_cnv["cell_num"],
                                  org_cluster_cnv["cell2clone"])


        # for tmp_norm_method in self.cluster_norm:
        #     norm_org_cluster_cnv = org_cluster_cnv[tmp_norm_method]
        #     cluster_cnv = {}
        #     cluster_cell_num = {}
        #     original_cluster_cell2clone = {}
        #
        #     for tmp_patient in norm_org_cluster_cnv.keys():
        #         cluster_cnv[tmp_patient] = norm_org_cluster_cnv[tmp_patient]["cnv"]
        #         cluster_cell_num[tmp_patient] = norm_org_cluster_cnv[tmp_patient]["cell_num"]
        #         original_cluster_cell2clone[tmp_patient] = norm_org_cluster_cnv[tmp_patient]["cell2clone"]
        #
        #     self.patient_agg_all_turn(cluster_cnv, cluster_cell_num, original_cluster_cell2clone, tmp_norm_method)


def cluster_distance(data_df):
    min_dis = float('inf')
    cluster_data = {}
    for cluster, group in data_df.groupby('clone'):
        cluster_data[cluster] = group.drop('clone', axis=1)
    cluster_id = list(cluster_data.keys())
    cluster_num = len(cluster_id)
    min_c1, min_c2 = None, None

    for i in range(cluster_num):
        for j in range(i + 1, cluster_num):
            c1_data = cluster_data[cluster_id[i]]
            c2_data = cluster_data[cluster_id[j]]
            dis = cdist(c1_data, c2_data, "euclidean")
            dis_sum = dis.sum() / (len(c1_data) + len(c2_data))

            if dis_sum < min_dis:
                min_dis = dis_sum
                min_c1, min_c2 = cluster_id[i], cluster_id[j]

    return min_c1, min_c2, min_dis




# agglomerative cluster and merge
def cluster_and_merge(data_df, cluster_num, clusterDir):
    if not os.path.isdir(clusterDir):
        os.mkdir(clusterDir)

    # Step 1: Perform agglomerative clustering to find four cell clusters
    data_df = data_df.T
    clustering = AgglomerativeClustering(n_clusters=cluster_num,metric="euclidean",linkage="ward")
    cluster_labels = clustering.fit_predict(data_df)

    # Save the cluster results
    data_df['clone'] = cluster_labels
    cluster_result = pd.DataFrame(list(zip(data_df.index.tolist(),
                                           ['clone_'+str(label) for label in data_df['clone']])),
                                  columns=["cell", "clone"])
    cluster_result.to_csv(clusterDir+'4_clusters.csv', index=False)

    # Step 2: Calculate distances between each pair of clusters and merge the two with the least distance
    while cluster_num > 2:
        min_c1, min_c2,min_dis = cluster_distance(data_df)

        data_df['clone'] = np.where(data_df['clone'] == min_c2, min_c1, data_df['clone'])
        data_df['clone'] = np.where(data_df['clone'] > min_c2, data_df['clone'] - 1, data_df['clone'])
        cluster_num -= 1

        cluster_result = pd.DataFrame(list(zip(data_df.index.tolist(),
                                               ['clone_' + str(label) for label in data_df['clone']])),
                                      columns=["cell", "clone"])
        # Save the cluster results
        cluster_result.to_csv(clusterDir + str(cluster_num) +'_clusters.csv', index=False)
