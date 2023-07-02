import os

import numpy as np
from scipy import stats
import gurobipy as gp
from gurobipy import *
import pandas as pd
from scipy import spatial
from scipy.stats import spearmanr
from numpy.linalg import norm
import math
from math import nan
import multiprocessing
import time
import random

class MaCroDNA:
    def __init__(self, rna_df=None, dna_df=None, dna_label=None):
        # INPUT
        # rna_df, dna_df -- dataframe, cols are cell ids, index are genes
        # dna_label -- dataframe, two cols, one is the clone id "clone", the other is the cell id "cell"
        self.rna_df = rna_df
        self.dna_df = dna_df
        self.dna_label = dna_label

    def cosine_similarity(self, x1, x2):
        return 1 - spatial.distance.cosine(x1 - x1.mean(), x2 - x2.mean())

    def cosine_similarity_np(self, x1, x2):
        return np.dot(x1 - x1.mean(), x2 - x2.mean()) / (1e-10 + norm(x1 - x1.mean()) * norm(x2 - x2.mean()))

    def ilp(self, rna_idx, dna_idx, corrs):

        n_min = min(rna_idx.shape[0], dna_idx.shape[0])
        # print("the smallest set has %s number of cells" % (n_min))
        # create the Gruobi model

        m = gp.Model("Pearson_GRB")

        # create a 2-D array of binary variables
        # x[i,j]=1 means that cell i in RNA corresponds to cell j in DNA
        x = []
        for i in range(rna_idx.shape[0]):
            x.append([])
            for j in range(dna_idx.shape[0]):
                x[i].append(m.addVar(vtype=GRB.BINARY, name="x[%d,%d]" % (i, j)))

        # set the contraints for rows and columns

        for i in range(rna_idx.shape[0]):
            # At most one cell per row
            m.addConstr(quicksum([x[i][j] for j in range(dna_idx.shape[0])]) <= 1, name="row" + str(i))

        for j in range(dna_idx.shape[0]):
            # At most one cell per column
            m.addConstr(quicksum([x[i][j] for i in range(rna_idx.shape[0])]) <= 1, name="col" + str(i))

        m.addConstr(quicksum([x[i][j] for i in range(rna_idx.shape[0]) for j in range(dna_idx.shape[0])]) == n_min,
                    name="global")

        # update the model
        m.update()

        # create the linear expression
        obj = LinExpr()

        for i in range(rna_idx.shape[0]):
            for j in range(dna_idx.shape[0]):
                obj += x[i][j] * corrs[rna_idx[i]][dna_idx[j]]

        # set the objective
        m.setObjective(obj, GRB.MAXIMIZE)

        # optimize the model
        m.optimize()

        print('IsMIP: %d' % m.IsMIP)
        if m.status == GRB.Status.INFEASIBLE:
            print("The model is infeasible")
        print("Solved with MIPFocus: %d" % m.Params.MIPFocus)
        print("The model has been optimized")
        print('Obj: %g' % m.objVal)

        results = np.empty((rna_idx.shape[0], dna_idx.shape[0]))
        for i in range(rna_idx.shape[0]):
            for j in range(dna_idx.shape[0]):
                results[i][j] = int(x[i][j].x)

        return results

    def cell2cell_assignment(self):
        dna_cells = list(self.dna_df.columns)
        rna_cells = list(self.rna_df.columns)
        genes = set(self.dna_df.index).intersection(self.rna_df.index.to_list())
        self.dna_df = self.dna_df.loc[genes, :]
        self.rna_df = self.rna_df.loc[genes, :]

        dna_np = self.dna_df.T.to_numpy()
        rna_np = self.rna_df.T.to_numpy()

        # create the array containing Pearson correlation coefficients between all possible pairs

        # the rows of the matrix correspond the rna cells and the columns correspond the dna cells
        corrs = np.empty((rna_np.shape[0], dna_np.shape[0]))

        for i in range(rna_np.shape[0]):
            for j in range(dna_np.shape[0]):
                corrs[i][j] = self.cosine_similarity_np(rna_np[i], dna_np[j])

        # create a matrix for storing the global correspondences
        global_correspondence = np.zeros((rna_np.shape[0], dna_np.shape[0]))
        # create a matrix of correspondence where the entries represent the steps
        tagged_correspondence = np.zeros((rna_np.shape[0], dna_np.shape[0]))
        # create index lists for dna and rna
        global_dna_idx = np.array([i for i in range(dna_np.shape[0])])
        global_rna_idx = np.array([i for i in range(rna_np.shape[0])])

        # calculate the number of steps
        quotient, remainder = divmod(global_rna_idx.shape[0],
                                     global_dna_idx.shape[0])
        
        n_iters = int(quotient)
        if remainder != 0:
            n_iters += 1

        rna_idx = np.copy(global_rna_idx)

        # iterations

        for step_ in range(n_iters):

            r = self.ilp(rna_idx, global_dna_idx, corrs)
            # identify the indices with assignments
            r = np.array(r)
            for i in range(r.shape[0]):
                for j in range(r.shape[1]):
                    if r[i][j] == 1:
                        global_correspondence[rna_idx[i]][global_dna_idx[j]] = 1  # binary correspondence matrix
                        tagged_correspondence[rna_idx[i]][global_dna_idx[j]] = step_ + 1  # tagged correspondence matrix

            sums = np.sum(global_correspondence, axis=1)
            idx_remove = np.argwhere(sums == 1)
            idx_remove = np.squeeze(idx_remove)
            rna_idx = np.delete(global_rna_idx,
                                idx_remove)  # the rna cells whose correspondence was found in the previous step are removed from the batch


        result_df = pd.DataFrame(data=global_correspondence, columns=dna_cells, index=rna_cells)
        tagged_df = pd.DataFrame(data=tagged_correspondence, columns=dna_cells, index=rna_cells)

        # change the format of the tagged correspondence matrix similar to the binary correspondence matrix
        result_rna_tagged = []
        result_dna_tagged = []
        result_tags = []
        rna_cells_tagged = list(tagged_df.index)
        dna_cells_tagged = list(tagged_df.columns)
        for idx_ in range(len(rna_cells_tagged)):
            tmp_rna = list(tagged_df.iloc[idx_, :])
            for step__ in range(n_iters):
                if (step__ + 1) in tmp_rna:
                    tmp_rna_dna_index = tmp_rna.index(step__ + 1)
                    tmp_rna_dna = dna_cells_tagged[tmp_rna_dna_index]
                    result_dna_tagged.append(tmp_rna_dna)
                    result_rna_tagged.append(rna_cells_tagged[idx_])
                    result_tags.append(step__ + 1)
        tmp_result_tagged = pd.DataFrame(list(zip(result_dna_tagged, result_rna_tagged, result_tags)),
                                         columns=["predicted_dna_cell", "rna_cell", "step"])

        return tmp_result_tagged

def func(rand_num, m_p, n_p, s_p, dna, rna, clonal_dfs, clus, names_dir):

    # ========= set the seeds =========
    np.random.seed(rand_num)

    # ========= draw new proportions =========
    new_props = np.random.multinomial(m_p, np.random.dirichlet(np.ones(n_p) * s_p))

    # ========= check the sum of proportions =========
    assert sum(new_props) == dna.shape[1], \
    'sum of the randomly drawn clonal proportions is not equal to the total number of DNA cells'

    # ========= aggregate the sampled names =========
    new_dna_names = pd.concat([clonal_dfs[i]['cell'].sample(n = prop, replace=True, random_state = rand_num) for i, prop in enumerate(new_props)], 
        ignore_index = True)

    new_dna = dna.loc[:,new_dna_names]

    model_dat = MaCroDNA(rna, new_dna)
    df_macrodna_res = model_dat.cell2cell_assignment()
    merged_df = df_macrodna_res.merge(clus, left_on = 'rna_cell', 
        right_on = 'cell', how = 'inner', validate = "many_to_one")
    merged_df.rename(columns = {'clone': 'rna_clone'}, inplace=True)
    merged_df = merged_df.merge(clus, left_on = 'predicted_dna_cell',
        right_on = 'cell', how='inner', validate = "many_to_one")
    merged_df.rename(columns = {'clone': 'dna_clone'}, inplace=True)
    if merged_df.shape[0] == 0:
        acc = nan
    else:
        acc = merged_df.query("dna_clone == rna_clone").shape[0]/merged_df.shape[0]
    print(f"accuracy on the sampled dna data set and original rna data: {acc}")

    pd.DataFrame({'cell_ID': new_dna.columns.tolist()}).to_csv(os.path.join(names_dir, f"resampled_dna_{rand_num}.csv"))

    return (new_props, acc)


if __name__ == "__main__":

    # ========= directory and option names ========= 

    dna_src_dir = "./all_genes_log/"
    rna_src_dir = "./all_genes_log/"

    biop_sample = ["crc04", "crc10", "crc11"]
    clustering_settings = ["agglomerative_log", "agglomerative_raw", "intNMF_log", "intNMF_raw"]

    mother_dir = "./all_props_exp_res/"
    clus_path = "./Clusters/"

    if not os.path.isdir(mother_dir):
        os.mkdir(mother_dir)

    # ========= loop over all options =========
    for d in biop_sample:
        for clustering_setting in clustering_settings:

            # ========= read the DNA and RNA data =========
            dna = pd.read_csv(os.path.join(dna_src_dir,d+"_dna.csv"), index_col=0)
            rna = pd.read_csv(os.path.join(rna_src_dir,d+"_rna.csv"), index_col=0)

            # ========= read the cluster info =========
            cluster_dir = clus_path+clustering_setting
            if clustering_setting == "agglomerative_log" or clustering_setting == "agglomerative_raw":
                clus = pd.read_csv(os.path.join(cluster_dir,"dna_cluster_"+d+"_4.csv"), index_col=0)
            elif clustering_setting == "intNMF_log" or clustering_setting == "intNMF_raw":
                filenames = [f for f in os.listdir(cluster_dir) if os.path.isfile(os.path.join(cluster_dir, f))]
                for filename in filenames:
                    if d in filename:
                        clus = pd.read_csv(os.path.join(cluster_dir, filename), index_col=0)

            # ========= check if the two sets of dna names are identical =========
            assert set(dna.columns.to_list()) == set(clus['cell'].values.tolist()),\
            "some cell names are missing in either copy number data or in the cluster assignment data frame"

            # ========= path handle =========
            tgt_dir = f"props_exp_res_{clustering_setting}_{d}/"
            tgt_dir = os.path.join(mother_dir, tgt_dir)
            names_dir = os.path.join(tgt_dir, 'names')

            if not os.path.isdir(tgt_dir):
                os.mkdir(tgt_dir)

            if not os.path.isdir(names_dir):
                os.mkdir(names_dir)

            # ========= group the cells by their clones =========
            clonal_ids, clonal_dfs = list(zip(*clus.groupby(['clone'])))

            # ========= create the columns of the stats data frame =========
            cols = ['dat_type', 'accuracy'] + [f"clone_{clone_id}" for clone_id in clonal_ids]
            dat_dict = {c : [] for c in cols}


            # ========= call macrodna on the original data =========
            model_orig_dat = MaCroDNA(rna, dna)
            df_macrodna_res = model_orig_dat.cell2cell_assignment()
            merged_df = df_macrodna_res.merge(clus, left_on = 'rna_cell', 
                right_on = 'cell', how = 'inner', validate = "many_to_one")
            merged_df.rename(columns = {'clone': 'rna_clone'}, inplace=True)
            merged_df = merged_df.merge(clus, left_on = 'predicted_dna_cell',
                right_on = 'cell', how='inner', validate = "many_to_one")
            merged_df.rename(columns = {'clone': 'dna_clone'}, inplace=True)
            if merged_df.shape[0] == 0:
                acc = nan
            else:
                acc = merged_df.query("dna_clone == rna_clone").shape[0]/merged_df.shape[0]
            print(f"accuracy on the original data set: {acc}")

            # ========= fill out the table ==========
            dat_dict['dat_type'].append('original')
            dat_dict['accuracy'].append(acc)
            for i, clone_id in enumerate(clonal_ids): dat_dict[f"clone_{clone_id}"].append(clonal_dfs[i].shape[0])


            # save the names of the cells into csv files in names_dir
            pd.DataFrame({'cell_ID': dna.columns.tolist()}).to_csv(os.path.join(names_dir, 'original_dna.csv'))
            pd.DataFrame({'cell_ID': rna.columns.tolist()}).to_csv(os.path.join(names_dir, 'original_rna.csv'))


            # ========= resampling on the dna cells =========
            s_p = 1.0
            n_p = len(clonal_dfs)
            m_p = dna.shape[1]

            # ========= number of replicates for the current biopsy =========
            N = 10000
            # ========= number of processes in the pool =========
            P = 20

            pool = multiprocessing.Pool(P)
            start_time = time.perf_counter()
            processes = [pool.apply_async(func,
                args = (k, m_p, n_p, s_p, dna, rna, clonal_dfs, clus, names_dir)) for k in range(N)]
            returns = [p.get() for p in processes]
            pool.close()
            pool.join()
            finish_time = time.perf_counter()
            print(f"{N} random replicates for patient {d} were processed in {finish_time-start_time} seconds")


            # returned = func(0, m_p, n_p, s_p, dna, rna, clonal_dfs, option, clus, names_dir)
            for returned in returns:
                # ========= fill out the table with the sampled data =========
                dat_dict['dat_type'].append('resampled')
                dat_dict['accuracy'].append(returned[1])
                for i, clone_id in enumerate(clonal_ids): dat_dict[f"clone_{clone_id}"].append(returned[0][i])

            # create the data frame of all results 
            df = pd.DataFrame.from_dict(dat_dict)
            df.to_csv(os.path.join(tgt_dir, 'dat.csv'))











