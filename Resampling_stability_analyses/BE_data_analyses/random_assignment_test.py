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
import random
import pickle
import multiprocessing
import time



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
        print("the smallest set has %s number of cells" % (n_min))
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

        return results, m.objVal

    def cell2cell_assignment(self):
        dna_cells = list(self.dna_df.columns)
        rna_cells = list(self.rna_df.columns)
        genes = set(self.dna_df.index).intersection(self.rna_df.index.to_list())
        self.dna_df = self.dna_df.loc[genes, :]
        self.rna_df = self.rna_df.loc[genes, :]

        dna_np = self.dna_df.T.to_numpy()
        rna_np = self.rna_df.T.to_numpy()
        print("number of cells in dna data %s" % (dna_np.shape[0]))
        print("number of cells in rna data %s" % (rna_np.shape[0]))
        print("number of genes in dna data %s" % (dna_np.shape[1]))
        print("number of genes in rna data %s" % (rna_np.shape[1]))

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
        print(quotient, remainder)
        n_iters = int(quotient)
        if remainder != 0:
            n_iters += 1
        print("MaCroDNA will be run for %s steps" % (n_iters))

        rna_idx = np.copy(global_rna_idx)

        obj_scores = []

        # iterations

        for step_ in range(n_iters):

            r, obj_s = self.ilp(rna_idx, global_dna_idx, corrs)
            obj_scores.append(obj_s)
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

        print("the number of associations in the correspondence matrix %s" % np.sum(global_correspondence))

        result_df = pd.DataFrame(data=global_correspondence, columns=dna_cells, index=rna_cells)
        tagged_df = pd.DataFrame(data=tagged_correspondence, columns=dna_cells, index=rna_cells)

        # find the dna that map to rna
        result_rna = []
        result_dna = []
        rna_cells = list(result_df.index)
        dna_cells = list(result_df.columns)
        test = result_df.sum(axis=1)
        for d in rna_cells:
            tmp_rna = list(result_df.loc[d, :])
            tmp_rna_dna_index = tmp_rna.index(1)
            tmp_rna_dna = dna_cells[tmp_rna_dna_index]
            result_dna.append(tmp_rna_dna)
            result_rna.append(d)
        tmp_result = pd.DataFrame(list(zip(result_dna, result_rna)), columns=["predict_cell", "cell"])
        tmp_result = tmp_result.set_index("cell")

        # change the format of the tagged correspondence matrix similar to the binary correspondence matrix
        result_rna_tagged = []
        result_dna_tagged = []
        result_tags = []
        rna_cells_tagged = list(tagged_df.index)
        dna_cells_tagged = list(tagged_df.columns)
        for d in rna_cells_tagged:
            tmp_rna = list(tagged_df.loc[d, :])
            for step__ in range(n_iters):
                if (step__ + 1) in tmp_rna:
                    tmp_rna_dna_index = tmp_rna.index(step__ + 1)
                    tmp_rna_dna = dna_cells_tagged[tmp_rna_dna_index]
                    result_dna_tagged.append(tmp_rna_dna)
                    result_rna_tagged.append(d)
                    result_tags.append(step__ + 1)
        tmp_result_tagged = pd.DataFrame(list(zip(result_dna_tagged, result_rna_tagged, result_tags)),
                                         columns=["predict_cell", "cell", "step"])
        tmp_result_tagged = tmp_result_tagged.set_index("cell")
        print(f"objective values from all steps {obj_scores}")
        return tmp_result, tmp_result_tagged, sum(obj_scores)

class random_test:
    def __init__(self, rna_df=None, dna_df=None):
        self.rna_df = rna_df
        self.dna_df = dna_df
        self.dna_cells = list(self.dna_df.columns)
        self.rna_cells = list(self.rna_df.columns)

        self.genes = set(self.dna_df.index).intersection(self.rna_df.index.to_list())
        self.dna_df = self.dna_df.loc[self.genes, :]
        self.rna_df = self.rna_df.loc[self.genes, :]


        self.dna_np = self.dna_df.T.to_numpy()
        self.rna_np = self.rna_df.T.to_numpy()
        # print("number of cells in dna data %s" % (self.dna_np.shape[0]))
        # print("number of cells in rna data %s" % (self.rna_np.shape[0]))
        # print("number of genes in dna data %s" % (self.dna_np.shape[1]))
        # print("number of genes in rna data %s" % (self.rna_np.shape[1]))
        # the rows of the matrix correspond the rna cells and the columns correspond the dna cells
        self.corrs = np.empty((self.rna_np.shape[0], self.dna_np.shape[0]))

        for i in range(self.rna_np.shape[0]):
            for j in range(self.dna_np.shape[0]):
                self.corrs[i][j] = self.cosine_similarity_np(self.rna_np[i], self.dna_np[j])



        # calculate the number of steps
        self.quotient, self.remainder = divmod(self.rna_np.shape[0],
                                     self.dna_np.shape[0])
        # print(self.quotient, self.remainder)
        self.n_iters = int(self.quotient)
        if self.remainder != 0:
            self.n_iters += 1

        # print(self.quotient, self.remainder)
        # print(f"number of steps {self.n_iters}")


    def cosine_similarity_np(self, x1, x2):
        return np.dot(x1 - x1.mean(), x2 - x2.mean()) / (1e-10 + norm(x1 - x1.mean()) * norm(x2 - x2.mean()))


    def assign(self):

        global_dna_idx = np.array([i for i in range(self.dna_np.shape[0])])
        global_rna_idx = np.array([i for i in range(self.rna_np.shape[0])])

        rna_idx = np.copy(global_rna_idx)
        dna_idx = np.copy(global_dna_idx)

        s = 0
        to_remove = []

        for step_ in range(self.n_iters):
            # print([len(rna_idx), len(dna_idx)])
            n_min = min([len(rna_idx), len(dna_idx)])
            # print(n_min)
            if n_min == len(rna_idx):
                selected_idx = np.random.choice(dna_idx, n_min, replace=False)
                selected_corr = self.corrs[rna_idx, selected_idx]
                # print(selected_corr.shape)
                # print(np.sum(selected_corr))
                s += np.sum(selected_corr)
            elif n_min == len(dna_idx):
                selected_idx = np.random.choice(rna_idx, n_min, replace=False)
                selected_corr = self.corrs[selected_idx, dna_idx]
                s += np.sum(selected_corr)
                to_remove.append(selected_idx)
                rna_idx = np.delete(global_rna_idx, np.concatenate(to_remove))
                # print(len(rna_idx))
            # selected_idx = np.random.choice(rna_idx, n_min, replace=False)

        return s

def func(dna_src_dir, rna_src_dir, tgt_dir, d, N = int(1e6)):

    dna = pd.read_csv(os.path.join(dna_src_dir,d+"_annotated_filtered_normed_count_table.csv"),index_col=0)
    rna = pd.read_csv(os.path.join(rna_src_dir,d+"_filtered_normed_count_table.csv"), index_col=0)

    # radnom assignments
    random_assigner = random_test(rna, dna)
    random_scores = []
    max_val = 0
    for i in range(N):
        sampled_s = random_assigner.assign()
        # print(f"sum of pearson correlation values of pairs, random sample number {i}: {sampled_s}")
        if sampled_s > max_val:
            max_val = sampled_s
        random_scores.append(sampled_s)
        # if i % int(1e5) == 0:
        #     print(f"maximum value of all samples until iteration {i}, {max_val}")
        # random_scores.append(random_assigner.assign())
    # print(f"maximum value of the all samples {max_val}")


    # with open(os.path.join(tgt_dir, d+"_random_samples.pkl"), 'wb') as filehandle:
    #     pickle.dump(random_scores, filehandle)

    random_scores = np.array(random_scores)
    np.save(os.path.join(tgt_dir, d+"_random_samples"), random_scores)


    return max_val


if __name__ == "__main__":

    biop_sample = ["PAT20_CARD", "PAT20_ESO", "PAT9_NDBE", "PAT14_NDBE", "PAT16_NDBE", "PAT6_LGD", "PAT19_LGD", "PAT6_HGD", "PAT14_HGD", "PAT20_HGD1", "PAT16_EAC"]
    # biop_sample = ["PAT20_ESO"]
    dna_src_dir = "./annotated_scDNAseq_filtered_cells_pseudocount_copynumber_log"
    rna_src_dir = "./scRNAseq_filtered_cells_genes_pseudocount_rpm_log"
    tgt_dir = "./macrodna_res_log_random_test/"
    if not os.path.isdir(tgt_dir):
        os.mkdir(tgt_dir)
    # set the seed
    seed = 2023
    random.seed(seed)
    np.random.seed(seed)

    dict_ = {}

    for d in biop_sample:
        dna = pd.read_csv(os.path.join(dna_src_dir,d+"_annotated_filtered_normed_count_table.csv"),index_col=0)
        rna = pd.read_csv(os.path.join(rna_src_dir,d+"_filtered_normed_count_table.csv"), index_col=0)
        # print(f"DNA cells: {dna.shape[1]}, RNA cells: {rna.shape[1]}, num of genes: {dna.shape[0]}")
        print(f"name of the biopsy: {d}")
        print(f"-------------------------------------")
        run_model = MaCroDNA(rna, dna)
        tmp_result, tmp_result_tagged, obj_val = run_model.cell2cell_assignment()
        print(f"sum of all preason correlation values of the assigned pairs {obj_val}")
        dict_[d] = obj_val
        tmp_result_tagged.to_csv(os.path.join(tgt_dir, d+"_cell2cell_assignment_indexed.csv"))
        # tmp_result.to_csv(os.path.join(tgt_dir, d+"_cell2cell.csv"))

    with open(os.path.join(tgt_dir, "macrodna_objvals.pkl"), 'wb') as file_:
        pickle.dump(dict_, file_)

    N = int(1e8)
    print(f"start generating random assignments for all biopsies")
    pool = multiprocessing.Pool(len(biop_sample))
    start_time = time.perf_counter()
    processes = [pool.apply_async(func, args=(dna_src_dir, rna_src_dir, tgt_dir, d, N)) for d in biop_sample]
    max_vals = [p.get() for p in processes]
    pool.close()
    pool.join()
    finish_time = time.perf_counter()
    print(f"random assignment was done for all biopsies in {finish_time-start_time} seconds")
    print(max_vals)











