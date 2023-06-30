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

        return results

    def cell2cell_assignment(self):
        dna_cells = list(self.dna_df.columns)
        rna_cells = list(self.rna_df.columns)
        genes = set(self.dna_df.index.to_list()).intersection(set(self.rna_df.index.to_list()))
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

        return tmp_result, tmp_result_tagged


if __name__ == "__main__":

    biop_sample = ["PAT20_CARD", "PAT20_ESO", "PAT9_NDBE", "PAT14_NDBE", "PAT16_NDBE", "PAT6_LGD", "PAT19_LGD", "PAT6_HGD", "PAT14_HGD", "PAT20_HGD1", "PAT16_EAC"]
    dna_src_dir = "./annotated_scDNAseq_filtered_cells_pseudocount_copynumber_log"
    rna_src_dir = "./scRNAseq_filtered_cells_genes_pseudocount_rpm_log"
    tgt_dir = "./macrodna_res_log_cross_cells/"
    if not os.path.isdir(tgt_dir):
        os.mkdir(tgt_dir)

    genes_ = []
    for d_biop in biop_sample:
        for r_biop in biop_sample:
            dna = pd.read_csv(os.path.join(dna_src_dir,d_biop+"_annotated_filtered_normed_count_table.csv"),index_col=0)
            rna = pd.read_csv(os.path.join(rna_src_dir,r_biop+"_filtered_normed_count_table.csv"), index_col=0)
            print(f"dna {d_biop} {dna.shape}, rna {r_biop}, {rna.shape}")
            genes_.append(dna.index.to_list())
            genes_.append(rna.index.to_list())
    shared_genes = set.intersection(*map(set,genes_))
    print(f"total number of shared genes across all biopsies {len(shared_genes)}")
    # aggregate all dna cells into one data frame 
    dna_dat = []
    new_dna_cols = []
    for d_biop in biop_sample:
        dna = pd.read_csv(os.path.join(dna_src_dir,d_biop+"_annotated_filtered_normed_count_table.csv"),index_col=0)
        dna_dat.append(dna.loc[shared_genes, :].T)
        new_dna_cols += [d_biop.replace("_", "(")+")_"+col for col in dna.columns.to_list()]
    entire_dna_df = pd.concat(dna_dat, axis=0).T
    entire_dna_df.columns = new_dna_cols
    print(entire_dna_df.shape)

    # aggregate all rna cells into one data frame 
    rna_dat = []
    new_rna_cols = []
    for r_biop in biop_sample:
        rna = pd.read_csv(os.path.join(rna_src_dir,r_biop+"_filtered_normed_count_table.csv"), index_col=0)
        rna_dat.append(rna.loc[shared_genes, :].T)
        new_rna_cols += [r_biop.replace("_", "(")+")_"+col for col in rna.columns.to_list()]
    entire_rna_df = pd.concat(rna_dat, axis=0).T
    entire_rna_df.columns = new_rna_cols
    print(entire_rna_df.shape)


    # check the uniqueness of the cell ID's in both modalities
    assert len(entire_dna_df.columns.to_list()) == len(set(entire_dna_df.columns.to_list())), f'cell IDs in DNA data are not unique: number of ids {len(entire_dna_df.columns.to_list())}, number of unique ids {len(set(entire_dna_df.columns.to_list()))}'
    assert len(entire_rna_df.columns.to_list()) == len(set(entire_rna_df.columns.to_list())), f'cell IDs in RNA data are not unique: number of ids {len(entire_rna_df.columns.to_list())}, number of unique ids {len(set(entire_rna_df.columns.to_list()))}'
    print(f"the cells IDs in both modalities are unique!")


    # run MaCroDNA
    run_model = MaCroDNA(entire_rna_df, entire_dna_df)
    tmp_result, tmp_result_tagged = run_model.cell2cell_assignment()
    tmp_result_tagged.to_csv(os.path.join(tgt_dir, "cross_cells_cell2cell_assignment_indexed.csv"))
    tmp_result.to_csv(os.path.join(tgt_dir, "cross_cells_cell2cell.csv"))











