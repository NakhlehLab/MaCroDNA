import numpy as np
import gurobipy as gp
from gurobipy import *
import pandas as pd
from scipy import spatial
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

        rna_genes = list(self.rna_df.index)
        dna_genes = list(self.dna_df.index)
        same_genes = set(rna_genes).intersection(dna_genes)

        self.dna_df = self.dna_df.loc[same_genes, :]
        self.rna_df = self.rna_df.loc[same_genes, :]


        dna_np = self.dna_df.T.to_numpy()
        rna_np = self.rna_df.T.to_numpy()
        print("After selecting the same genes in RNA and DNA")
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
        print(corrs)

        # create a matrix for storing the global correspondences
        global_correspondence = np.zeros((rna_np.shape[0], dna_np.shape[0]))
        # create index lists for dna and rna
        global_dna_idx = np.array([i for i in range(dna_np.shape[0])])
        global_rna_idx = np.array([i for i in range(rna_np.shape[0])])

        # calculate the number of steps
        quotient, remainder = divmod(max(global_rna_idx.shape[0], global_dna_idx.shape[0]),
                                     min(global_rna_idx.shape[0], global_dna_idx.shape[0]))
        print(quotient, remainder)
        n_iters = int(math.ceil(quotient))
        print("MaCroDNA will be run for %s steps" % (n_iters))

        # iterations
        rna_idx = np.copy(global_rna_idx)

        for _ in range(n_iters):

            r = self.ilp(rna_idx, global_dna_idx, corrs)
            # identify the indices with assignments
            r = np.array(r)
            for i in range(r.shape[0]):
                for j in range(r.shape[1]):
                    if r[i][j] == 1:
                        global_correspondence[rna_idx[i]][global_dna_idx[j]] = 1

            sums = np.sum(global_correspondence, axis=1)
            idx_remove = np.argwhere(sums == 1)
            idx_remove = np.squeeze(idx_remove)
            rna_idx = np.delete(global_rna_idx, idx_remove)

        print("the number of associations in the correspondence matrix %s" % np.sum(global_correspondence))

        result_df = pd.DataFrame(data=global_correspondence, columns=dna_cells, index=rna_cells)

        # find the dna that map to rna
        result_rna = []
        result_dna = []
        rna_cells = list(result_df.index)
        dna_cells = list(result_df.columns)
        for d in rna_cells:
            tmp_rna = list(result_df.loc[d, :])
            tmp_rna_dna_index = tmp_rna.index(1)
            tmp_rna_dna = dna_cells[tmp_rna_dna_index]
            result_dna.append(tmp_rna_dna)
            result_rna.append(d)
        tmp_result = pd.DataFrame(list(zip(result_dna, result_rna)), columns=["predict_cell", "cell"])
        tmp_result = tmp_result.set_index("cell")

        return tmp_result

    def cell2clone_assignment(self):

        #     get cells in both rna and dna
        rna_result = self.cell2cell_assignment()
        dna_cells = list(self.dna_df.columns)
        same_cell = list(set(list(rna_result.index)).intersection(dna_cells))
        result_same = rna_result.loc[same_cell, :]

        # map to clone
        dna_label = self.dna_label
        dna_label = dna_label.set_index("cell")

        dna_result_label = dna_label.loc[result_same["predict_cell"]]["clone"].tolist()
        result_same["predict_clone"] = dna_result_label
        return result_same

    def tiny_test(self):
        dna_data = pd.DataFrame.from_dict({"cell1": [2, 2, 3, 1, 6, 2],
                                           "cell2": [2, 2, 2, 2, 2, 2],
                                           "cell3": [1, 1, 2, 2, 2, 3],
                                           "cell4": [2, 2, 2, 2, 2, 6],
                                           "gene": ["g1", "g2", "g3", "g4", "g5", "g6"]})
        dna_data = dna_data.set_index("gene")

        rna_data = pd.DataFrame.from_dict({"cell1": [0, 0, 10, 0, 20, 0, 0],
                                           "cell2": [2, 2, 2, 2, 2, 2, 0],
                                           "cell3": [0, 0, 2, 2, 0, 5, 0],
                                           "cell4": [1, 1, 1, 1, 1, 20, 0],
                                           "gene": ["g1", "g2", "g3", "g4", "g5", "g6", "g7"]})
        rna_data = rna_data.set_index("gene")

        dna_cluster = pd.DataFrame.from_dict({"clone":[0, 1, 2, 3],
                                              "cell": ["cell1", "cell2","cell3", "cell4"]})

        print("******Test DNA data is:")
        print(dna_data)
        print("******Test RNA data is:")
        print(rna_data)
        print("******Clone id for each DNA cell is:")
        print(dna_cluster)
        print("**********")
        print("Start Mapping RNA cells to DNA clones")
        print("**********")

        self.dna_df = dna_data
        self.rna_df = rna_data
        self.dna_label = dna_cluster
        self.cell2clone_assignment()

        print("**********")
        print("Finish Mapping")
        print("Test Success")
        print("**********")


if __name__ == '__main__':
    test = MaCroDNA()
    test.tiny_test()