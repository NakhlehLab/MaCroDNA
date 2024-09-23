import numpy as np
import pandas as pd
from scipy import spatial
from numpy.linalg import norm
import math


class MaCroDNA:
    def __init__(self, rna_df=None, dna_df=None, dna_label=None, solver="gurobi"):
        # INPUT
        # rna_df, dna_df -- dataframe, cols are cell ids, index are genes
        # dna_label -- dataframe, two cols, one is the clone id "clone", the other is the cell id "cell"
        self.rna_df = rna_df
        self.dna_df = dna_df
        self.dna_label = dna_label

        if solver not in ["gurobi", "glpk"]:
            raise ValueError("Solver must be either 'gurobi' or 'glpk'")
        

        if solver == "gurobi":
            self.ilp = self.ilp_gurobi
    
        else:
            self.ilp = self.ilp_pyomo

            



    def cosine_similarity(self, x1, x2):
        return 1 - spatial.distance.cosine(x1 - x1.mean(), x2 - x2.mean())

    def cosine_similarity_np(self, x1, x2):
        return np.dot(x1 - x1.mean(), x2 - x2.mean()) / (1e-10 + norm(x1 - x1.mean()) * norm(x2 - x2.mean()))


    def ilp_pyomo(self, rna_idx, dna_idx, corrs):
        from pyomo.environ import ConcreteModel, Var, Objective
        from pyomo.environ import Constraint, SolverFactory, Binary, maximize, RangeSet, value

    
        model = ConcreteModel()


        num_rna = rna_idx.shape[0]
        num_dna = dna_idx.shape[0]
        n_min = min(num_rna, num_dna)

        print("The smallest set has %s number of cells" % (n_min))

        model.rna = RangeSet(0, num_rna - 1)
        model.dna = RangeSet(0, num_dna - 1)

        # Binary variables x[i, j] (correspondence between RNA and DNA)
        model.x = Var(model.rna, model.dna, domain=Binary)



        # Each RNA cell is assigned to at most one DNA cell (row constraints)
        def row_constraint(model, i):
            return sum(model.x[i, j] for j in model.dna) <= 1
        model.row_constraints = Constraint(model.rna, rule=row_constraint)

        # Each DNA cell is assigned to at most one RNA cell (column constraints)
        def col_constraint(model, j):
            return sum(model.x[i, j] for i in model.rna) <= 1
        model.col_constraints = Constraint(model.dna, rule=col_constraint)

        # Total number of assignments should be equal to n_min
        def global_constraint(model):
            return sum(model.x[i, j] for i in model.rna for j in model.dna) == n_min
        model.global_constraint = Constraint(rule=global_constraint)

        # Objective function (maximize the correlation between RNA and DNA assignments)
        def objective_function(model):
            return sum(model.x[i, j] * corrs[rna_idx[i]][dna_idx[j]] for i in model.rna for j in model.dna)
        model.obj = Objective(rule=objective_function, sense=maximize)

        #Solve the model
        solver = SolverFactory('glpk')  # or use another solver like 'cbc' or 'gurobi'
        results = solver.solve(model, tee=False, keepfiles=False, options={'output': '/dev/null'})
        # Extract the objective value
        objective_value = value(model.obj)

        print(f"The optimal objective value is: {objective_value}")

        # Extract the solution
        solution = np.empty((num_rna, num_dna))
        for i in model.rna:
            for j in model.dna:
                solution[i][j] = int(model.x[i, j].value)

        return solution
    

    def ilp_gurobi(self, rna_idx, dna_idx, corrs):
        import gurobipy as gp


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
                x[i].append(m.addVar(vtype=gp.GRB.BINARY, name="x[%d,%d]" % (i, j)))

        # set the contraints for rows and columns

        for i in range(rna_idx.shape[0]):
            # At most one cell per row
            m.addConstr(gp.quicksum([x[i][j] for j in range(dna_idx.shape[0])]) <= 1, name="row" + str(i))

        for j in range(dna_idx.shape[0]):
            # At most one cell per column
            m.addConstr(gp.quicksum([x[i][j] for i in range(rna_idx.shape[0])]) <= 1, name="col" + str(i))

        m.addConstr(gp.quicksum([x[i][j] for i in range(rna_idx.shape[0]) for j in range(dna_idx.shape[0])]) == n_min,
                    name="global")

        # update the model
        m.update()

        # create the linear expression
        obj = gp.LinExpr()

        for i in range(rna_idx.shape[0]):
            for j in range(dna_idx.shape[0]):
                obj += x[i][j] * corrs[rna_idx[i]][dna_idx[j]]

        # set the objective
        m.setObjective(obj, gp.GRB.MAXIMIZE)

        # optimize the model
        m.optimize()

        print('IsMIP: %d' % m.IsMIP)
        if m.status == gp.GRB.Status.INFEASIBLE:
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
        genes = list(set(self.dna_df.index).intersection(self.rna_df.index.to_list()))
        self.dna_df = self.dna_df.loc[genes,:]
        self.rna_df = self.rna_df.loc[genes,:]

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

    def cell2clone_assignment(self):

        rna_result, _ = self.cell2cell_assignment()
        dna_cells = list(self.dna_df.columns)

        # map to clone
        dna_label = self.dna_label
        dna_label = dna_label.set_index("cell")

        dna_result_label = dna_label.loc[rna_result["predict_cell"]]["clone"].tolist()
        rna_result["predict"] = dna_result_label
        return rna_result

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
        res = self.cell2clone_assignment()

        print("**********")
        print("Finish Mapping")
        print("Test Success")
        print("**********")


if __name__ == '__main__':
    test = MaCroDNA()
    test.tiny_test()