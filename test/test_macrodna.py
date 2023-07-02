from MaCroDNA import MaCroDNA
import pandas as pd

if __name__ == '__main__':
    DataFolder = "test/"
    rna_df = pd.read_csv(DataFolder+"rna_data.csv", index_col=0)
    dna_df = pd.read_csv(DataFolder + "dna_data.csv", index_col=0)
    dna_cluster = pd.read_csv(DataFolder + "dna_cluster.csv", index_col=0)
    run_model = MaCroDNA(rna_df, dna_df, dna_cluster)
    cell2cell = run_model.cell2cell_assignment()
    cell2clone = run_model.cell2clone_assignment()
    print("**********************")
    print("TEST FINISHED")
    print("**********************")