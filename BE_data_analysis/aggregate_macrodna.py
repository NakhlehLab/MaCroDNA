 import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
sns.set_style("darkgrid", {"axes.facecolor": ".9"})
import os
import shutil

if __name__=="__main__":

	dna_src_dir = "./data/Busslinger_data/scDNAseq_filtered_cells_pseudocount_copynumber_log/"
	rna_src_dir = "./data/Busslinger_data/scRNAseq_filtered_cells_genes_pseudocount_rpm_log/"
	res_dir = "./data/Busslinger_data/macrodna_res_log/"
	tgt_dir = "./data/Busslinger_data/macrodna_paired/"
	if os.path.isdir(tgt_dir):
		shutil.rmtree(tgt_dir)
	os.mkdir(tgt_dir)
	biopsies = ["PAT20_CARD", "PAT20_ESO", "PAT9_NDBE", "PAT14_NDBE", "PAT16_NDBE", "PAT6_LGD", "PAT19_LGD", "PAT6_HGD", "PAT14_HGD", "PAT20_HGD1", "PAT16_EAC"]

	for biopsy in biopsies:
		dna = pd.read_csv(os.path.join(dna_src_dir,biopsy+"_filtered_normed_count_table.csv"),index_col=0, header=0)
		dna.drop(['reference_name', 'start', 'end'], axis=1, inplace=True)
		rna = pd.read_csv(os.path.join(rna_src_dir,biopsy+"_filtered_normed_count_table.csv"),index_col=0, header=0)
		res = pd.read_csv(os.path.join(res_dir,biopsy+"_cell2cell_assignment_indexed.csv"), header=0)
		print(f"number of dna cells with matched rna in {biopsy}: {len(set(res['predict_cell'].values.tolist()))}")
		print(f'number of dna cells in the original data of {biopsy}: {len(dna.columns.values.tolist())}')
		print(f'number of rna cells in the original data of {biopsy}: {len(rna.columns.values.tolist())}')

		# create a list for each dna cell in the result and its best pair
		res_dnas = []
		res_rnas = []
		num_steps = res['step'].values.tolist()
		max_step = max(num_steps)
		dna_cell_names = list(set(res['predict_cell'].values.tolist()))
		for s in range(max_step):
			df_subset = res.loc[res['step']==s+1]
			if not df_subset.empty:
				for index, row in df_subset.iterrows():
					if not row['predict_cell'] in res_dnas:
						res_dnas.append(row['predict_cell'])
						res_rnas.append(row['cell'])
		


		dna_subset = dna.loc[:, res_dnas]

		rna_subset_all = rna.loc[:, res_rnas]
		rna_subset_all.columns = res_dnas
		rna_subset_all = rna_subset_all.T


		dna_subset.to_csv(os.path.join(tgt_dir, biopsy+"_dna.csv"))
		rna_subset_all.to_csv(os.path.join(tgt_dir, biopsy+"_rna.csv"))
		print("Done!")






















