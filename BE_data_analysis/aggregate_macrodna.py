import pickle 
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
	# biopsies = pd.read_csv("./data/biopsy_sample_list.csv",index_col=0)["biopsy"].tolist()
	biopsies = ["PAT20_CARD", "PAT20_ESO", "PAT9_NDBE", "PAT14_NDBE", "PAT16_NDBE", "PAT6_LGD", "PAT19_LGD", "PAT6_HGD", "PAT14_HGD", "PAT20_HGD1", "PAT16_EAC"]
	df_cancer = pd.read_csv("./data/tier1_gene.csv", index_col=0)
	cancer_genes = df_cancer.values.tolist()
	cancer_genes = [x[0] for x in cancer_genes]

	# biopsies = ["PAT14_HGD"]

	for biopsy in biopsies:
		dna = pd.read_csv(os.path.join(dna_src_dir,biopsy+"_filtered_normed_count_table.csv"),index_col=0, header=0)
		dna.drop(['reference_name', 'start', 'end'], axis=1, inplace=True)
		rna = pd.read_csv(os.path.join(rna_src_dir,biopsy+"_filtered_normed_count_table.csv"),index_col=0, header=0)
		# genes = list(set(dna.index.values.tolist()).intersection(set(rna.index.values.tolist())))
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
		# now use the lists of the pairs to extract the subset of the original dataframes 
		# make sure the list of genes is passed too
		# for caluclating the pairwise distance matrices, we need the set of all genes
		# cosmic_genes = list(set(cancer_genes).intersection(set(genes)))
		# print(f"number of COSMIC cancer genes: {len(cancer_genes)}")
		# print(f"number of COSMIC genes found among the genes in the biopsy {biopsy}: {len(cosmic_genes)}")
		# export two dataframes for dna, one containing all genes and the other with only COSMIC genes
		dna_subset = dna.loc[:, res_dnas]
		# # for the rna data, we only need the COSMIC genes
		# rna_subset = rna.loc[cosmic_genes, res_rnas]
		# rna_subset.columns = res_dnas
		# rna_subset = rna_subset.T

		rna_subset_all = rna.loc[:, res_rnas]
		rna_subset_all.columns = res_dnas
		rna_subset_all = rna_subset_all.T


		dna_subset.to_csv(os.path.join(tgt_dir, biopsy+"_dna.csv"))
		# rna_subset.to_csv(os.path.join(tgt_dir, biopsy+"_rna_rpm_paired.csv"))
		rna_subset_all.to_csv(os.path.join(tgt_dir, biopsy+"_rna.csv"))
		print("Done!")






















