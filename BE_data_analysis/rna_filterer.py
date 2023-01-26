import numpy as np
import pandas as pd
import os
import shutil

if __name__=="__main__":
	sample_list = ["PAT20_CARD", "PAT20_ESO", "PAT9_NDBE", "PAT14_NDBE", "PAT16_NDBE", "PAT6_LGD", "PAT19_LGD", "PAT14_HGD", "PAT20_HGD1", "PAT6_HGD", "PAT16_EAC"]
	src_dir = "./data/Busslinger_data/scRNAseq_chromosome_count_tables/"
	tgt_dir = "./data/Busslinger_data/scRNAseq_filtered_cells_genes_pseudocount_rpm_log/"

	# biopsy = "PAT20_CARD"
	for biopsy in sample_list:

		df = pd.read_csv(os.path.join(src_dir, biopsy+"_count_table.csv"), index_col=0, header=0, sep=",", low_memory=False)
		# print(f"number of cells before filtering based on less than 3000 transcripts {df.shape[1]}")
		if df.isnull().values.any():
			print(f"there are nans in the biopsy, replacing them with zeros")
			df = df.fillna(0)
		# print(f"apply the threshold of 3000 on each cell\'s coverage")
		df = df[df.columns[df.sum()>3000]]
		# print(f"number of cells after filtering out the low quality cells in {biopsy}: {df.shape[1]}")
		df = df.T
		# print(f"number of genes before filtering them based on their expression in {biopsy}: {df.shape[1]}")

		df = df[df.columns[(df>=3).any()]]
		# print(f"number of genes after filtering them based on their expression in {biopsy}: {df.shape[1]}")
		# transpose the df again so that the columns are the cells and rows are genes 
		df = df.T
		# add the pseudocount (1) to all the counts after filtering the cells and genes 
		df = df + 1
		# RPM normalization
		df = df.div(df.sum())
		df = df.mul(1e6)
		df = np.log1p(df)

		# remove the suffix __chr from the gene names 
		old_indices = df.index.values.tolist()
		new_indices = [x.split("__chr")[0] for x in old_indices]
		d_  =dict(zip(old_indices, new_indices))
		df = df.rename(index=d_)
		# save the filtered, normalized scRNAseq data into a new file
		df.to_csv(os.path.join(tgt_dir, biopsy+"_filtered_normed_count_table.csv"))
		if (df < 0).values.any():
			print("there are negative values in the normalized data frame")
		print(f"{biopsy} is done!")
