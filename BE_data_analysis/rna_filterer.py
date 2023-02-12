import numpy as np
import pandas as pd
import os
import shutil

if __name__=="__main__":
	sample_list = ["PAT20_CARD", "PAT20_ESO", "PAT9_NDBE", "PAT14_NDBE", "PAT16_NDBE", "PAT6_LGD", "PAT19_LGD", "PAT14_HGD", "PAT20_HGD1", "PAT6_HGD", "PAT16_EAC"]
	src_dir = "./data/Busslinger_data/scRNAseq_chromosome_count_tables/"
	tgt_dir = "./data/Busslinger_data/scRNAseq_filtered_cells_genes_pseudocount_rpm_log/"

	for biopsy in sample_list:

		# read the rna read count table of the current biopsy
		# the columns are the cell id's and the rows are gene expressions
		df = pd.read_csv(os.path.join(src_dir, biopsy+"_count_table.csv"), index_col=0, header=0, sep=",", low_memory=False)

		# check for the nan values in the table, if exists, replace them with 0's
		if df.isnull().values.any():
			print(f"there are nans in the biopsy, replacing them with zeros")
			df = df.fillna(0)
		# filter out the cells whose total number of transcripts is less than 3000
		df = df[df.columns[df.sum()>3000]]
		# transpose the data frame for the next filtering 
		df = df.T
		# keep genes that are expressed in at least one cell with at least three transcripts
		# remove the rest of the genes from the data frame
		df = df[df.columns[(df>=3).any()]]
		# transpose the df again so that the columns are the cells and rows are genes 
		df = df.T
		# add the pseudocount (1) to all the counts after filtering the cells and genes 
		df = df + 1
		# RPM normalization
		df = df.div(df.sum())
		df = df.mul(1e6)
		# perform log-normalization 
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
