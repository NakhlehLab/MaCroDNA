import numpy as np
import pandas as pd
import os
import shutil

if __name__=="__main__":
	# list of the biopsies 
	sample_list = ["PAT20_CARD", "PAT20_ESO", "PAT9_NDBE", "PAT14_NDBE", "PAT16_NDBE", "PAT6_LGD", "PAT19_LGD", "PAT6_HGD", "PAT14_HGD", "PAT20_HGD1", "PAT16_EAC"]
	# path to the directory containing the read count tables
	src_dir = "./data/Busslinger_data/scDNAseq_chromosome_count_tables/"
	# path to the destination directory where we store the filtered read count tables 
	tgt_dir = "./data/Busslinger_data/scDNAseq_filtered_cells_pseudocount_copynumber_log/"

	for biopsy in sample_list:

		# read the dna read count table of the current biopsy 
		df = pd.read_csv(os.path.join(src_dir, biopsy+"_count_table.csv"), index_col=0, header=0, sep=",", low_memory=False)
		# there are two rows with headers in the tables, we merge them in the following lines
		df = df.reset_index()
		df.columns = ['reference_name', 'start', 'end'] + df.columns.values.tolist()[3:]
		df.drop(0, axis=0, inplace=True)
		df['reference_name'] = 'chr' + df['reference_name'].astype(str)
		df['start'] = df['start'].astype(int)
		df['end'] = df['end'].astype(int)
		df_ = df.drop(['reference_name', 'start', 'end'], axis=1, inplace=False)

		# check if there are nan or null values in the data, if exists, replace them with 0's
		if df_.isnull().values.any()  or df_.isna().values.any():
			print(f"there are nans in the biopsy, replacing them with zeros")
			df_ = df_.fillna(0)
		# remove the cells with coverage less than 3000
		df_ = df_[df_.columns[df_.sum()>3000]]
		print(f"number of cells after filtering out the low quality cells in {biopsy}: {df_.shape[1]}")
		# add pseudocounts to all counts 
		df_ = df_ + 1
		# calculate the rough estimations of the absolute copy number values
		df_ = df_.div(df_.median())
		df_ = df_.mul(2)
		# take the logarithm of the absolute copy number estimations
		df_ = np.log1p(df_)

		if (df_ < 0).values.any():
			print("there are negative values in the normalized data frame")

		# insert the reference, start, and end columns from the original dataframe 
		df_.insert(0, "end", df['end'].values, True)
		df_.insert(0, "start", df['start'].values, True)
		df_.insert(0, "reference_name", df['reference_name'].values, True)
		df_.to_csv(os.path.join(tgt_dir, biopsy+"_filtered_normed_count_table.csv"))
		print(f"biopsy {biopsy} is done!")

