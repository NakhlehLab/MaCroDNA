import numpy as np
import pandas as pd
import os
import shutil

if __name__=="__main__":
	sample_list = ["PAT20_CARD", "PAT20_ESO", "PAT9_NDBE", "PAT14_NDBE", "PAT16_NDBE", "PAT6_LGD", "PAT19_LGD", "PAT6_HGD", "PAT14_HGD", "PAT20_HGD1", "PAT20_HGD2", "PAT16_EAC"]
	src_dir = "./data/Busslinger_data/scDNAseq_chromosome_count_tables/"
	tgt_dir = "./data/Busslinger_data/scDNAseq_filtered_cells_pseudocount_copynumber_log/"

	# biopsy = "PAT20_CARD"
	for biopsy in sample_list:

		df = pd.read_csv(os.path.join(src_dir, biopsy+"_count_table.csv"), index_col=0, header=0, sep=",", low_memory=False)
		df = df.reset_index()
		df.columns = ['reference_name', 'start', 'end'] + df.columns.values.tolist()[3:]
		df.drop(0, axis=0, inplace=True)
		df['reference_name'] = 'chr' + df['reference_name'].astype(str)
		df['start'] = df['start'].astype(int)
		df['end'] = df['end'].astype(int)

		df_ = df.drop(['reference_name', 'start', 'end'], axis=1, inplace=False)

		# print(f"number of cells before filtering based on less than 3000 transcripts {df_.shape[1]}")
		if df_.isnull().values.any()  or df_.isna().values.any():
			print(f"there are nans in the biopsy, replacing them with zeros")
			df_ = df_.fillna(0)
		# print(f"apply the threshold of 3000 on each cell\'s coverage")
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

