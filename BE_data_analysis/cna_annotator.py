import pandas as pd
import numpy as np 
import os
import shutil

def gene_info(x):
# Extract gene names, gene_type, gene_status and level
	g_name = list(filter(lambda x: 'gene_name' in x,  x.split(";")))[0].split("=")[1]
	g_type = list(filter(lambda x: 'gene_type' in x,  x.split(";")))[0].split("=")[1]
	g_status = list(filter(lambda x: 'gene_status' in x,  x.split(";")))[0].split("=")[1]
	g_level = int(list(filter(lambda x: 'level' in x,  x.split(";")))[0].split("=")[1])
	return (g_name, g_type, g_status, g_level)

if __name__=="__main__":
	sample_list = ["PAT20_CARD", "PAT20_ESO", "PAT9_NDBE", "PAT14_NDBE", "PAT16_NDBE", "PAT6_LGD", "PAT19_LGD", "PAT6_HGD", "PAT14_HGD", "PAT20_HGD1", "PAT20_HGD2", "PAT16_EAC"]
	src_dir = "./data/Busslinger_data/scDNAseq_filtered_cells_pseudocount_copynumber_log/"
	tgt_dir = "./data/Busslinger_data/annotated_scDNAseq_filtered_cells_pseudocount_copynumber_log/"
	gff_f = "./gene_table_gencode.v19.annotation.gff3.csv"

	df_genes = pd.read_csv(gff_f, index_col=0, header=0, sep=",")

	for biopsy in sample_list:
		# biopsy = "PAT20_ESO"
		df = pd.read_csv(os.path.join(src_dir, biopsy+"_filtered_normed_count_table.csv"), index_col=0, header=0, sep=",", low_memory=False)

		df_genes['start'] = df_genes['start'].astype(int)
		df_genes['end'] = df_genes['end'].astype(int)
		df['start'] = df['start'].astype(int)
		df['end'] = df['end'].astype(int)

		new_header = df.columns.values.tolist()[3:]
		new_indices = []
		vals = []

		# iterate over rows of the biopsy
		for index, row in df.iterrows():
			x = df_genes[ (row['reference_name']==df_genes['seqname']) & (row['start']<=df_genes['start']) & (row['end']>=df_genes['end']) ]
			new_indices = new_indices + x['gene_name'].values.tolist()
			for _ in range(len(x['gene_name'].values)):
				vals.append(row.values[3:])
		vals = np.array(vals)
		new_df = pd.DataFrame(vals, index=new_indices, columns=new_header)
		print(new_df.shape)
		# check the uniqueness of the genes across the biopsy 
		if len(new_indices)!=len(set(new_indices)):
			print(f"there are repeared genes across the indices of biopsy {biopsy}")
		else:
			print(f"None of the genes are repeated.")
		new_df.to_csv(os.path.join(tgt_dir, biopsy+"_annotated_filtered_normed_count_table.csv"))
		print(f"{biopsy} is done!")

