import numpy as np 
import pandas as pd
import os
import shutil
from collections import Counter



if __name__=="__main__":
	sample_list = ["PAT20_CARD", "PAT20_ESO", "PAT9_NDBE", "PAT14_NDBE", "PAT16_NDBE", "PAT6_LGD", "PAT19_LGD", "PAT14_HGD", "PAT20_HGD1", "PAT20_HGD2", "PAT16_EAC"]
	src_dir = "./data/Busslinger_data/"
	tgt_dir = "./data/Busslinger_data/scRNAseq_chromosome_count_tables/"

	biopsy = "PAT20_CARD"

	biopsies_codes = []
	biopsies_names = []
	history_ = open(os.path.join(src_dir, "scRNAseq_read_count_all_samples_history.txt"), "r").readlines()
	for line in history_:
		if "/PAT" in line:
			biopsies_codes.append(line.strip().split(":")[0].split(" ")[-1].replace("(", "").replace(")", ""))
			biopsies_names.append(line.strip().split(":")[1].split("/")[-2])

	print(biopsies_names)
	print(biopsies_codes)

	for biopsy in sample_list:
		print(f"extracting the scRNA data of the biopsy {biopsy}")

		index_ = []
		header_ = []
		cell_idx = []
		cell_names = None
		dat = []
		codes_in_rna = []
		scRNA_file = open(os.path.join(src_dir, "scRNAseq_read_count_all_samples.csv"), "r")
		flag = False
		for line in scRNA_file:
			if not flag:
				x = line.strip().split("\t")
				cell_names = x[1:]
				codes_in_rna = [w[:2] for w in cell_names]
				code = biopsies_codes[biopsies_names.index(biopsy.replace("_", ""))]
				for i in range(len(cell_names)):
					if cell_names[i].startswith(code):
						cell_idx.append(i)
				# index_.append(x[0])
				header_ = [cell_names[j] for j in cell_idx]
				flag = True
			else:
				x = line.strip().split("\t")
				counts = x[1:]
				dat.append([counts[k] for k in cell_idx])
				index_.append(x[0])
		dat = np.array(dat)
		df = pd.DataFrame(dat, columns = header_, index = index_)
		if biopsy == "PAT20_HGD2":
			biopsy = "PAT6_HGD"
		df.to_csv(os.path.join(tgt_dir, biopsy+"_count_table.csv"))
		print(df.shape)
		gene_names = [q.split("__chr")[0] for q in index_]
		print(len(gene_names), len(set(gene_names)))
		c = Counter(gene_names)
		r_f = open(os.path.join(tgt_dir, biopsy+"_repeated_genes.txt"), "w")
		for k in c:
			if c[k]>1:
				r_f.write(k+"\n")
		r_f.close()
		scRNA_file.close()



	print(set(codes_in_rna))
	print(len(set(codes_in_rna)), len(biopsies_codes))
	print(set(codes_in_rna).difference(set(biopsies_codes)))
	print(set(biopsies_codes).difference(set(codes_in_rna)))


