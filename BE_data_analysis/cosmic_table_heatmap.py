import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("darkgrid", {"axes.facecolor": ".87"})
# sns.set_style("darkgrid")
# sns.set_style("whitegrid")
import os

def argsort(seq):
	# http://stackoverflow.com/questions/3071415/efficient-method-to-calculate-the-rank-vector-of-a-list-in-python
	return sorted(range(len(seq)), key=seq.__getitem__)


if __name__=="__main__":
	# df_cancer = pd.read_csv("./data/tier1_gene.csv", index_col=0)
	df_cancer = pd.read_csv("./data/Census_all.csv",header=0)
	print(df_cancer)
	print(df_cancer.columns.tolist())
	cancer_genes = df_cancer['Gene Symbol'].values.tolist()
	print(cancer_genes)
	cancer_genes = [x for x in cancer_genes]
	sample_list = ["PAT20_CARD", "PAT20_ESO", "PAT9_NDBE", "PAT14_NDBE", "PAT16_NDBE", "PAT6_LGD", "PAT19_LGD", "PAT6_HGD", "PAT14_HGD", "PAT20_HGD1", "PAT16_EAC"]
	src_dir = f"./data/Busslinger_data/phylosignal_res_K.star"
	dna_src_dir = "./data/Busslinger_data/macrodna_paired/"
	q = []
	r = []
	cellvalues_ = []
	cellvalues_cosmic = []
	gene_set_arr = []
	gene_set_arr_names = []
	overlap_ = []
	max_genes = []
	max_indices = []
	n = len(sample_list)
	colors = sns.husl_palette(l=0.5, s=1.0, h=0.8, as_cmap=True)(np.linspace(0,1,n))
	for sample_name in sample_list:
		dna_df = pd.read_csv(os.path.join(dna_src_dir,sample_name+"_dna.csv"),index_col=0, header=0)
		N = len(dna_df.columns.values.tolist())
		pvals_df = pd.read_csv(os.path.join(src_dir, sample_name+"_pval.csv"), header=0, index_col=0, sep=",")
		stats_df = pd.read_csv(os.path.join(src_dir, sample_name+"_stat.csv"), header=0, index_col=0, sep=",")
		cols = stats_df.columns.values.tolist()
		s = stats_df[cols[0]].to_numpy()
		p = pvals_df[cols[0]].to_numpy()
		print(f"min of pvalues {np.amin(p)}")
		names = stats_df.index.values.tolist()
		cutoff = 0.8
		p_cutoff = 0.05
		indices = []
		tmp = []
		for i in range(8):
			idx = np.argwhere((s>cutoff+i*0.1) & (p<p_cutoff))
			selected_genes = [names[x] for x in idx[:,0]]
			selected_genes_bool = [1 if g in cancer_genes else 0 for g in selected_genes]

			tmp.append(sum(selected_genes_bool))

		idx = np.argwhere((s>1.0) & (p<p_cutoff))
		max_index = 0
		max_gene = "---"
		max_index_all = 0
		max_gene_all = "---"
		g = []
		gene_set_arr_names.append(sample_name)
		for x in idx[:,0]:
			if names[x] in cancer_genes:
				g.append(names[x])
				if s[x] > max_index:
					max_index = s[x]
					max_gene = names[x]
			if s[x] > max_index_all:
				max_index_all = s[x]
				max_gene_all = names[x]

		cosmic_genes_ = []
		cosmic_Ks = []
		if max_index!=0:
			f = open(f"{sample_name}_cosmic_genes_sorted.txt", "w")
			for x in idx[:,0]:
				if names[x] in cancer_genes:
					cosmic_genes_.append(names[x])
					cosmic_Ks.append(s[x])
			sorted_idx = list(reversed(argsort(cosmic_Ks)))
			cosmic_genes__ = [cosmic_genes_[i] for i in sorted_idx]
			cosmic_Ks_ = [cosmic_Ks[j] for j in sorted_idx]
			f.write(", ".join(cosmic_genes__))
			f.close()



		gene_set_arr.append(g)
		w = np.argwhere((s>1.0) & (p<p_cutoff))
		cellvalues_.append([len(w[:,0]), tmp[2], max_index, max_gene])
		cellvalues_cosmic.append(tmp)

	for i in range(len(gene_set_arr)):
		overlap_.append([])
		for j in range(len(gene_set_arr)):
			print(gene_set_arr_names[i], gene_set_arr_names[j])
			overlap_[i].append(len(set(gene_set_arr[i]).intersection(set(gene_set_arr[j]))))

	##########################################################################
	# plot the summary table
	collabels = [" > 0.8", "> 0.9", "> 1.0", "> 1.1", "> 1.2", "> 1.3", "> 1.4", "> 1.5"]
	rownames = [x.replace("_", " (")+")" for x in sample_list]
	df = pd.DataFrame(cellvalues_cosmic, index = rownames, columns=collabels)
	fig, ax = plt.subplots()
	cbar_kws = {"orientation":"vertical", 
           } # color bar keyword arguments
 
	g = sns.heatmap(df, cmap="magma", annot = True, linewidth = 1, fmt='g', cbar_kws=cbar_kws, ax=ax)
	ax.set_title('Table of COSMIC genes',
             fontweight ="bold")
	plt.tight_layout()
	plt.savefig("./cosmic_genes_table.pdf")
	# plt.show()

	plt.clf()
	plt.cla()

	collabels = ["Genes", "COSMIC genes", "highest index", "gene with highest index"]
	fig, ax = plt.subplots()
	ax.set_axis_off()
	table = ax.table(
	    cellText = cellvalues_, 
	    rowLabels = rownames, 
	    colLabels = collabels,
	    rowColours =["skyblue"] * len(rownames),
	    colColours = ["skyblue"] * len(collabels),
	    cellLoc ='center', 
	    loc ='upper left')
	# table.set_fontsize(20)
	ax.set_title('Table of COSMIC genes',
             fontweight ="bold")
	table.auto_set_font_size(False)
	table.set_fontsize(7)
	plt.tight_layout()
	plt.savefig("./cosmic_genes_table2.pdf")

	plt.clf()
	plt.cla()
	collabels = gene_set_arr_names
	rownames = gene_set_arr_names
	fig, ax = plt.subplots()
	ax.set_axis_off()
	table = ax.table(
	    cellText = overlap_, 
	    rowLabels = rownames, 
	    colLabels = collabels,
	    rowColours =["skyblue"] * len(rownames),
	    colColours = ["skyblue"] * len(collabels),
	    cellLoc ='center', 
	    loc ='upper left')
	# table.set_fontsize(20)
	ax.set_title('Overlaps between COSMIC genes',
             fontweight ="bold")
	table.auto_set_font_size(False)
	table.set_fontsize(7)
	plt.tight_layout()
	plt.savefig("./cosmic_genes_overlaps.pdf")





