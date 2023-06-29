import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import seaborn as sns
import os

import numpy as np
from scipy import stats
import pandas as pd
from scipy import spatial
from scipy.stats import spearmanr
from numpy.linalg import norm
from scipy.spatial.distance import pdist
import math
import random
import pickle
import multiprocessing
import time
import heapq
import sklearn
from sklearn.linear_model import LinearRegression
from matplotlib.cbook import boxplot_stats
sns.set_style("whitegrid")
sns.set_context('talk')
sns.despine(left=False)

def obj_vals_hists(rna_cells, sample, src_dir, original_val):
	# read the LOO files corresponding to each RNA cells
	loo_obj_vals = []
	for cell_idx in range(len(rna_cells)):
		df = pd.read_csv(os.path.join(src_dir, f"{sample}_loo_RNA_idx_{cell_idx}.csv"))
		loo_obj_vals.append(np.sum(df['corr_val'].values))

	# plot the histogram of the objective values
	fig, ax = plt.subplots()
	ax.hist(loo_obj_vals, bins=50, density=False)
	ax.axvline(original_val, ls="--", color="r")
	_ = ax.set_ylabel("Counts")
	ax.set_title(sample.replace("_", " (").replace("PAT", "Patient ")+")")
	plt.savefig(os.path.join("./", sample)+"_loo_obj_vals_hist.pdf"
		, transparent=True)

	return 

def E_loo_error(df_entire, src_dir):
	# calculate the original value of the objective function
	original_val = np.sum(df_entire['corr_val'].values)
	# print(original_val)
	D_vals = df_entire['corr_val'].values
	D_ids = df_entire['rna_cell'].values
	f_D_dict = dict(zip(D_ids, D_vals))

	f_Di_dict = {}
	for c in D_ids:
		f_Di_dict[c] = []

	for cell_idx in range(len(rna_cells)):
		df = pd.read_csv(os.path.join(src_dir, f"{sample}_loo_RNA_idx_{cell_idx}.csv"))
		current_cell_ids = df['rna_cell'].values
		current_cell_vals = df['corr_val'].values

		for idx in range(len(current_cell_ids)):
			f_Di_dict[current_cell_ids[idx]].append(current_cell_vals[idx])
			
	abs_dict = {}
	a = []
	for c in D_ids:
		abs_dict[c] = [abs(x - f_D_dict[c]) for x in f_Di_dict[c]]
		a += [abs(x - f_D_dict[c]) for x in f_Di_dict[c]]
	print(len(a))
	total_std = np.std(np.array(a))
	total_avg = np.mean(np.array(a))
	print(f"average LOO error in sample {sample}: {total_avg}")
	print(f"std of LOO error in sample {sample}: {total_std}")

	return np.std(np.array(a))

def calculate_pdist(dna_src_dir, sample):

	dna = pd.read_csv(os.path.join(dna_src_dir,sample+"_annotated_filtered_normed_count_table.csv"),index_col=0)
	dna_np = dna.T.to_numpy()
	# pdv = pdist(dna_np, 'euclidean')
	# pdv = pdist(dna_np, 'cityblock')
	pdv = pdist(dna_np, 'minkowski', p=1.)
	return pdv

def label_point(x, y, val, ax):
	a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
	for i, point in a.iterrows():
		ax.text(point['x']-.02, point['y']-0.1, str(point['val']))



if __name__ == "__main__":

	biop_sample = ["PAT20_CARD", "PAT20_ESO", "PAT9_NDBE", "PAT14_NDBE", "PAT16_NDBE", "PAT6_LGD", "PAT19_LGD", "PAT6_HGD", "PAT14_HGD", "PAT20_HGD1", "PAT16_EAC"]
	# biop_sample = ["PAT20_ESO"]
	src_dir = "./macrodna_res_log_rna_stability/"
	rna_src_dir = "/Users/edrisi/Documents/ongoing_projects/MaCroDNA/phylosignal_analysis/Busslinger_data/scRNAseq_filtered_cells_genes_pseudocount_rpm_log"
	dna_src_dir = "/Users/edrisi/Documents/ongoing_projects/MaCroDNA/phylosignal_analysis/Busslinger_data/annotated_scDNAseq_filtered_cells_pseudocount_copynumber_log"

	x = []
	y = []
	z = []
	df_diversity_idx = {'sample': [], 'diversity_idx': []}
	for sample in biop_sample:
		# read the rna data
		rna = pd.read_csv(os.path.join(rna_src_dir, sample+"_filtered_normed_count_table.csv"), index_col=0)
		rna_cells = rna.columns.tolist()

		dna = pd.read_csv(os.path.join(dna_src_dir,sample+"_annotated_filtered_normed_count_table.csv"),index_col=0)
		dna_cells = dna.columns.tolist()

		# read the result of MaCroDNA on the entire data set 
		df_entire = pd.read_csv(os.path.join(src_dir, sample+"_cell2cell_assignment_indexed.csv"))

		original_assignments = dict(zip(df_entire['rna_cell'], df_entire['predicted_dna_cell']))
		# print(original_assignments)

		all_pairs_counts = np.zeros((len(rna_cells), len(dna_cells)))


		loo_assignments = {}
		loo_different_assign_counts = {}
		loo_different_assign_ratio = {}
		loo_different_assign_types = {}
		loo_different_assign_types_counts = {}
		for c in rna_cells:
			loo_assignments[c] = []
			loo_different_assign_counts[c] = 0
			# loo_different_assign_types_counts[c] = 0
			loo_different_assign_types[c] = set()

		for cell_idx in range(len(rna_cells)):
			df = pd.read_csv(os.path.join(src_dir, f"{sample}_loo_RNA_idx_{cell_idx}.csv"))
			current_cell_ids = df['rna_cell']
			current_cell_assignments = df['predicted_dna_cell']
			for i in range(len(current_cell_ids)):
				loo_assignments[current_cell_ids[i]].append(current_cell_assignments[i])
				if current_cell_assignments[i] != original_assignments[current_cell_ids[i]]:
					loo_different_assign_counts[current_cell_ids[i]] += 1
					loo_different_assign_types[current_cell_ids[i]].add(current_cell_assignments[i])


				rna_cell_idx = rna_cells.index(current_cell_ids[i])
				dna_cell_idx = dna_cells.index(current_cell_assignments[i])
				all_pairs_counts[rna_cell_idx][dna_cell_idx] += 1

		for i in range(len(all_pairs_counts)):
			all_pairs_counts[i] /= all_pairs_counts.shape[1]
		# g = sns.heatmap(all_pairs_counts, cmap="Blues")
		# plt.show()
		for key in loo_different_assign_counts:
			# loo_different_assign_ratio[key] = loo_different_assign_counts[key]/len(rna_cells)
			loo_different_assign_ratio[key] = loo_different_assign_counts[key]
			# loo_different_assign_types_counts[key] = len(loo_different_assign_types[key])/len(dna_cells)
			# if len(loo_different_assign_types[key]) != 0:
			loo_different_assign_types_counts[key] = len(loo_different_assign_types[key])


		pdv = calculate_pdist(dna_src_dir = dna_src_dir, sample = sample)
		print(f"median of the pairwise distances in sample {sample}: {np.median(pdv)}")
		print(f"max number of types of different assignments in sample {sample}: {np.mean(np.array(heapq.nlargest(10, list(loo_different_assign_types_counts.values()))))}")

		cells = []
		assigns = []
		for key in loo_different_assign_types_counts:
			cells.append(key)
			assigns.append(loo_different_assign_types_counts[key])


		# x.append(np.mean(pdv) + np.std(pdv))
		x.append(np.median(pdv))

		current_div_idx = list(loo_different_assign_types_counts.values())
		data = np.array(current_div_idx)
		outliers_boxplot = [y for stat in boxplot_stats(data) for y in stat['fliers']]
		print("outliers from the box plot are calculated")
		print(outliers_boxplot)
		q1, q3 = np.percentile(data, [25, 75])
		iqr = q3 - q1
		lower_bound = q1 - 1.5*iqr
		upper_bound = q3 + 1.5*iqr
		# outliers = data[(data < lower_bound) | (data > upper_bound)]
		outliers = data[(data > upper_bound)]
		print(outliers, upper_bound, lower_bound)
		# y.append(np.mean(np.array(heapq.nlargest(10, list(loo_different_assign_types_counts.values())))))
		y.append(np.mean(outliers))

		for val_idx in range(len(current_div_idx)):
			df_diversity_idx['sample'].append(sample)
			df_diversity_idx['diversity_idx'].append(current_div_idx[val_idx])

	# spearman correlation
	x = np.array(x)
	y = np.array(y)
	x = (x - x.mean())/x.std()
	y = (y - y.mean())/y.mean()
	res = stats.spearmanr(x, y)
	print(res.statistic)
	print(res.pvalue)
	# pearson correlation
	res_pearson = stats.pearsonr(x, y)
	print(res_pearson.statistic)
	print(res_pearson.pvalue)

	X = x.reshape(-1,1)
	Y = y.reshape(-1,1)
	reg = LinearRegression().fit(X, Y)
	print(reg.score(X, Y))
	print(reg.coef_)
	print(reg.intercept_)
	# plt.cla()
	# plt.clf()
	colors = sns.husl_palette(l=0.5, s=1.0, h=0.8, as_cmap=True)(np.linspace(0,1,len(biop_sample)))
	df_div = pd.DataFrame.from_dict(df_diversity_idx)
	fig, ax = plt.subplots(figsize=(7, 10))
	g = sns.boxplot(data = df_div, y = 'sample', x = 'diversity_idx', palette = colors)
	plt.tight_layout()
	plt.savefig(os.path.join("./","diversity_idx_boxplots.pdf"))





















		# obj_vals_hists(rna_cells = rna_cells, sample = sample, src_dir = src_dir, original_val = original_val)
		# E_loo_error(df_entire = df_entire, src_dir = src_dir)




