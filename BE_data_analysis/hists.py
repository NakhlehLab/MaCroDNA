import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rc,rcParams
sns.set_style("darkgrid", {"axes.facecolor": ".87"})
# sns.set_style("darkgrid")
# sns.set_style("whitegrid")
import os

def argsort(seq):
	# http://stackoverflow.com/questions/3071415/efficient-method-to-calculate-the-rank-vector-of-a-list-in-python
	return sorted(range(len(seq)), key=seq.__getitem__)


if __name__=="__main__":
	df_cancer = pd.read_csv("./data/tier1_gene.csv", index_col=0)
	cancer_genes = df_cancer.values.tolist()
	cancer_genes = [x[0] for x in cancer_genes]
	sample_list = ["PAT20_CARD", "PAT20_ESO", "PAT9_NDBE", "PAT14_NDBE", "PAT16_NDBE", "PAT6_LGD", "PAT19_LGD", "PAT6_HGD", "PAT14_HGD", "PAT20_HGD1", "PAT16_EAC"]
	src_dir = f"./data/Busslinger_data/phylosignal_res_K.star"
	dna_src_dir = "./data/Busslinger_data/macrodna_paired/"
	q = []
	r = []
	cellvalues = []
	cellvalues_cosmic = []
	gene_set_arr = []
	gene_set_arr_names = []
	overlap_ = []
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
		idx = np.argwhere(p<0.05)

		q.append(s)
		# q.append(np.log10(s))
		# q.append(s[idx[:,0]])

	r = q[0]
	for i in range(len(sample_list[1:])):
		r = np.hstack((r, q[i]))

	##############################################################################
	# plot the histograms
	plt.cla()
	plt.clf()

	n = len(sample_list)
	# colors = sns.dark_palette("royalblue", as_cmap=True)(np.linspace(0,1,n))
	# colors = sns.color_palette("husl", 11, as_cmap = True)(np.linspace(0,1,n))
	# colors = sns.cubehelix_palette(hue=1.0, as_cmap=True)(np.linspace(0,1,n))
	# colors = sns.husl_palette(l=0.5, s=1.0, h=0.8, as_cmap=True)(np.linspace(0,1,n))
	lw = 8
	title_fontsize = 24
	stat_ = "probability"
	ylabel_ = "Density"
	xlabel_ = r"${\mathrm{K}^{\ast}}$"
	ylabel_fontsize = 25
	xlabel_fontsize = 25
	tick_size = 23
	tick_size_y = 23
	ylim_ = 0.91
	bins_arr = [15,5,5,10,5,10,5,30,20,20,20]
	# bins_arr = [10 for _ in range(len(sample_list))]
	plt.rcParams["font.weight"] = "bold"
	plt.rcParams["axes.labelweight"] = "bold"
	# ylim_ = 0.4
	# bins_arr = [50 for _ in range(len(sample_list))]
	f, ax = plt.subplots(7, 2, figsize=(16, 19), sharex=True, constrained_layout = True, sharey = False)
	f.delaxes(ax[4,1])
	f.delaxes(ax[5,1])
	f.delaxes(ax[6,1])
	# f, ax = plt.subplots(11, 1, figsize=(8, 25), sharex=True)
	# g = sns.kdeplot(q[0], color=colors[0],label=sample_list[0], ax=ax[0,0])
	g = sns.histplot(q[0], stat=stat_, color=colors[0], kde=True, label=sample_list[0], ax=ax[0,0], line_kws=dict(linewidth=lw), bins=bins_arr[0])
	ax[0,0].set_title(sample_list[0].replace("_", " (").replace("PAT", "Patient ")+")", fontweight="bold", fontsize=title_fontsize)
	ax[0,0].tick_params(axis='x', rotation=0, labelsize=tick_size)
	ax[0,0].tick_params(axis='y', rotation=0, labelsize=tick_size_y)
	ax[0,0].set_ylim(0, ylim_)
	# g.lines[0].set_color('black')
	g.set_ylabel(ylabel_, fontsize=ylabel_fontsize)
	g = sns.histplot(q[1], stat=stat_, color=colors[1], kde=True, label=sample_list[1], ax=ax[1,0], line_kws=dict(linewidth=lw), bins=bins_arr[1])
	ax[1,0].set_title(sample_list[1].replace("_", " (").replace("PAT", "Patient ")+")", fontweight="bold", fontsize=title_fontsize)
	ax[1,0].tick_params(axis='x', rotation=0, labelsize=tick_size)
	ax[1,0].tick_params(axis='y', rotation=0, labelsize=tick_size_y)
	ax[1,0].set_ylim(0, ylim_)
	g.set_ylabel(ylabel_, fontsize=ylabel_fontsize)
	g = sns.histplot(q[2], stat=stat_, color=colors[2], kde=True, label=sample_list[2], ax=ax[2,0], line_kws=dict(linewidth=lw), bins=bins_arr[2])
	ax[2,0].set_title(sample_list[2].replace("_", " (").replace("PAT", "Patient ")+")", fontweight="bold", fontsize=title_fontsize)
	ax[2,0].tick_params(axis='x', rotation=0, labelsize=tick_size)
	ax[2,0].tick_params(axis='y', rotation=0, labelsize=tick_size_y)
	ax[2,0].set_ylim(0, ylim_)
	g.set_ylabel(ylabel_, fontsize=ylabel_fontsize)
	g = sns.histplot(q[3], stat=stat_, color=colors[3], kde=True, label=sample_list[3], ax=ax[3,0], line_kws=dict(linewidth=lw), bins=bins_arr[3])
	ax[3,0].set_title(sample_list[3].replace("_", " (").replace("PAT", "Patient ")+")", fontweight="bold", fontsize=title_fontsize)
	ax[3,0].tick_params(axis='x', rotation=0, labelsize=tick_size)
	ax[3,0].tick_params(axis='y', rotation=0, labelsize=tick_size_y)
	ax[3,0].set_ylim(0, ylim_)
	g.set_ylabel(ylabel_, fontsize=ylabel_fontsize)
	g = sns.histplot(q[4], stat=stat_, color=colors[4], kde=True, label=sample_list[4], ax=ax[4,0], line_kws=dict(linewidth=lw), bins=bins_arr[4])
	ax[4,0].set_title(sample_list[4].replace("_", " (").replace("PAT", "Patient ")+")", fontweight="bold", fontsize=title_fontsize)
	ax[4,0].tick_params(axis='x', rotation=0, labelsize=tick_size)
	ax[4,0].tick_params(axis='y', rotation=0, labelsize=tick_size_y)
	ax[4,0].set_ylim(0, ylim_)
	g.set_ylabel(ylabel_, fontsize=ylabel_fontsize)
	g = sns.histplot(q[5], stat=stat_, color=colors[5], kde=True, label=sample_list[5], ax=ax[5,0], line_kws=dict(linewidth=lw), bins=bins_arr[5])
	ax[5,0].set_title(sample_list[5].replace("_", " (").replace("PAT", "Patient ")+")", fontweight="bold", fontsize=title_fontsize)
	ax[5,0].tick_params(axis='x', rotation=0, labelsize=tick_size)
	ax[5,0].tick_params(axis='y', rotation=0, labelsize=tick_size_y)
	ax[5,0].set_ylim(0, ylim_)
	g.set_ylabel(ylabel_, fontsize=ylabel_fontsize)
	g = sns.histplot(q[6], stat=stat_, color=colors[6], kde=True, label=sample_list[6], ax=ax[6,0], line_kws=dict(linewidth=lw), bins=bins_arr[6])
	ax[6,0].set_title(sample_list[6].replace("_", " (").replace("PAT", "Patient ")+")", fontweight="bold", fontsize=title_fontsize)
	ax[6,0].tick_params(axis='x', rotation=0, labelsize=tick_size)
	ax[6,0].tick_params(axis='y', rotation=0, labelsize=tick_size_y)
	ax[6,0].set_ylim(0, ylim_)
	g.set_ylabel(ylabel_, fontsize=ylabel_fontsize)
	g.set_xlabel(xlabel_, fontsize=xlabel_fontsize, labelpad=8, fontweight="bold")
	g = sns.histplot(q[7], stat=stat_, color=colors[7], kde=True, label=sample_list[7], ax=ax[0,1], line_kws=dict(linewidth=lw), bins=bins_arr[7])
	ax[0,1].set_title(sample_list[7].replace("_", " (").replace("PAT", "Patient ")+")", fontweight="bold", fontsize=title_fontsize)
	ax[0,1].tick_params(axis='x', rotation=0, labelsize=tick_size)
	ax[0,1].tick_params(axis='y', rotation=0, labelsize=tick_size_y)
	ax[0,1].set_ylim(0, ylim_)
	g.set_ylabel(ylabel_, fontsize=ylabel_fontsize)
	g = sns.histplot(q[8], stat=stat_, color=colors[8], kde=True, label=sample_list[8], ax=ax[1,1], line_kws=dict(linewidth=lw), bins=bins_arr[8])
	ax[1,1].set_title(sample_list[8].replace("_", " (").replace("PAT", "Patient ")+")", fontweight="bold", fontsize=title_fontsize)
	ax[1,1].tick_params(axis='x', rotation=0, labelsize=tick_size)
	ax[1,1].tick_params(axis='y', rotation=0, labelsize=tick_size_y)
	ax[1,1].set_ylim(0, ylim_)
	g.set_ylabel(ylabel_, fontsize=ylabel_fontsize)
	g = sns.histplot(q[9], stat=stat_, color=colors[9], kde=True, label=sample_list[9], ax=ax[2,1], line_kws=dict(linewidth=lw), bins=bins_arr[9])
	ax[2,1].set_title(sample_list[9].replace("_", " (").replace("PAT", "Patient ")+")", fontweight="bold", fontsize=title_fontsize)
	ax[2,1].tick_params(axis='x', rotation=0, labelsize=tick_size)
	ax[2,1].tick_params(axis='y', rotation=0, labelsize=tick_size_y)
	ax[2,1].set_ylim(0, ylim_)
	g.set_ylabel(ylabel_, fontsize=ylabel_fontsize)
	g = sns.histplot(q[10], stat=stat_, color=colors[10], kde=True, label=sample_list[10], ax=ax[3,1], line_kws=dict(linewidth=lw), bins=bins_arr[10])
	ax[3,1].set_title(sample_list[10].replace("_", " (").replace("PAT", "Patient ")+")", fontweight="bold", fontsize=title_fontsize)
	ax[3,1].xaxis.set_tick_params(which='both', labelbottom=True)
	ax[3,1].tick_params(axis='x', rotation=0, labelsize=tick_size)
	ax[3,1].tick_params(axis='y', rotation=0, labelsize=tick_size_y)
	ax[3,1].set_ylim(0, ylim_)
	g.set_ylabel(ylabel_, fontsize=ylabel_fontsize)
	g.set_xlabel(xlabel_, fontsize=xlabel_fontsize)
	# plt.setp(ax[3, 1], xlabel=xlabel_)
	# ax[3,1].get_xaxis().set_visible(True)
	plt.savefig(f'dists_Kstar_idx_sep.pdf')
	# plt.show()

	plt.cla()
	plt.clf()
	sns.set_style("whitegrid")

	lw = 5
	ylabel_fontsize = 26
	xlabel_fontsize = 26
	tick_size_y = 25
	tick_size_x = 25
	# ubs = [1.3, 1.025, 1.05, 1.025, 1.025, 1.1, 1.05, 2.2, 1.45, 1.5, 1.65]
	ubs = [2.2, 2.2, 2.2, 2.2, 2.2, 2.2, 2.2, 2.2, 2.2, 2.2, 2.2]
	ubs_y = [0.03, 0.005, 0.06, 0.005, 0.01, 0.02, 0.004, 0.15, 0.1, 0.15, 0.04]
	# bins_ = [20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20]
	bins_ = [10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10]
	for i in range(len(sample_list)):
		# f, ax = plt.subplots(1, 1, figsize=(4, 4), constrained_layout=False)
		f, ax = plt.subplots(1, 1, figsize=(5.5, 2), constrained_layout=False)
		g = sns.histplot(q[i], stat=stat_, color=colors[i], kde=True, label=sample_list[i], ax=ax, line_kws=dict(linewidth=lw), bins=bins_[i])
		# ax.set_title(sample_list[i].replace("_", " (").replace("PAT", "Patient ")+")", fontweight="bold", fontsize=title_fontsize)
		ax.tick_params(axis='x', rotation=0, labelsize=tick_size_x)
		ax.tick_params(axis='y', rotation=0, labelsize=tick_size_y)
		# plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
		g.set_ylabel(ylabel_, fontsize=ylabel_fontsize)
		g.set_xlabel(xlabel_, fontsize=xlabel_fontsize)
		ax.xaxis.set_label_coords(1.1, 0.1)

		# g.lines[0].set_color('black')
		# g.set_ylabel(ylabel_, fontsize=ylabel_fontsize)
		# g.set(yticklabels=[])
		# g.set(ylabel=None)
		plt.xlim(0.99, 1.5)
		plt.ylim(0, 0.01)
		plt.margins(x=0)
		plt.margins(y=0)
		plt.tight_layout()
		plt.savefig(f'{sample_list[i]}_magnified.pdf')









