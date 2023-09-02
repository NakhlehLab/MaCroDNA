import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import numpy as np 
import seaborn as sns
import os
import pickle

from scipy import stats
import pandas as pd
from scipy import spatial
from scipy.stats import spearmanr
# sns.set_style("whitegrid", {"axes.facecolor": ".87"})

def label_point(x, y, val, ax):
	a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
	for i, point in a.iterrows():
		ax.text(point['x'], point['y'], str(point['val']), size=1)

if __name__=="__main__":

	biop_sample = ["PAT20_CARD", "PAT20_ESO", "PAT9_NDBE", "PAT14_NDBE", "PAT16_NDBE", "PAT6_LGD", "PAT19_LGD", "PAT6_HGD", "PAT14_HGD", "PAT20_HGD1", "PAT16_EAC"]

	src_dir_meidan = "./macrodna_res_log_random_test_median"
	src_dir_sum = "./macrodna_res_log_random_test"

	sum_arr = []
	median_arr = []

	for bio in biop_sample:

		a_sum = np.load(os.path.join(src_dir_sum, bio)+"_random_samples.npy")
		a_median = np.load(os.path.join(src_dir_meidan, bio)+"_random_medians.npy")

		sum_arr.append(np.mean(a_sum))
		median_arr.append(np.mean(a_median))

	# spearman correlation
	x = np.array(sum_arr)
	y = np.array(median_arr)
	res = stats.spearmanr(x, y)
	
	print(res.correlation)
	print(res.pvalue)
	# pearson correlation
	res_pearson = stats.pearsonr(x, y)
	print(res_pearson)
	print(res_pearson[0])
	print(res_pearson[1])


	df_ = pd.DataFrame({'sum': x, 'median': y, 'sample': biop_sample})
	
	g = sns.lmplot(x="sum", y="median", data=df_, fit_reg=True)
	label_point(df_['sum'], df_['median'], df_['sample'], plt.gca())
	g.set(title = f"pearson: {res_pearson[0]:.3f}, {res_pearson[1]:.3f}\nspearman: {res.correlation:.3f}, {res.pvalue:.3f}")
	plt.tight_layout()

	plt.savefig(os.path.join("./","sum_median_correlation_plot.pdf")
		, transparent=True)

	
	df_ = df_.sort_values('sum', ascending=False)
	print(df_)







