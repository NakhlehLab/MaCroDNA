import pandas as pd
import numpy as np 
import os
import scipy
from scipy import stats
from scipy.stats import wasserstein_distance
import sklearn
from sklearn.linear_model import LinearRegression
from numpy.linalg import norm
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("whitegrid")
sns.set_context('talk')



if __name__=="__main__":
	biop_sample = ["crc04", "crc10", "crc11"]
	options = ['constrained']
	clustering_settings = ["agglomerative_log", "agglomerative_raw", "intNMF_log", "intNMF_raw"]
	src_dir = "/Users/edrisi/Documents/ongoing_projects/MaCroDNA/phylosignal_analysis/Bian_scripts/all_props_exp_res/"

	# for loop on all possible options
	
	for setting in clustering_settings:
		for patient_biop in biop_sample:
			for option in options:
				dists = []
				df = pd.read_csv(os.path.join(src_dir, f"props_exp_res_{setting}_{option}_{patient_biop}/dat.csv"), index_col=0)
				vals = df[df['dat_type'] == "resampled"]['accuracy'].values.tolist()
				accs = vals

				dna_cols = [x for x in df.columns if x.startswith("clone_")]

				dna_vals = df.loc[1:,dna_cols].values
				ref_vals = df.loc[0, dna_cols].values

				print(dna_vals.shape)
				print(ref_vals)
				for idx in range(dna_vals.shape[0]):
					dist = wasserstein_distance(np.arange(len(dna_vals[idx,:])), np.arange(len(ref_vals)),
						dna_vals[idx,:], ref_vals)

					dists.append(dist)

				df_ = pd.DataFrame({'Accuracy': accs, 'Distance': dists})
				X = np.array(dists).reshape(-1,1)
				y = np.array(accs).reshape(-1,1)
				reg = LinearRegression().fit(X, y)

				# spearman correlation
				res = stats.spearmanr(dists, accs)
				# pearson correlation
				res_pearson = stats.pearsonr(dists, accs)
				# rank the data 
				df_ranked = df_.rank()

				plt.cla()
				plt.clf()
				g = sns.scatterplot(data = df_, x="Distance", y="Accuracy", s = 10, linewidth=0.0)
				g.set(title = f"{setting}\n {res.statistic:.3f}\n{res.pvalue:.3e}")
				g.set_ylim(-0.05, 1.05)
				plt.tight_layout()
				plt.savefig(f"./props_exp_emd_vs_acc_{setting}_{option}_{patient_biop}.pdf")


























