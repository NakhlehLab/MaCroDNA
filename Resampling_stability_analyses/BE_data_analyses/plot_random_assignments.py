import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import numpy as np 
import seaborn as sns
import os
import pickle



if __name__=="__main__":
	biop_sample = ["PAT20_CARD", "PAT20_ESO", "PAT9_NDBE", "PAT14_NDBE", "PAT16_NDBE", "PAT6_LGD", "PAT19_LGD", "PAT6_HGD", "PAT14_HGD", "PAT20_HGD1", "PAT16_EAC"]
	src_dir = "./macrodna_res_log_random_test"

	# read the actual objective functions 
	with open(os.path.join(src_dir, "macrodna_objvals.pkl"), "rb") as f:
		d = pickle.load(f)
	print(d)

	for bio in biop_sample:

		a = np.load(os.path.join(src_dir, bio)+"_random_samples.npy")
		print(a.shape)
		# calculate the pvalue of the real score 
		x_ = len(np.argwhere(a >= d[bio]))
		print(f"number of random scores greater than or equal to the real score {x_}")
		pval = x_ / len(a)
		fig, ax = plt.subplots()

		ax.hist(a, bins=50, density=False)
		ax.axvline(d[bio], ls="--", color="r")
		score_label = f"Score on original\ndata: {d[bio]:.3f}\n(p-value: {pval:.3f})"
		ymin, ymax = ax.get_ylim()
		ax.text(d[bio] - 2.6, ymax/2, score_label, fontsize=12)
		ax.set_xlabel("Score (sum of the Pearson correlation values of pairs)")
		_ = ax.set_ylabel("Counts")
		ax.set_title(bio.replace("_", " (").replace("PAT", "Patient ")+")")
		plt.savefig(os.path.join(src_dir, bio)+"random_hist.pdf"
			, transparent=True)



