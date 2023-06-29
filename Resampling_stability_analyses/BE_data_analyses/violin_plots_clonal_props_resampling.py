import pandas as pd
import numpy as np 
import os
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("whitegrid")
sns.set_context('talk')
sns.despine(left=True)



if __name__=="__main__":
	biop_sample = ["crc04", "crc10", "crc11"]
	options = ['constrained']
	clustering_settings = ["agglomerative_log", "agglomerative_raw", "intNMF_log", "intNMF_raw"]
	src_dir = "/Users/edrisi/Documents/ongoing_projects/MaCroDNA/phylosignal_analysis/Bian_scripts/all_props_exp_res/"

	# for loop on all possible options
	
	for setting in clustering_settings:
		accs = []
		patient_labels = []
		option_labels = []
		for patient_biop in biop_sample:
			for option in options:
				df = pd.read_csv(os.path.join(src_dir, f"props_exp_res_{setting}_{option}_{patient_biop}/dat.csv"), index_col=0)
				df = df.dropna(thresh=1)
				vals = df[df['dat_type'] == "resampled"]['accuracy'].values.tolist()
				accs += vals
				patient_labels += [patient_biop for _ in range(len(vals))]
				option_labels += [option for _ in range(len(vals))]
		print(len(accs))
		df_ = pd.DataFrame({'Accuracy': accs, 'Patient': patient_labels, 'Option': option_labels})
		plt.cla()
		plt.clf()
		g = sns.violinplot(data=df_, y="Accuracy", x="Patient", cut=0, scale = "width")
		g.set(title = setting)
		g.set_ylim(-0.05, 1.05)
		plt.tight_layout()
		plt.savefig(f"./{setting}_patients_options_prop_exp.pdf")

