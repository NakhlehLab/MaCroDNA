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
	src_dir = "/Users/edrisi/Documents/ongoing_projects/MaCroDNA/phylosignal_analysis/Bian_scripts/all_dna_removal_exp/"

	# for loop on all possible options
	
	for setting in clustering_settings:
		accs = []
		accs_labels  = []
		patient_labels = []
		option_labels = []

		diffs = []
		diffs_patient_labels = []
		diffs_option_labels = []

		for patient_biop in biop_sample:
			for option in options:
				df = pd.read_csv(os.path.join(src_dir, f"dna_removal_exp_res_{setting}_{option}_{patient_biop}/dat.csv"), index_col=0)
				df = df.dropna(thresh=1)

				vals = df['original_acc'].values.tolist()
				accs += vals
				accs_labels += ['original' for _ in range(len(vals))]
				patient_labels += [patient_biop for _ in range(len(vals))]
				option_labels += [option for _ in range(len(vals))]

				vals = df['dropped_acc'].values.tolist()
				accs += vals
				accs_labels += ['subsampled' for _ in range(len(vals))]
				patient_labels += [patient_biop for _ in range(len(vals))]
				option_labels += [option for _ in range(len(vals))]


				# ========= the differences between the two columns =========
				diff_df = df['original_acc'] - df['dropped_acc']
				diffs += diff_df.values.tolist()
				diffs_patient_labels += [patient_biop for _ in range(len(diff_df.values.tolist()))]
				diffs_option_labels += [option for _ in range(len(diff_df.values.tolist()))]

		# ========= plot the violins of the accuracies of the dropped batches in the original and resampled data =========
		df_ = pd.DataFrame({'Accuracy': accs, 'Data_label': accs_labels, 'Patient': patient_labels, 'Option': option_labels})
		plt.cla()
		plt.clf()
		g = sns.boxplot(data = df_, y = "Accuracy", x = "Patient", hue = "Data_label")
		g.set(title = setting)
		plt.tight_layout()
		plt.savefig(f"./{setting}_patients_options_dna_batch_removal_exp_boxplot.pdf")

		# ========= plot the violins showing the difference in the accuracies of the dropped batches before and after drop out =========
		plt.cla()
		plt.clf()
		df = pd.DataFrame({'Change in accuracy': diffs, 'Patient': diffs_patient_labels, 'Option': diffs_option_labels})
		g = sns.boxplot(data = df, y = 'Change in accuracy', x = 'Patient')
		g.set(title = setting)
		plt.tight_layout()
		plt.savefig(f"./{setting}_patients_options_dna_batch_removal_exp_changes_boxplot.pdf")




















