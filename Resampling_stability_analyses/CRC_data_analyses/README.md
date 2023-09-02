# Resampling and stability analyses on colorectal cancer (CRC) data set

## Investigating the effect of imbalance in clonal proportions across modalities in CRC patients
We investigated the effect of deviation from this assumption on the performance of MaCroDNA on the CRC data set. In particular, for each CRC patient, we performed two experiments that involved resampling the original data sets, resulting in variations in clonal proportions between scDNA-seq and scRNA-seq data.

### Random removal of scDNA-seq cells from DNA data
In this experiment, for each CRC patient, we performed the following test:

1. From the scDNA-seq cells whose joint RNA information is available, randomly select 10 cells and remove them from the original DNA data.
2. Run MaCroDNA on the new data set and measure the accuracy of clonal assignments only for the scRNA-seq cells whose actual scDNA-seq pairs were removed from the DNA data in the previous step.

We repeated the above test 10,000 times on the preprocessed CRC data named all_genes_log with different clustering settings. For the details on how to prepare the all_genes_log data, see the Jupyter notebook [macrodna_analysis.ipynb](https://github.com/NakhlehLab/MaCroDNA/blob/main/CRC_data_analysis/macrodna_analysis.ipynb). We also have provided the details for obtaining the results of different clustering algorithms in the Jupyter notebook [clone_generation.ipynb](https://github.com/NakhlehLab/MaCroDNA/blob/main/CRC_data_analysis/clone_generation.ipynb).

The script named `run_dna_batch_removal_exp.py` performs the above experiment on all CRC patients. To run this code:

1. Set `dna_src_dir` to the path to the directory containing the DNA data of all CRC patients.
2. Set `rna_src_dir` to the path to the directory containing the RNA data of all CRC patients. Here, we have stored both DNA and RNA data in a directory named `all_genes_log`.
3. Set the variable `mother_dir` to a path you would like to store the outputs.
4. Set `clus_path` to the directory where the different clustering results are stored. Here, this variable is set to `Clusters`.

Note that we assume that each DNA or RNA data is in the form of a CSV file named `<patient ID>_dna.csv` or `<patient ID>_rna.csv`. The specific patient IDs we used are provided in the Python list named `biop_sample`. Also, the names of the clustering data can be found in the Python list named `clustering_settings`.

After running the code, for each clustering setting and each patient, it creates a directory in the `mother_dir` with the format of `dna_removal_exp_res_<clustering setting>_<patient ID>`. Inside each subdirectory, a CSV file named `dat.csv` can be found that contains the accuracy of the assignments of the randomly selected RNA cells, before and after the removal of their DNA pairs.

The sample outputs of this code for 10 trials on each CRC patient and each clustering setting can be found in the directory [all_dna_removal_exp](https://github.com/NakhlehLab/MaCroDNA/tree/main/Resampling_stability_analyses/CRC_data_analyses/sample_outputs/all_dna_removal_exp).

#### Visualizing the results of the random removal experiment
We drew two sets of plots for this experiment: one was the box plots showing the accuracy of clonal assignments of the randomly selected RNA cells whose true DNA pairs were removed, before and after their removal. The second was the box plots showing the distribution of the changes in accuracy before and after the removal of DNA cells in each trial. 
The script named `box_plots_dna_batch_removal.py` draws these plots for each clustering setting. To run this code, set `src_dir` to the same path as `mother_dir` in `run_dna_batch_removal_exp.py`. 


### Resampling clonal proportions in DNA data

The previous experiment by removal of scDNA-seq cells does not provide insights into the behavior of MaCroDNA under extreme conditions where the proportions of clones deviate from the original ones significantly. To address this, we designed another experiment on the CRC patients, in which we randomly drew the clonal proportions in the DNA data and resampled the scDNA-seq cells with replacement within each clone to achieve the drawn proportions. We repeated this trial 10,000 times on each CRC patient.

The script `clonal_proportions_resampling.py` performs the above experiment. To run this code (similar to `run_dna_batch_removal_exp.py`):

1. Set `dna_src_dir` to the path to the directory containing the DNA data of all CRC patients.
2. Set `rna_src_dir` to the path to the directory containing the RNA data of all CRC patients. Here, we have stored both DNA and RNA data in a directory named `all_genes_log`.
3. Set the variable `mother_dir` to a path you would like to store the outputs.
4. Set `clus_path` to the directory where the different clustering results are stored. Here, this variable is set to `Clusters`.

After running the code, for each clustering setting and each patient, it creates a directory in the `mother_dir` with the format of `prop_exp_res_<clustering setting>_<patient ID>`. Inside each subdirectory, a CSV file named `dat.csv` can be found that contains the accuracy of the assignments of the RNA cells (for which we know the true DNA pairs) along with the randomly drawn clonal proportions. Also, you can find the all DNA cell IDs present in each clone in another subdirectory named `names` in the same directory.

The sample outputs of this code for 10 trials on all CRC patients and clustering settings can be found in the directory [all_props_exp_res](https://github.com/NakhlehLab/MaCroDNA/tree/main/Resampling_stability_analyses/CRC_data_analyses/sample_outputs/all_props_exp_res).

#### Visualizing the accuracy measurements of the clonal proportions' resampling experiment

The script `violin_plots_clonal_props_resampling.py` plots the violin plots of the accuracy measurements for each clustering setting separately. To run this code, set `src_dir` to the same path as `mother_dir` in `clonal_proportions_resampling.py`.

#### Visualizing the effect of the deviation from the original clonal proportions and the accuracy

To better understand the impact of the difference between the clonal proportions in the two modalities on the clonal assignment accuracy, we calculated the Earth mover's distance (EMD) between the original clonal proportions and the sampled ones for each random test. Moreover, we calculated the Spearman correlation between the EMDs and the corresponding accuracy values for each patient and clustering setting in the script named `emd_vs_acc_clonal_props_resampling.py`. To run this code, set `src_dir` to the same path as `mother_dir` in `clonal_proportions_resampling.py`.








