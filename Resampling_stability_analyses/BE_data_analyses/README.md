# Resampling and stability analyses on Barretts' esophagus data set

For all the following experiments, we prepared and preprocessed the DNA and RNA input data according to the instructions that we described in [Instructions for reproducing the results of our analysis on Barretts' esophagus data](https://github.com/NakhlehLab/MaCroDNA/tree/main/BE_data_analysis#instructions-for-reproducing-the-results-of-our-analysis-on-barretts-esophagus-data-set-originally-introduced-in).
Here, we describe the experiments briefly and provide the instructions for running their corresponding codes.

## Random assignment test for BE biopsies
The first experiment is the random assignment of scRNA-seq cells to scDNA-seq cells for each biopsy. For each biopsy, we repeated the random assignment for 100 million times to ensure we have enough number of samples. For each random assignment, we computed the sum of the Pearson correlation coefficients between the paired cells as the *score* of the random trial. Next, we compared the the scores obtained from radnom assignments with that of MaCroDNA's assignment and calculated the p-value of MaCroDNA assignment and plotted the figures.
The script `random_assignment_test.py` performs this test. To run `random_assignment_test.py`:

1. Set `dna_src_dir` to the path to the directory containing the filtered (from the section [Filtering on the scDNA-seq read count tables](https://github.com/NakhlehLab/MaCroDNA/tree/main/BE_data_analysis#filtering-on-the-scdna-seq-read-count-tables)).
2. Set `rna_src_dir` to the path to the directory containing the filtered gene expression data (from the section [Filtering on the scRNA-seq gene expression tables](https://github.com/NakhlehLab/MaCroDNA/blob/main/BE_data_analysis/README.md#filtering-on-the-scrna-seq-gene-expression-tables)).
3. Set `tgt_dir` to the desired path for storing the results.

Please note that, here, we assume the DNA data of each biopsy is a CSV file stored in the `dna_src_dir` named as `<biopsy name>_annotated_filtered_normed_count_table.csv`. Similarly, for the RNA data, all of them are stored in `rna_src_dir` as CSV files named `<biopsy name>__filtered_normed_count_table.csv`.

The outputs of the above script consist of a dictionary containing the scores obtained from MaCroDNA on each biopsy, named `macrodna_objvals.pkl`, the MaCroDNA assignments of cells stored as `<biopsy name>_cell2cell_assignment_indexed.csv`, and all the random scores from the random trials for each biopsy saved as `<biopsy name>_random_samples.npy`. All these files can be found in the `tgt_dir`.

### plotting the results of random assignments 
The code named `plot_random_assignments.py` plots the histograms of the BE biopsies along with the red vertical line indicating the score obtained from MaCroDNA. To run this code, set `src_dir` to the path to the directory where the results of the random assignments are stored (`tgt_dir` in `random_assignment_test.py` code). 

## Stability analysis of MaCroDNA's assignments for BE biopsies
To measure the stability of scRNA-seq cells’ assignments in the BE biopsies, we designed a leave-one-out experiment that involved introducing a small perturbation into the input data by leaving out one of the scRNA-seq cells.
Such perturbation in the RNA inputs might change the assignment of not only the left-out scRNA-seq cell but also the remaining scRNA-seq cells. Therefore, for each scRNA-seq cell, we aggregated the scDNA-seq cells’ IDs assigned to it across all the leave-one-out trials. Next, for each scRNA-seq cell, we collected the scDNA-seq cell IDs that were different from the assigned scDNA-seq cell ID when running MaCroDNA on the entire RNA data and defined the AII as the number of such different scDNA-seq cell IDs. The higher the AII, the more unstable and indefinitive the assignment of a scRNA-seq is in the presence of perturbation.

The script named `run_loo_experiment.py` performs the leave-one-out trials on all the BE biopsies. To run this code:

1. Set `dna_src_dir` to the path to the directory containing the filtered (from the section [Filtering on the scDNA-seq read count tables](https://github.com/NakhlehLab/MaCroDNA/tree/main/BE_data_analysis#filtering-on-the-scdna-seq-read-count-tables)).
2. Set `rna_src_dir` to the path to the directory containing the filtered gene expression data (from the section [Filtering on the scRNA-seq gene expression tables](https://github.com/NakhlehLab/MaCroDNA/blob/main/BE_data_analysis/README.md#filtering-on-the-scrna-seq-gene-expression-tables)).
3. Set `tgt_dir` to the desired path for storing the results. Here, this is set to `./macrodna_res_log_rna_stability/`.

The names of the input CSV files are assumed to follow the same format as in the experiment [Random assignment test for BE biopsies](https://github.com/NakhlehLab/MaCroDNA/blob/main/Resampling_stability_analyses/BE_data_analyses/README.md#random-assignment-test-for-be-biopsies).

### Visualizing the box plots of the assignment instability indices 
The script named `box_plots_from_loo_errors.py` draws the box plots showing the range of values for the assignment instability indices across all scRNA-seq cells for each biopsy. To run this code, set the variables `dna_src_dir` and `rna_src_dir` as for `run_loo_experiment.py` and set `src_dir` to the same path for `tgt_dir` in `run_loo_experiment.py`. The output is a figure named `diversity_idx_boxplots.pdf`.

### Visualizing the correlation between AII and intratumor heterogeneity



