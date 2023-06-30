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
