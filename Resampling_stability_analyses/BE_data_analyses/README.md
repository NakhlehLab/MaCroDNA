# Resampling and stability analyses on Barretts' esophagus data set

For all the following experiments, we prepared and preprocessed the DNA and RNA input data according to the instructions that we described in [Instructions for reproducing the results of our analysis on Barretts' esophagus data](https://github.com/NakhlehLab/MaCroDNA/tree/main/BE_data_analysis#instructions-for-reproducing-the-results-of-our-analysis-on-barretts-esophagus-data-set-originally-introduced-in).

## Random assignment test for BE biopsies
In this experiment, we randomly assigned the scRNA-seq cells to scDNA-seq cells in each biopsy. For each biopsy, we repeated the random assignment for 100 million times. 
The script `random_assignment_test.py` performs this test. To `random_assignment_test.py`:

1. Set `dna_src_dir` to the path to the directory containing the filtered (from the section [Filtering on the scDNA-seq read count tables](https://github.com/NakhlehLab/MaCroDNA/tree/main/BE_data_analysis#filtering-on-the-scdna-seq-read-count-tables)).
2. Set `rna_src_dir` to the path to the directory containing the filtered gene expression data (from the section [Filtering on the scRNA-seq gene expression tables](https://github.com/NakhlehLab/MaCroDNA/blob/main/BE_data_analysis/README.md#filtering-on-the-scrna-seq-gene-expression-tables)).
