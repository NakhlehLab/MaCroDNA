# CRC data analysis

The aim of this project is to assess the performance of integration methods for analyzing CRC data.

## Introduction

Colorectal cancer is a prevalent and deadly disease, and understanding its molecular characteristics is crucial for developing effective treatment strategies. This project evaluates different integration methods by analyzing a CRC dataset introduced by [Bian et al](https://www.science.org/doi/10.1126/science.aao3791?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed). The dataset provides both RNA and DNA measurements for certain cells, which are used as ground truth data for evaluation.


## Repository Structure

The repository is structured as follows:
1. `Data/`: the test data to run the scripts
2. `macrodna_analysis.ipynb`: a script for preprocessing input data for clonealign, Seurat, and MaCroDNA, as well as running MaCroDNA step by step. The expected outputs are located in the ExpectedOutput/ directory.
3. `ExpectedOutput/`: the expected outputs from `macrodna_analysis.ipynb`
4. `clone_generation.ipynb`: a script explaining how clones were generated using the intNMF and agglomerative clustering methods on untransformed and log-transformed data. The expected generated clusters can be found in the ExpectedCluster/ directory.
5. `ExpectedCluster/`: the expected clusters generated from  `clone_generation.ipynb`
6. `macrodna.py` `seurat_process.R` `utils.py` `run_IntNMF.R`: The functions used in `macrodna_analysis.ipynb` and `clone_generation.ipynb`.
6. `OtherIntegrationMethods/`: includes the scripts of Seurat, clonealign, and baseline that we used in the paper
8. `clonealign_performance_analysis/`:  includes the script used to analyze the performance of CloneAlign, as mentioned in section 2 of the supplementary materials titled  **2 Analysis of clonealign performance**


## Usage
1. Run the macrodna_analysis.ipynb script to preprocess the data and run MaCroDNA. The expected outputs can be found in the `ExpectedOutput/` directory.
2. Run the clone_generation.ipynb script to generate clones using intNMF and agglomerative clustering methods. The expected clusters can be found in the `ExpectedCluster/` directory.
