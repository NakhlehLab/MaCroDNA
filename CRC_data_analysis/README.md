# MaCroDNA: Accurate integration of single-cell DNA and RNA data for a deeper understanding of tumor heterogeneity

MaCroDNA (Matching Cross-Domain Nucleic Acids) is a correlation-based method to perform mapping between scRNA-seq gene expression and scDNA-seq copy number values.


## Authors

- Mohammadamin Edrisi
- Xiru Huang
- Huw A. Ogilvie
- Luay Nakhleh

## Usage

1. `Data/`: the test data used in the scripts
2. `macrodna_analysis.ipynb`: a script for preprocessing input data for clonealign, Seurat, and MaCroDNA, as well as running MaCroDNA step by step. The expected outputs are located in the ExpectedOutput/ directory.
3. `ExpectedOutput/`: the expected outputs from `macrodna_analysis.ipynb`
4. `clone_generation.ipynb`: a script explaining how clones were generated using the intNMF and agglomerative clustering methods on untransformed and log-transformed data. The expected generated clusters can be found in the ExpectedCluster/ directory.
5. `ExpectedCluster/`: the expected clusters generated from  `clone_generation.ipynb`
6. `macrodna.py` `seurat_process.R` `utils.py` `run_IntNMF.R`: The functions used in `macrodna_analysis.ipynb` and `clone_generation.ipynb`.
7. `clonealign_performance_analysis/`:  includes the script used to analyze the performance of CloneAlign, as mentioned in section 2 of the supplementary materials titled  **2 Analysis of clonealign performance**



## Installation

**Python package requirements**

MaCroDNA is written on python3 (tested on verision 3.7.7) and the following packages are needed to run MaCroDNA.
- numpy (tested on verision 1.21.5)
- pandas (tested on verision 1.3.5)
- scipy (tested on verision 1.7.3)
- gurobipy (tested on verision 9.5.1)

To install gurobipy, Gurobi need to be installed. The academic license is free.The installation instructions for
Gurobi are

| OS      | Instruction |
|:--------|:------------|
| Linux | [Gurobi Installation Guide: Linux](https://youtu.be/yNmeG6Wom1o) |
| Windows | [Gurobi Installation Guide: Windows](https://youtu.be/fQVxuWOiPpI) |
| macOS | [Gurobi Installation Guide: macOS](https://youtu.be/ZcL-NmckTxQ) |

After the installation of Gurobi, gurobipy can be installed using **pip**

`python -m pip install gurobipy`

Other installation methods can be found 
[How do I install Gurobi for Python](https://support.gurobi.com/hc/en-us/articles/360044290292-How-do-I-install-Gurobi-for-Python-)
    

**R libraries requirements**

MaCroDNA doesn't contain any R code, but the data processing procedure needs to use the functions in Seurat. The R library needed is
- Seurat

**Typical install time**

The install time mainly comes from installing gurobipy (~10 mins) and Seurat (~8 mins). The total install time is around 20 minutes.


