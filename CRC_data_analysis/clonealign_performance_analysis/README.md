# clonealign performance analysis


This repository contains the analysis of clonealign's performance in the context of clones generated from agglomerative clustering on log-transformed data(agg-log clones). The results indicate that clonealign's performance was unsatisfactory when using agg-log clones, as it assigned all cells to a single clone, disregarding the inherent relation between copy number and gene expression.

## Background
In agg-log clones in the CRC data, it was observed that each patient had a copy-number-2 clone, where the copy number for every gene is 2. This is the clone that clonealign prefers assigning all the cells to.

It is important to note that the input of clonealign is different from that of Seurat and MaCroDNA. The input of Seurat and MaCroDNA consists of three parts: a gene expression read count matrix, a cell copy number matrix, and the clone labels for all cells in the copy number matrix. However, clonealign only needs a gene expression read count matrix and a clone copy number matrix. The gene expression read count matrix is a cell-by-gene matrix containing the gene expression read count for each gene in each RNA cell. The cell copy number matrix is also a cell-by-gene matrix, and each entry is the copy number of one gene in one DNA cell. The clone copy number matrix, on the other hand, is a clone-by-gene matrix, with each entry representing the copy number of one gene in one clone. To obtain the copy number of one gene in one clone, we calculated the median value of that gene's copy numbers in all the cells in that clone. Even though the copy number is not 2 for each gene for each cell in that clone, the calculation of the median number results in the copy-number-2 clone, which serves as the input of clonealign.

## Further analysis

To investigate the performance issues observed with clonealign, the following experiments were conducted in the `clonealign_investigate.R` script:

1. **Experiment 0**: Original clonealign results and results after removing the copy-number-2 clones.
2. **Experiment 1.1 single uniform-copy-number clone**: After removing the copy-number-2 clones, a copy-number-1 clone was added to each patient.
2. **Experiment 1.2 single uniform-copy-number clone**: After removing the copy-number-2 clones, a copy-number-5 clone was added to each patient.
3. **Experiment 2.1 dominant-copy-number**: After removing the copy-number-2 clones, add a dominant-copy-number clone where 90% genes have copy number 2
4. **Experiment 2.2 dominant-copy-number**: After removing the copy-number-2 clones, add a dominant-copy-number clone where 80% genes have copy number 2
4. **Experiment 2.2 dominant-copy-number**: After removing the copy-number-2 clones, add a dominant-copy-number clone where 70% genes have copy number 2
4. **Experiment 3 multiple uniform-copy-number clones**: Copy-number-2 clones were retained, and a copy-number-5 clone was added. This test was repeated multiple times with different seeds.

, where
- uniform-copy-number clone: the clone where all the genes have the same copy number
- copy-number-x clone: the clone where all the genes have copy number x
- dominant-copy-number clone: the clone that has a significantly high proportion of its most common copy number

## Usage
To run the analysis, execute the `clonealign_investigate.R`.

For a more detailed understanding of the experiments and their results, please refer to the section 2 in the supplementary material, titled as "2 Analysis of clonealign performance"