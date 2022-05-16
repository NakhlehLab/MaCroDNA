# MaCroDNA

By Mohammadamin Edrisi, Xiru Huang, Huw A. Ogilvie and Luay Nakhleh

Department of Computer Science, Rice University

### Table of Contents
[1. Introduction](#Introduction)  

[2. Citation](#Citation)  

[3. Installation](#Installation)  

[4. Usage](#Usage)

[5. Contacts](#contacts)



### Introduction
MaCroDNA (**Ma**tching **Cro**ss-**D**omain **N**ucleic **A**cids) is a correlation-based 
method to perform mapping between scRNA-seq gene expression and scDNA-seq copy number values. 
This repository contains the source code for MaCroDNA described in the paper "MaCroDNA: Accurate integration of single-cell DNA and RNA
data for a deeper understanding of tumor heterogeneity".

### Citation


### Installation
**Package Requirements**
>Python >= 3.7
> 
>numpy
> 
>pandas
> 
>scipy
> 
>math
> 
> gurobipy

Here, only gurobipy need to be installed manually before installing MaCroDNA. Other packages can be installed automatically
 while installing MaCroDNA

**gurobipy Installation**

First, Gurobi need to be installed. The academic license is free. The installation instructions for
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

**MaCroDNA Installation**

`pip install -i https://test.pypi.org/simple/ MaCroDNA==0.0.5 `

**Installation Test**

After installing MaCroDNA, you can test the installation using the following code in python

````
$ python
> from MaCroDNA import MaCroDNA
> test_macrodna = MaCroDNA()
> test_macrodna.tiny_test()
````
And the output should be like

````
******Test DNA data is:
      cell1  cell2  cell3  cell4
gene                            
g1        2      2      1      2
g2        2      2      1      2
g3        3      2      2      2
g4        1      2      2      2
g5        6      2      2      2
g6        2      2      3      6
******Test RNA data is:
      cell1  cell2  cell3  cell4
gene                            
g1        0      2      0      1
g2        0      2      0      1
g3       10      2      2      1
g4        0      2      2      1
g5       20      2      0      1
g6        0      2      5     20
g7        0      0      0      0
******Clone id for each DNA cell is:
   clone   cell
0      0  cell1
1      1  cell2
2      2  cell3
3      3  cell4
**********
Start Mapping RNA cells to DNA clones
**********
After selecting the same genes in RNA and DNA
number of cells in dna data 4
number of cells in rna data 4
number of genes in dna data 6
number of genes in rna data 6
[[ 0.95553309  0.          0.15877684 -0.29277002]
 [ 0.          0.          0.          0.        ]
 [-0.34698896  0.          0.87447463  0.86824314]
 [-0.18650096  0.          0.7592566   1.        ]]
1 0
MaCroDNA will be run for 1 steps
the smallest set has 4 number of cells
Set parameter Username
Academic license - for non-commercial use only - expires 2023-04-18
Gurobi Optimizer version 9.5.1 build v9.5.1rc2 (linux64)
Thread count: 10 physical cores, 20 logical processors, using up to 20 threads
Optimize a model with 9 rows, 16 columns and 48 nonzeros
Model fingerprint: 0xbe36aead
Variable types: 0 continuous, 16 integer (16 binary)
Coefficient statistics:
  Matrix range     [1e+00, 1e+00]
  Objective range  [2e-01, 1e+00]
  Bounds range     [1e+00, 1e+00]
  RHS range        [1e+00, 4e+00]
Found heuristic solution: objective 2.8300077
Presolve removed 9 rows and 16 columns
Presolve time: 0.00s
Presolve: All rows and columns removed

Explored 0 nodes (0 simplex iterations) in 0.00 seconds (0.00 work units)
Thread count was 1 (of 20 available processors)

Solution count 1: 2.83001 

Optimal solution found (tolerance 1.00e-04)
Best objective 2.830007718086e+00, best bound 2.830007718086e+00, gap 0.0000%
IsMIP: 1
Solved with MIPFocus: 0
The model has been optimized
Obj: 2.83001
the number of associations in the correspondence matrix 4.0
**********
Finish Mapping
Test Success
**********
````

For more complicated test, you can use `test_macrodna.py` under `test/` directory

### Usage

**cell-to-cell mapping**

For mapping RNA cells to DNA cells, the input should be two dataframes `dna_df` and 
`rna_df`. In these two dataframes, the row ids should be the genes and the column ids 
should be the cells. The genes on RNA data and those in DNA data can be different, because MaCroDNA can select the same genes
 from the data by itself.

To get the cell-to-cell mapping

````
$ python
> from MaCroDNA import MaCroDNA
> macrodna = MaCroDNA(rna_df, dna_df)
> cell2cell = macrodna.cell2cell_assignment()
````

The output `cell2cell` is also a dataframe.
The index ids are the RNA cell ids. 
And it only has on column "predict_cell", which is the DNA cell assigned to the corresponding RNA cell.


**cell-to-clone mapping**

For mapping RNA cells to DNA clones, the input needs three dataframes `dna_df`,
`rna_df` and `dna_label`. `dna_df` and `rna_df` are same as cell-to-cell mapping.
`dna_label` has two columns "clone" and "cell", where "cell" is the DNA cells and "clone" is the corresponding clone id
for that DNA cell.

To get the cell-to-clone mapping

````
$ python
> from MaCroDNA import MaCroDNA
> macrodna = MaCroDNA(rna_df, dna_df, dna_label)
> cell2clone = macrodna.cell2clone_assignment()
````

The output `cell2clone` is also a dataframe.
The index ids are the RNA cell ids. 
It has two columns. One is "predict_cell", which is the corresponding DNA cell for that RNA cell.
The other column is "predict_clone", which is the predict clone id for that RNA cell.

### Contacts
If you have any questions, please contact us via edrisi@rice.edu or xiru.huang@rice.edu
