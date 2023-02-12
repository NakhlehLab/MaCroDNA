# Instructions for reproducing the results of our analysis on Barretts' esophagus data set originally introduced in:

> Busslinger, Georg A., et al. "Molecular characterization of Barrett’s esophagus at single-cell resolution." Proceedings of the National Academy of Sciences 118.47 (2021): e2113061118.


## Download and prepare the data set

All the data sets from this study are available at European Genome-Phenome Archive (https://ega-archive.org) and can be accessed by the accession number EGAS00001005221
The accession number for the single-cell DNA sequencing data is EGAD00001007521 and the accession number of the single-cell RNA sequencing data is EGAD00001007523
To map the reads of the scDNA-seq data to the human reference genome use the NIaIII mapping pipeline of SingleCellMultiOmics package: https://github.com/BuysDB/SingleCellMultiOmics/tree/master/singlecellmultiomics/snakemake_workflows/nlaIII
The reads of the scRNA-seq data were mapped to the human genome using the SingleCellMultiOmics pipeline:
https://github.com/BuysDB/SingleCellMultiOmics/tree/master/singlecellmultiomics/snakemake_workflows/cs2_scmo
From the outputs of these pipelines we use the read count tables and gene expression tables of the scDNA-seq and scRNA-seq data, respectively.


## Installing the required packages

**Package requirements:**
Python (tested on version 3.7.7)
- numpy (tested on version 1.21.5)
- pandas (tested on version 1.3.5)
- scipy (tested on version 1.7.3)
- math 
- gurobipy (tested on version 9.5.1)

To install gurobipy, one needs to install Gurobi optimizer. We used the free academic license of Gurobi. The installation instructions for Linux, Windows, and MacOS are as follows:
| OS      | Instruction |
|:--------|:------------|
| Linux | [Gurobi Installation Guide: Linux](https://youtu.be/yNmeG6Wom1o) |
| Windows | [Gurobi Installation Guide: Windows](https://youtu.be/fQVxuWOiPpI) |
| macOS | [Gurobi Installation Guide: macOS](https://youtu.be/ZcL-NmckTxQ) |

After installation of Gurobi optimizer, `gurobipy` can be installed using `pip`:
`Python -m pip install gurobipy`

**Installation of required R libraries:** 

Install phylosignal using the command:
`devtools::install_github("fkeck/phylosignal")`

Install phangorn R package for phylogenetic inference and the dependencies:
```
install.packages("phangorn")
install.packages("ape")
install.packages(“phylobase”)
install.packages(“”adephylo)
```


## Filtering on the scDNA-seq read count tables

For each biopsy, we have a read count table where the columns are the cell ID’s and the rows match the number of mapped reads in genomic bins with size of 250,000 base pairs. We filtered and normalized the DNA data using the script named cna_filterer.py which reads each read count table, removes the low-quality cells with less than 3000 coverage, adds pseudo count of 1 to all values, estimates the rough copy number values by multiplying the counts by 2 and dividing by the median for each cell, and finally performs log1p normalization.
To run `cna_filterer.py`:

1. Change the path to the source directory containing the original read count tables (`src_dir`)
2. Change the path to the destination directory to store the filtered and normalized tables (`tgt_dir`)
3. Run: `python cna_filterer.py`

The code does not need any argument, and the complete list of biopsy names are provided in the code in a list. The results are saved in the `tgt_dir`. For each biopsy, there will be a filtered CSV file named with `<biopsy name>_filtered_normed_count_table.csv`

## Filtering on the scRNA-seq gene expression tables:

For each biopsy, we have a gene expression table where the columns are the cell ID’s and the rows correspond the number of transcripts mapped to each gene. We filtered and normalized the RNA data using the script named `rna_filterer.py` which reads each gene expression table, removes the low-quality cells with total number of transcripts less than 3000, keeps the genes that are expressed in at least one cell with at least three transcripts, adds the pseudo count of 1 to all values, performs RPM normalization, and finally log1p normalization.
To run `rna_filterer.py`: 
1. Change the path to the source directory containing the original gene expression tables (`src_dir`)
2. Change the path to the destination directory to store the filtered and normalized tables (`tgt_dir`)
3. Run: `python rna_filterer.py`
The code does not need any argument, and the complete list of biopsy names are provided in the code in a list. The results are saved in the `tgt_dir`. For each biopsy, there will be a filtered CSV file named `<biopsy name>_filtered_normed_count_table.csv`


## Annotation of the scDNA-seq data

Prior to running MaCroDNA on the scDNA-seq and scRNA-seq data, we annotated the genomic bins in the copy number data by the genes that lie within them. The gencode data was downloaded the GFF3 annotation file for GRCh37.p13 from https://www.gencodegenes.org/human/release_19.html. The GFF3 file is named `gencode.v19.annotation.gff3` in the current folder. The protein-coding genes were extracted from the gencode data following the commands in https://medium.com/intothegenomics/annotate-genes-and-genomic-coordinates-using-python-9259efa6ffc2 . The commands for extracting the protein-coding genes are written in `gencode_reader.py`
In the python script, `gencode_reader.py` set the path to the GFF3 file (`gff_f`) and run:

`python gencode_reader.py`

It generates a CSV file named `gene_table_gencode.v19.annotation.gff3.csv` in the same directory that contains all the protein-coding genes with their corresponding genomic coordinates.

Next, to annotate the genomic bins in the copy number data with their overlapping genes, run `cna_annotator.py`.  Set the variable `src_dir` to the path to the directory containing the filtered copy number data and set `tgt_dir` to the path to the target directory for storing the annotated data. The variable `gff_f` points to the CSV file of the protein-coding genes.

Run:
`python cna_annotator.py`

The outputs are written in the `tgt_dir`. For each biopsy, the annotated file is named with `<biopsy name>_annotated_filtered_normed_count_table.csv`.


## Running MaCroDNA:


The script for running MaCroDNA on the BE data set is provided in `macrodna_BE.py`. Before running this script, set the variables `dna_src_dir` and `rna_src_dir` to the source directory of the scDNA-seq and scRNA-seq data, respectively. Note, the scDNA-seq data used at this step are the annotated data from the previous steps. 

Run: 
`python macrodna_BE.py`

For each biopsy, the output can be found in the `tgt_dir` in a CSV file named as `<biopsy name>_cell2cell_assignment_indexed.csv`. This file contains the names of the paired cells from RNA and DNA, and also the MaCroDNA round that they were assigned to each other. 


## Post-processing on MaCroDNA’s results for phylogenetic analysis


We used `aggergate_macrodna.py` to prepare the inputs to phylosignal (https://cran.r-project.org/web/packages/phylosignal/index.html) R package. This script takes the MaCroDNA results from the previous section and outputs the copy number and gene expression profiles for the paired cells. 
To run `aggregate_macrodna.py`:

1. Set `dna_src_dir` to the path to the directory containing the filtered (from the section [Filtering on the scDNA-seq read count tables]).
2. Set `rna_src_dir` to the path to the directory containing the filtered gene expression data (from the section [Filtering on the scRNA-seq gene expression tables]).
3. Set `res_dir` to the path to the directory containing the MaCroDNA results from the previous section.
4. Set `tgt_dir` to the desired path for storing the results.
Run the following command:
`python aggregate_macrodna.py` 

The outputs are saved in `tgt_dir` directory. For each biopsy, the copy number profiles are saved into CSV files named as `<biopsy name>_dna.csv`, and the gene expression profiles of the paired cells are saved in CSV files named as `<biopsy name>_rna.csv`


## Running phylosignal R package 

The codes for generating the results of our phylogenetic signal analysis are in `phylosignal.R`.
This script reads the files from the previous step and calculates the K* values for all the genes in each biopsy.
To run this code:

1. Set the working directory in line 7
2. Set the name of the biopsy in line 8
3. Change the path to the target directory from the previous step in line 9
Run the R script using the following command:
`Rscript phylosignal.R`

The outputs are the K* values (written in`<biopsy name>_stat.csv`), the corresponding p-values (written in `<biopsy name>_pval.csv`) and the newick file of the phylogenetic tree (written in `<biopsy name>_tre.nwk`) in the working directory.

