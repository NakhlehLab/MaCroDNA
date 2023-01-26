import pandas as pd
import numpy as np 
import os
import shutil

def gene_info(x):
# Extract gene names, gene_type, gene_status and level
	g_name = list(filter(lambda x: 'gene_name' in x,  x.split(";")))[0].split("=")[1]
	g_type = list(filter(lambda x: 'gene_type' in x,  x.split(";")))[0].split("=")[1]
	g_status = list(filter(lambda x: 'gene_status' in x,  x.split(";")))[0].split("=")[1]
	g_level = int(list(filter(lambda x: 'level' in x,  x.split(";")))[0].split("=")[1])
	return (g_name, g_type, g_status, g_level)

if __name__=="__main__":
	gff_f = "./gencode.v19.annotation.gff3"

	# Read the gff3 file into a pandas dataframe
	gencode = pd.read_table(gff_f, comment="#", sep = "\t", names = ['seqname', 'source', 'feature', 'start' , 'end', 'score', 'strand', 'frame', 'attribute'])
	print(gencode.head())
	print(gencode.info())

	# Extract Genes in the gff3 file “feature = gene”
	gencode_genes = gencode[(gencode.feature == "gene")][['seqname', 'start', 'end', 'attribute']].copy().reset_index().drop('index', axis=1)
	# Extract gene_name, gene_type, gene_status, level of each gene
	gencode_genes["gene_name"], gencode_genes["gene_type"], gencode_genes["gene_status"], gencode_genes["gene_level"] = zip(*gencode_genes.attribute.apply(lambda x: gene_info(x)))
	# Extract all known protein_coding genes
	gencode_genes = gencode_genes[gencode_genes['gene_status'] == 'KNOWN'].reset_index().drop('index', axis=1)
	gencode_genes = gencode_genes[gencode_genes['gene_type'] == 'protein_coding'].reset_index().drop('index', axis=1)

	print(gencode_genes)

	print(set(gencode_genes['gene_status'].values.tolist()))

	# # Remove duplicates of the genes with multiple levels 
	gencode_genes = gencode_genes.sort_values(['gene_level', 'seqname'], ascending=True).drop_duplicates('gene_name', keep='first').reset_index().drop('index', axis=1)

	gencode_genes.to_csv("./gene_table_gencode.v19.annotation.gff3.csv")