---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: '0.8'
    jupytext_version: 1.4.2
kernelspec:
  display_name: Python 3
  name: python3
---

# TransposonRead_Profile_Compare 

+++

## Introduction
One of the most important visual results of the analysis of SATAY data is a plot that show the number of transposon insertions and the number of reads per location in the genome.
This script visualizes this using a bar plot for each individual chromosome.
These plots can also be useful for visualy comparing two datasets with each other as these bar plots can show deviations between two datasets.

This notebook includes two functions that are closely related to the functions in TransposonRead_Profile_Plot.py, but is altered to automatically show two datasets with each other and show the differences between the two.
One function is to show the location of transposon insertions and the second function shows the number of reads for each transposon insertion.

This can be used as a initial visual check for large differences between two datasets in terms of transposon insertion locations and number of reads per transposon.

Below are the packages that needs to be loaded for both functions.

```{code-cell} ipython3
import os, sys
import numpy as np
import matplotlib.pyplot as plt

file_dirname = os.path.dirname(os.path.abspath('__file__'))
# sys.path.insert(1,os.path.join(file_dirname,'python_modules'))
from Python_scripts.python_modules.chromosome_and_gene_positions import chromosome_position, chromosomename_roman_to_arabic, gene_position
from Python_scripts.python_modules.essential_genes_names import list_known_essentials
from Python_scripts.python_modules.gene_names import gene_aliases
from Python_scripts.python_modules.chromosome_names_in_files import chromosome_name_bedfile, chromosome_name_wigfile
```

## The python script for comparing transposon insertion locations

This function creates a bar plot with the number of transposons in a chromosome. The background of the bar plot has colored regions that indicate the locations of the genes. Green regions indicate annotated essential genes and red regions are the genes not annotated as essential. The names of the essential genes are given as well.

This graph can be used to check how well the transposon data corresponds with the location of essential genes, where only few transposon insertions are expected, which can be directly compared between two datasets.

+++

### Input
The function inputs the path to a bed file (`bed_file`, type=string) which is created by the Matlab code provided from the Kornman lab [https://sites.google.com/site/satayusers/complete-protocol/bioinformatics-analysis/matlab-script].
Next it inputs the chromosome number given as a roman numeral (`chrom_user_set`, type=string, list of strings or None-type, default=`None`).
When `chrom_user_set` is not given or set a None, the scripts automatically considers all chromosomes.
When a list of chromosomes is given, it considers all chromosomes given in the list.
The bar width that indicates the width of the final bar plot (`bar_width`, type=int, default=chromosome length divided by 500).
The `bar_width` variable takes as default value the length of the chromosome and divides this N bins or equal length.
When both `savefigure_path` (type=string, default=`None`) and `savefigure_name` (type=string, default=`None`) are given, the scripts automatically saves all figures to the given location and adds the extension `_chromX` where `X` is changed to the respective chromosome number.

If either one of the two are not given or set to None, the figures won't be saved.

```{code-cell} ipython3
bed_files = [os.path.join(file_dirname,'satay_analysis_testdata','Output_Processing','Cerevisiae_WT2_Michel2017_trimmed1.bam.bed'),
os.path.join(file_dirname,'satay_analysis_testdata','Output_Processing','Cerevisiae_WT2_Michel2017_trimmed2.bam.bed')] # CHANGE THIS TO ANY .BED FILES YOU WANT TO ANALYSE.
chrom_user_set = 'I'
bar_width_user_set=None
savefigure_path=None
savefigure_name=None
```

### Loading additional files
Next additional files are loaded. Change these to your local paths leading to the gff-file (for example downloaded from SGD [https://www.yeastgenome.org/] or get a copy from the docs folder on Github [https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis]) and yeast_protein_Names file.
Also two lists of essential genes are loaded. Some essential genes are present only in a single file, hence both files are used simultaneously.

```{code-cell} ipython3
gff_file = os.path.join(file_dirname,'Python_scripts','Data_Files','Saccharomyces_cerevisiae.R64-1-1.99.gff3')

essential_genes_files = [os.path.join(file_dirname,'Python_scripts','Data_Files','Cervisiae_EssentialGenes_List_1.txt'),
                        os.path.join(file_dirname,'Python_scripts','Data_Files','Cervisiae_EssentialGenes_List_2.txt')]
gene_information_file = os.path.join(file_dirname,'Python_scripts','Data_Files','Yeast_Protein_Names.txt')

```

### Get chromosome length and essential genes

Determine the length and position of all chromosomes and get a list of all genes, all essential genes and the aliases of all the genes.

```{code-cell} ipython3
# chr_length_dict, chr_start_pos_dict, chr_end_pos_dict = chromosome_position(gff_file)

# gene_pos_dict = gene_position(gff_file)
genes_essential_list = list_known_essentials(essential_genes_files)
gene_alias_list = gene_aliases(gene_information_file)[0]
```

