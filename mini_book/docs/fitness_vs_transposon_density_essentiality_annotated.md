---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: '0.8'
    jupytext_version: 1.4.1+dev
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Fitness vs Transposon density 

**Is there a fitness difference between low density and high density?**


This script creates a scatterplot for transposon_density vs fitness-SGA.
All genes that are used to create the plot are searched for is they are annotated as essential genes.
This is also done for any potential aliases the genes might have.

This code is a part from 'fitness_vs_transposon_density.py' created by Leila.


```{code-cell} ipython3 

import numpy as np
import scipy.io
import seaborn as sns
from scipy import stats, optimize, interpolate
import pandas as pd
from collections import defaultdict 
import math
import matplotlib.pyplot as plt
from scipy.stats import norm, lognorm
from scipy import stats
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import os, fnmatch
import sys
```
```{code-cell} ipython3 

file_dirname = os.path.dirname('__file__')

from Python_scripts.python_modules.essential_genes_names import list_known_essentials
from Python_scripts.python_modules.gene_names import gene_aliases

```
## Importing datasets 
### Published data in the cellmap.org 2016 for the fitness data

```{code-cell} ipython3 


transposon_density_file=os.path.join(file_dirname,'Python_scripts','Data_Files','transposon-density-per-gene-benoit.xlsx')
essential_list_1= os.path.join(file_dirname,'Python_scripts','Data_Files','Cervisiae_EssentialGenes_List_1.txt')
essential_list_2= os.path.join(file_dirname,'Python_scripts','Data_Files','Cervisiae_EssentialGenes_List_2.txt')
yeast_protein_names_list=os.path.join(file_dirname,'Python_scripts','Data_Files','Yeast_Protein_Names.txt')



data_transposon=pd.read_excel(transposon_density_file,header=0)
data_transposon=data_transposon.drop(['Unnamed: 0'],axis=1)
data_transposon=data_transposon.apply(lambda x: x.astype(str).str.lower()) # make everything lowercase

```



## Transform back these columns into float variables

```{code-cell} ipython3
data_transposon['Transposon_density_per_gene']=data_transposon['Transposon_density_per_gene'].apply(lambda x: (float(x)))
data_transposon['Read_density_per_gene']=data_transposon['Read_density_per_gene'].apply(lambda x: (float(x)))
```



## How transposon density varies with the reads per gene 

```{code-cell} ipython3
sns.pairplot(data=data_transposon,hue='Essential_gene',diag_kind='kde',kind='reg',palette='colorblind')
```

## USE FUNCTIONS TO GET NAMES OF ESSENTIAL GENES. FILES USED ARE:

```{code-cell} ipython3

known_essential_gene_list = list_known_essentials([essential_list_1,essential_list_2])
aliases_designation_dict, aliases_sgd_dict, aliases_swissprot_dict = gene_aliases(yeast_protein_names_list)

known_essential_gene_list
```


