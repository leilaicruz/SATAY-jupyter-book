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

# Comparison of the statistics per chromosome

+++

## Introduction

This script gives a numerical overview of the two datasets in the form of a textfile.
This can be used for comparing statistical values of two datasets.

+++

## The python script
The script inputs two .bed files in a list that were generated using the Matlab code provided by the Kornmann lab and the path and name of a text file where the results will be written.
Currently 11 different statistical values are determined, but this can be easily extended (see explainatory text between the codes below).

1. Number of transposon insertions.
2. Percentage of the chromosome that is covered by transposons
3. Mean distance between transposon insertions.
4. Median distance between transposon insertions.
5. 25th percentile of the distance between transposon insertions.
6. 75th percentile of the distance between transposon insertions.
7. Largest area devoid of transposons
8. Mean number of reads per transposon
9. Median numbr of reads per transposon
10. 25th percentile reads per transposon
11. 75th percentile reads per transposon


+++

### Input
The function inputs either a gene name (`gene_name`, type=string) or a region (`region`, type=list) specified as a list with three entries (chromosome number as a roman numeral, start position of the region and the end position respectively).
The variable `gene_name` can be set to any gene name or `holocus` or `ho-locus`.
Next it requires the bed-file (`bed_file`, type=string) which is created by the Matlab code provided from [the Kornman lab](https://sites.google.com/site/satayusers/complete-protocol/bioinformatics-analysis/matlab-script).
Finally, the figure can be automatically saved (at a location specified in the beginning of the function) by setting `savefigure` to `True`.

The custom build functions (stored in the ['python modules' folder on Github](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis/tree/master/python_modules)) that are required are:
- chromosome_and_gene_positions.chromosome_position
- chromosome_and_gene_positions.chromosome_roman_to_arabic
- chromosome_names_in_fies.chromosome_name_bedfile

```{code-cell} ipython3
import os, sys
import numpy as np
from tabulate import tabulate

file_dirname = os.path.dirname('__file__')
# sys.path.insert(1,os.path.join(file_dirname,'python_modules'))
from Python_scripts.python_modules.chromosome_and_gene_positions import chromosome_position, chromosomename_roman_to_arabic
from Python_scripts.python_modules.chromosome_names_in_files import chromosome_name_bedfile


#FUNCTION INPUTS
bed_files = [os.path.join(file_dirname,'satay_analysis_testdata','Output_Processing','Cerevisiae_WT2_Michel2017_trimmed1.bam.bed'),
            os.path.join(file_dirname,'satay_analysis_testdata','Output_Processing','Cerevisiae_WT2_Michel2017_trimmed2.bam.bed')]
text_file = os.path.join(file_dirname,'satay_analysis_testdata','Output_Processing','Cerevisiae_WT2_Michel2017_trimmed1-2_Compare.txt')

for test_path in bed_files:
    print(os.path.isfile(test_path))
```

### Loading files
Next additional files are loaded. Change this to your local paths leading to the gff-file (for example downloaded from SGD [https://www.yeastgenome.org/] or get a copy from the docs folder on Github [https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis]).

```{code-cell} ipython3
gff_file = os.path.join(file_dirname,'Python_scripts','Data_Files','Saccharomyces_cerevisiae.R64-1-1.99.gff3')
```

### Set name of text file and create header

The path and name that are given as input for saving the text file is extended with 'DataCompare.txt' to make sure the file is recognizable and have the right extension.

Secondly, the header for the text file is created here.

```{code-cell} ipython3
if text_file == True:
    text_file = os.path.splitext(bed_files[0])[0] + '_DataCompare.txt'

if text_file is not None:
    t = open(text_file, 'w+')
    if len(bed_files) == 1:
        t.write('Statistical values for datafile: ' + bed_files[0].split('\\')[-1] + '\n')
    elif len(bed_files) == 2:
        t.write('Statistical values for datafile: ' + bed_files[0].split('\\')[-1] + ' and ' + bed_files[1].split('\\')[-1] + '\n')
    t.close()
```

### Get chromosome information

Determine the lengths and position of the different chromosomes and get a list of the chromosome names in terms of roman numerals.

```{code-cell} ipython3
chr_length_dict, chr_start_pos_dict, chr_end_pos_dict = chromosome_position(gff_file)

roman_to_arabic_dict = chromosomename_roman_to_arabic()[1]
chromosome_romannames_list = []
for roman in roman_to_arabic_dict:
    chromosome_romannames_list.append(roman)
```

### Create lists for the variables

All the variables that are determined are put in individual lists.
When new statistics needs to be determined, add a new list for each new value.

```{code-cell} ipython3
N_Tn_Ins = []
Cov_Percent = []
Mean_Tn_Dist = []
Median_Tn_Dist = []
Tn_Dist_25th_Percent = []
Tn_Dist_75th_Percent = []
Tn_Gap = []
Mean_Reads_PerTn = []
Median_Reads_PerTn = []
Reads_25th_percent = []
Reads_75th_percent = []
```

### Determine statistics

For each chromosome in each .bed file the statistical values are determined and stored in the respective lists.
The values are also determined for the entire genome.

First a for loop is done over all bed files and for each bed file the names of the chromosomes as they are stored in the bed files are determined.
Next, some statistical values are determined and depending whether it is first file written a text file is generated or the values are appended to the already exsiting text file.
After all chromosomes the statistics are determined for the entire genome as well.

```{code-cell} ipython3
bed_file_counter = 0
for bed_file in bed_files:
    with open(bed_file) as f:
        lines = f.readlines()

    chrom_names_dict, chrom_start_index_dict, chrom_end_index_dict = chromosome_name_bedfile(lines)
    chrom_loop = chrom_names_dict

#DETERMINE STATISTICS
    bp_between_tn_insertions_dict = {}
    reads_per_tn_dict = {}
    for chrom in chrom_loop:
        tn_insertion_position_list = []
        reads_per_tn_list = []
        for line in lines[chrom_start_index_dict.get(chrom):chrom_end_index_dict.get(chrom)+1]:
            line = line.strip('\n').split()
            tn_insertion_position_list.append(int(line[1]))
            reads_per_tn_list.append((int(line[4])-100)/20)
        bp_between_tn_insertions = [abs(y-x) for x, y in zip(tn_insertion_position_list[:-1], tn_insertion_position_list[1:])]
        bp_between_tn_insertions.insert(0,tn_insertion_position_list[0]) #ADD START OF GENE (bp=0)
        bp_between_tn_insertions.append(chr_length_dict.get(chrom) - tn_insertion_position_list[-1]) #ADD END OF GENE (bp=GENE_LENGTH-TRANSPOSON INSERTION LAST GENE)
        bp_between_tn_insertions_dict[chrom] = bp_between_tn_insertions
        reads_per_tn_dict[chrom] = reads_per_tn_list

        tn_insertion_meanfrequency = np.nanmean(bp_between_tn_insertions)
        tn_insertion_25percentilefrequency = np.percentile(bp_between_tn_insertions,25)
        tn_insertion_medianfrequency = np.nanmedian(bp_between_tn_insertions)
        tn_insertion_75percentilefrequency = np.percentile(bp_between_tn_insertions,75)

#IF IT IS THE FIRST BED FILE, CREATE TITLES FOR EACH STATISTIC TO WRITE IN THE TEXT FILE.
        if bed_file_counter == 0:
            print('Print information chromosome ' + chrom + ' with length ' + str(chr_length_dict.get(chrom)))
            N_Tn_Ins.append([chrom, 'Number of transposon insertions', len(reads_per_tn_list), ''])
            Cov_Percent.append([chrom, 'Coverage percent', len(tn_insertion_position_list)/chr_length_dict.get(chrom)*100, ''])
            
            Mean_Tn_Dist.append([chrom, 'Mean distance between insertions', tn_insertion_meanfrequency, ''])
            Median_Tn_Dist.append([chrom, 'Median distance between insertions', tn_insertion_medianfrequency, ''])
            Tn_Dist_25th_Percent.append([chrom, '25th percentile distance between insertions', tn_insertion_25percentilefrequency, ''])
            Tn_Dist_75th_Percent.append([chrom, '75th percentile distance between insertions', tn_insertion_75percentilefrequency, ''])
            
            Tn_Gap.append([chrom, 'Largest area devoid of transposons', max(bp_between_tn_insertions), ''])
            
            Mean_Reads_PerTn.append([chrom, 'Mean number of reads per transposon', np.nanmean(reads_per_tn_list), ''])
            Median_Reads_PerTn.append([chrom, 'Median number of reads per transposon', np.nanmedian(reads_per_tn_list), ''])
            Reads_25th_percent.append([chrom, '25th percentile reads per transposon', np.percentile(reads_per_tn_list,25), ''])
            Reads_75th_percent.append([chrom, '75th percentile reads per transposon', np.percentile(reads_per_tn_list,75), ''])

#IF IT IS NOT THE FIRST BED FILE, APPEND THE STATISTICS FOR THE CURRENT BED FILE TO THE EXISTING STATISTICS OF THE PREVIOUS BED FILE.
        elif bed_file_counter == 1:
            N_Tn_Ins[chromosome_romannames_list.index(chrom)][-1] = len(reads_per_tn_list)
            Cov_Percent[chromosome_romannames_list.index(chrom)][-1] = len(tn_insertion_position_list)/chr_length_dict.get(chrom)*100
            
            Mean_Tn_Dist[chromosome_romannames_list.index(chrom)][-1] = tn_insertion_meanfrequency
            Median_Tn_Dist[chromosome_romannames_list.index(chrom)][-1] = tn_insertion_medianfrequency
            Tn_Dist_25th_Percent[chromosome_romannames_list.index(chrom)][-1] = tn_insertion_25percentilefrequency
            Tn_Dist_75th_Percent[chromosome_romannames_list.index(chrom)][-1] = tn_insertion_75percentilefrequency
            
            Tn_Gap[chromosome_romannames_list.index(chrom)][-1] = max(bp_between_tn_insertions)
            
            Mean_Reads_PerTn[chromosome_romannames_list.index(chrom)][-1] = np.nanmean(reads_per_tn_list)
            Median_Reads_PerTn[chromosome_romannames_list.index(chrom)][-1] = np.nanmedian(reads_per_tn_list)
            Reads_25th_percent[chromosome_romannames_list.index(chrom)][-1] =  np.percentile(reads_per_tn_list,25)
            Reads_75th_percent[chromosome_romannames_list.index(chrom)][-1] = np.percentile(reads_per_tn_list,75)

#DETERMINE STATISTICS FOR THE ENTIRE GENOME
    bp_between_tn_insertions_genome = []
    number_tn_insertions_list = []
    reads_per_tn_genome = []
    number_tn_insertions_genome = 0
    for chrom in chrom_loop:
        number_tn_insertions_genome += len(reads_per_tn_dict.get(chrom))
        for bp_between in bp_between_tn_insertions_dict.get(chrom):
            bp_between_tn_insertions_genome.append(bp_between)
        number_tn_insertions_list.append(len(bp_between_tn_insertions_dict.get(chrom)))
        for reads_tn in reads_per_tn_dict.get(chrom):
            reads_per_tn_genome.append(reads_tn)

    if bed_file_counter == 0:
        N_Tn_Ins.append(['Genome','Number of insertions',number_tn_insertions_genome, ''])
        Cov_Percent.append(['Genome', 'Coverage percent', sum(number_tn_insertions_list)/sum(chr_length_dict.values())*100, ''])
        Mean_Tn_Dist.append(['Genome', 'Mean distance between insertions', np.nanmean(bp_between_tn_insertions_genome), ''])
        Median_Tn_Dist.append(['Genome', 'Median distance between insertions', np.nanmedian(bp_between_tn_insertions_genome), ''])
        Tn_Dist_25th_Percent.append(['Genome', '25th percentile distance between insertions', np.percentile(bp_between_tn_insertions_genome,25), ''])
        Tn_Dist_75th_Percent.append(['Genome', '75th percentile distance between insertions', np.percentile(bp_between_tn_insertions_genome,75), ''])
        Tn_Gap.append(['','','', ''])
        Mean_Reads_PerTn.append(['Genome', 'Mean number of reads per transposon', np.nanmean(reads_per_tn_genome), ''])
        Median_Reads_PerTn.append(['Genome', 'Median number of reads per transposon', np.nanmedian(reads_per_tn_genome), ''])
        Reads_25th_percent.append(['Genome', '25th percentile reads per transposon', np.percentile(reads_per_tn_genome,25), ''])
        Reads_75th_percent.append(['Genome', '75th percentile reads per transposon', np.percentile(reads_per_tn_genome,75), ''])

    elif bed_file_counter == 1:
        N_Tn_Ins[-1][-1] = number_tn_insertions_genome
        Cov_Percent[-1][-1] = sum(number_tn_insertions_list)/sum(chr_length_dict.values())*100
        Mean_Tn_Dist[-1][-1] = np.nanmean(bp_between_tn_insertions_genome)
        Median_Tn_Dist[-1][-1] = np.nanmedian(bp_between_tn_insertions_genome)
        Tn_Dist_25th_Percent[-1][-1] = np.percentile(bp_between_tn_insertions_genome,25)
        Tn_Dist_75th_Percent[-1][-1] = np.percentile(bp_between_tn_insertions_genome,75)
        Tn_Gap[-1][-1] = ''
        Mean_Reads_PerTn[-1][-1] = np.nanmean(reads_per_tn_genome)
        Median_Reads_PerTn[-1][-1] = np.nanmedian(reads_per_tn_genome)
        Reads_25th_percent[-1][-1] = np.percentile(reads_per_tn_genome,25)
        Reads_75th_percent[-1][-1] = np.percentile(reads_per_tn_genome,75)

    bed_file_counter += 1
```

### Creating text file

Write the stored statistical value to the text file.

```{code-cell} ipython3
print('Writing to text file...')

header0 = ['chromosome','item','Dataset 1','Dataset 2']
header  = ['          ','    ','         ','         ']

with open(text_file,'a') as t:
    for i in range(0,len(N_Tn_Ins)):
        table = [N_Tn_Ins[i], Cov_Percent[i], Mean_Tn_Dist[i], Median_Tn_Dist[i], Tn_Dist_25th_Percent[i],
                    Tn_Dist_75th_Percent[i], Tn_Gap[i], Mean_Reads_PerTn[i], Median_Reads_PerTn[i],
                    Reads_25th_percent[i], Reads_75th_percent[i]]
        if i == 0:
            t.write(tabulate(table,tablefmt='github',headers=header0))
        else:
            t.write(tabulate(table,tablefmt='github',headers=header))
        t.write('\n')
    print('...Writing complete.')
```

### Printing result

Showing part of the text file to give an example of the created text file.
(This is just for showing the result in this notebook and is not important for the code itself.
Can be removed).

```{code-cell} ipython3
:tags: [outputPrepend]

with open(text_file) as t:
    for i in range(0,20):
        print(t.readline())
```

## Interpretation

The statistical values that are determined here for two datasets give an indication of some of the properties of the datasets.
Together with the TransposonRead_Profile_Compare.py script, this can be helpful when comparing two datasets with each other to possibly improve the preprocessing steps.
This script can relatively easy be extended with more statistical values in the future.

The `Coverage percentage` is the number of transposons divided by the number of basepairs of the chromosome.
The distance between transposon insertions is determined by the taking the absolute difference between all subsequent transposon and for the first and last transposon the distance is determined from the beginning and until the end of the chromosome, respectively.
The same goes for the median and the percentiles.
The largest distance between subsequent transposons is displayed as the largest area devoid of transposons.
The values regarding the number of reads are directly extracted from the bed files.

+++

## Bibliography
- Michel, A. H., Hatakeyama, R., Kimmig, P., Arter, M., Peter, M., Matos, J., ... & Kornmann, B. (2017). Functional mapping of yeast genomes by saturated transposition. Elife, 6, e23570.
