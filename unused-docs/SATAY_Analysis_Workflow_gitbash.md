---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: '0.8'
    jupytext_version: 1.4.2
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Workflow SATAY Analysis

This is an example on how the workflow for SATAY analysis can be done using the commandline.
A full description of the different programs and their options can be found at M/tnw/bn/ll/Shared/Gregory/Notes/ProjectNotes.docx.
This workflow is based on Windows.
However, some of the software is not available for Windows, so a virtual machine is created where a part of the workflow is completed.
If wanted, it should be possible to complete the whole workflow in Linux (without using windows).
In this case follow the same workflow and commands, just ignore everything that has to do with copying data to or from a shared folder.

The command line arguments below are written to be used in git-bash (or similar command line tool).
Running the same commands in the standard Windows command line can be troublesome.

## Overview
**Performed on Windows:**
Sequencing results in FASTQ files that contains all the reads in random order.
This needs to be checked and trimmed to remove low quality reads, adapter sequences etc.

**Performed on Linux Virtual Machine:**
After trimming the data needs to be aligned to a reference sequence.
If using the reference sequence for the first time, it is necessary to index this sequence.
This results in a SAM file that includes all the reads including where they are mapped to the reference sequence and the mapping quality.
The SAM file needs to be converted to its binary equivalent BAM file which is then index as well.

**Performed on Windows:**
After alignment the BAM file can be used for visualization using IGV.
Also this BAM file is the input for the Matlab code (created by the Kornmann lab) for counting the number of reads and transposon insertions per gene.


## Tips on using this notebook
The Initialization sections contain variables that need to be set manually.
This includes names of the files to be used, paths to these files and paths to directories where the results need to be stored.
Once this is set, the other commands should be able to simply be copy-pasted in the terminal and they should work.

For detailed information about the different programs, look at the SATAY_Analysis_notes.docx file.

+++

## Initializing
For easier use, variables can be used for setting paths and file names.
To set this up, enter the following command:

```{code-cell} ipython3
#!/bin/bash
```

Define all paths and file names.
These are

1. pathdata: Path to the FASTQ file where also all results will be stored.
2. filename: Name of the FASTQ file to be analyzed. This can be left uncompressed (i.e. with the extension .fastq.gz)

3. path_sharedfolder: Use this to copy the needed files to a shared folder for the virtual machine (leave empty or type 'None' if not needed).

The next two commands are for calling the fastq and trimmomatic programs. Enter here the paths as you would use to call the programs in the command line.

4. path_fastqc: Path to the location where the Fastqc program is saved.
5. path_trimmomatic: Path to the location where the Trimmomatic program is saved.

The last two commands are for the file containing the reference sequence.

6. Path_reference: Path to the reference sequence, typically stored as a fasta file.
7. file_reference: Name of the file containing the reference sequence.
This is needed for copying the refence to the shared folder for using it in the Virtual Linux Machine.
If a reference is already stored in the Virtual Linux Machine, than this variable can be left empty.

```{code-cell} ipython3
import os
script_dir = os.path.dirname('__file__') #<-- absolute dir the script is in
rel_path_fastq = "datasets/SRR062634.filt.fastq.gz"

# pathname='N:\tnw\BN\LL\Shared\Gregory\Sequence_Alignment_TestData\temp_test_Michel_WT1\'
# filename='SRR062634.filt.fastq.gz'
# path_sharedfolder="/C/Users/gregoryvanbeek/VirtualBox VMs/VM_SharedFolder_Ubuntu64_Large/"
# path_fastqc='/C/Programs/FastQC/'
# path_trimmomatic='/C/Programs/Trimmomatic-0.39/'
# path_reference='N:\tnw\BN\LL\Shared\Gregory\Reference_Genome_Cerevisiae\Cerevisiae_S288C_RefR64-2-1\'
rel_path_fsa='datasets/S288C_reference_sequence_R64-2-1_20150113.fsa'

abs_file_fastq= os.path.join(script_dir, rel_path_fastq)
abs_file_fsa= os.path.join(script_dir, rel_path_fsa)

os.chdir('mini_book/docs/') #<-- for binder os.chdir('../')
my_path_fastq= abs_file_fastq
my_path_fsa=abs_file_fsa
```

## Make directories

Create all the necessary directories within the directory specified in the 'pathdata' variable.
To make folders with clear names, the extension is removed from the file name.
If the file ends with something else than .gz and/or .fastq, change the commands below accordingly.

```{code-cell} ipython3
filename_without_extension=${filename//.gz}
filename_without_extension=${filename_without_extension//.fastq}

mkdir ${pathname}${filename_without_extension}_QC
path_qc=${pathname}${filename_without_extension}'_QC\'
mkdir ${pathname}${filename_without_extension}_Trimmed
path_trimmed=${pathname}${filename_without_extension}'_Trimmed\'
mkdir ${pathname}${filename_without_extension}_Aligned
path_aligned=${pathname}${filename_without_extension}'_Aligned\'
```

## Backup and uncompress data

Make a backup of the fastq file and if it needs to be unzipped.

```{code-cell} ipython3
filepath=${pathname}${filename}
filepath_backup=${pathname}'Backup_'${filename}


cp ${filepath} ${filepath_backup}
chmod -w ${filepath_backup}


if [ ${filepath: -3} == '.gz' ]; then gunzip ${filepath}; fi
filepath=${pathname}${filename//.gz}
```

## Fastqc

Open Fastqc for checking the quality of the sequencing results.
The input is the fastq file and output is stored in the _QC-folder.
(This can take a few minutes to run).
Use the 'start' command to checking the quality report.

```{code-cell} ipython3
${path_fastqc}'fastqc' --outdir ${path_qc} ${filepath} --extract

start ${path_qc}${filename_without_extension}'_fastqc.html'
```

## Trimmomatic

Open Trimmomatic to trim the sequencing results.
The input the fastq file and the output is a fastq file with the trimmed reads.
For trimming, adapter sequences might be needed that have to be stored at the same location where the fastq file is stored.
Some standard adapter sequences are stored in the 'adapter' directory within of Trimmomatic.
These can be copied to the directory where the fastq file is stored.
If wanted, custom adapter sequences can be stored in the form of a .fasta file.
In this case 'adapter_file' might need to be changed.

```{code-cell} ipython3
adapter_file='TruSeq3-SE.fa'
cp ${path_trimmomatic}'adapters\'${adapter_file} ${pathname}${adapter_file}


java -jar ${path_trimmomatic}'trimmomatic.jar' SE -phred33 ${filepath} ${path_trimmed}${filename_without_extension}'_Trimmed.fastq' ILLUMINACLIP:${adapter_file}:2:15:30:8 SLIDINGWINDOW:5:20 TRAILING:10
```

## Fastqc on trimmed data

Redo the Fastqc quality check on the trimmed data to see if the quality improved as expected.
The results are stored in the _QC folder.

```{code-cell} ipython3
${path_fastqc}'fastqc' --outdir ${path_qc} ${path_trimmed}${filename_without_extension}'_Trimmed.fastq' --extract
trimmed_filename=${filename_without_extension}'_Trimmed.fastq'

start ${path_qc}${filename_without_extension}'_Trimmed_fastqc.html'
```

## Copy results to shared folder

The results are copied to the shared folder for continuing in the Virtual Linux Environment.

```{code-cell} ipython3
cp ${path_trimmed}${trimmed_filename} "${path_sharedfolder}"${trimmed_filename}

cp ${path_reference}${file_reference} "${path_sharedfolder}"${file_reference}
```

## Initializing and setting up Linux terminal

Since now the process continues on a Linux machine, the terminal on the Linux machine needs to be initialized, similarly to what is done in Windows.

```{code-cell} ipython3
#!/bin/bash
path_sharedfolder_vm='/media/sf_VM_sharedFolder_Ubuntu64_Large/'

trimmed_filename_vm='SRR062634.filt_Trimmed.fastq'

path_reference_vm='~/Documents/Reference_Genomes/Reference_Sequence_CerevisiaeS288C/'
file_reference_vm='S288C_reference_sequence_R64-2-1_20150113.fsa'

path_datafolder_vm = '~/Documents/Michel2017_WT1/'
```

## Index reference sequence

In case this is not already done before, create a dictionary for the reference sequence and index the reference.
This only has to be done once for each reference sequence.

```{code-cell} ipython3
mkdir ${path_refence_vm}

cp ${path_sharedfolder_vm}${file_reference_vm} ${path_reference_vm}${file_reference_vm}

bwa index ${path_reference_vm}${file_reference_vm}
```

## Perform alignment

Align the sequence with respect to the reference sequence.
The aligned file needs to have the '.sam' extension.
The inputs of the bwa mem funtion are the scores and penalties used for aligning the reads.
These are respectively
* A: matching score (default=1)
* B: mismatch score (default=4)
* O: gap open penalty (default=6)
* E: gap extension penalty (default=1)
* U: penalty for unpaired reads in case of paired end data (default=9)

```{code-cell} ipython3
mkdir ${path_datafolder_vm}
path_datafolder_vm_aligned=${path_datafolder_vm}'AlignmentOutput/'
mkdir path_datafolder_vm_aligned

cp ${path_sharedfolder_vm}${trimmed_filename_vm} ${path_datafolder_vm}${trimmed_filename_vm}

bwa mem -A 1 -B 4 -O 6 -E 1 ${path_reference_vm}${file_reference_vm} ${path_datafolder_vm}${trimmed_filename_vm} > ${path_datafolder_vm_aligned}'Aligned_'${trimmed_filename_vm}'.sam'
```

## Converting SAM to BAM files
The human readable SAM file needs to be converted to the binary BAM file for downstream analysis.
Also, the BAM file needs to be indexed like what was done for the reference sequence.

```{code-cell} ipython3
samtools view -b ${path_datafolder_vm_aligned}'Aligned_'${trimmed_filename_vm}'.sam' > ${path_datafolder_vm_aligned}'Aligned_'${trimmed_filename_vm}'.bam'

sambamba-0.7.1 sort -m 500MB ${path_datafolder_vm_aligned}'Aligned_'${trimmed_filename_vm}'.bam'
sambamba-0.7.1 index ${path_datafolder_vm_aligned}'Aligned_'${trimmed_filename_vm}'.bam'
```

## Copying results to shared folder for continuing on Windows machine
The results folder should now contain the following files:
* .sam
* .bam
* .sorted.bam
* .bam.bai

```{code-cell} ipython3
cp ${path_datafolder_vm_aligned}'Aligned_'${trimmed_filename_vm}'.*' ${path_sharedfolder_vm}
```

## Copying on Windows machine and make backup

This uses some of the variables that we initialized in the beginning of the workflow.

```{code-cell} ipython3
aligned_filenames='Aligned_'${trimmed_filename}
cp "${path_sharedfolder}"${aligned_filenames}'.*' ${path_aligned}

sam_filename=${aligned_filenames}'.sam'
cp ${path_aligned}${sam_filename} ${path_aligned}'Backup_'${sam_filename}
chmod -w ${path_aligned}'Backup_'${sam_filename}
```

## Visualizing .sorted.bam in IGV viewer
Open the IGV viewer and use the field ‘genomes’ > ‘Load genome from file’ to enter a reference sequence (e.g. .fasta). Use the field ‘File’ > ‘Load from file’ to enter a .sorted.bam file (the folder where this file is located should also contain the .bam.bai file). It might be necessary to zoom in quite a bit in order to the reads and any potential differences with respect to the reference sequence (indicated with letters). If a specific locus or gene needs to be visualized, then the header in the top of the screen can be used. An example of an input would be chr1:335-649. Note that the name of the chromosome can be different. This is typically given in the header of the .sam file (e.g. ref|NC_001133| instead of chr.1). To get the locations of known genes, use the .gff file that comes with downloading the reference sequence.

+++

## Transposon count using Matlab code Benoit
Open the Matlab code from Benoit and run the code. It will ask to open a .bam file. Load the aligned and sorted bam file. The program runs for quite a while and creates a number of output files that are stored in the same folder as where the bam file is located.
