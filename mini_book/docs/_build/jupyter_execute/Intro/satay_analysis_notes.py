# All what you need to know from SATAY

This file discusses the general outline of the experiments and interpretation of experiments using SAturated Transposon Analysis in Yeast (SATAY).
The introduction explains the purpose of the experiments and what kind of results we expect and how to interpret these results.
The Methods and File types section explains the different kind of files used during processing and analysis and how to read and use them.
It also explains some custom made software for analyzing the data.
This software is stored in the [Github page](https://github.com/Gregory94/LaanLab-SATAY-DataAnalysis) of the this research.
The steps needed for the processing and initial analysis of the data is explained in the Data Analysis section.
This includes a detailed step-by-step tutorial on how to use the different software packages and what steps need to be taken for the processing to go from the raw sequencing reads to the location of all the transposon insertion and the number of reads at each location.

```{toggle} To get the context of the technique, click the buttom to reveal :) 

About 20% of the genes in wild type *Saccharomyces Cerevisiae* are essential, meaning that they cannot be deleted without crippling the cell to such an extent that it either cannot survive (lethality) or multiply
(sterility).
(Non)-essentiality of genes is not constant over different genetic backgrounds, but genes can gain or lose essentiality when other genes are mutated.
Moreover, it is expected that the interactions between genes changes in mutants (changes in the interaction map).
This raises a number of questions:

- *If there exists a relation between genes which gain or lose essentiality after mutations and the changes in the interaction map?*

- If a gene x gains or loses essentiality after a mutation in gene y, does the essentiality of gene y also changes if a mutation in gene x is provoked?

- After a mutation that reduces the fitness of a population of cells, the population is sometimes able to increase its fitness by mutating other genes (e.g. dBem1 eventually result in mutations in Bem3).
Can these mutations, that are initiated by cells themselves, be predicted based on the interaction maps (i.e. predict survival of the fittest)?

- If a gene x is suppressed, it will possibly change the essentiality of another gene.
It is expected that most changes in essentiality will occur in the same subnetwork of the mutated gene.
If a gene y is suppressed that is part of the same network as gene x, does this invoke similar changes in this subnetwork?
In other words, are there common changes in the subnetwork when a random change is made within this subnetwork?

- Are there relations between the changes in the interaction network after a mutation and the Genetic Ontology (GO-)terms of changed genes?

To check the essentiality of genes, SATAY (SAturated Transposon Analysis in Yeast) experiments will be performed on different genetic backgrounds [Michel et.al., 2017] [Segal et.al., 2018].
This uses transposons to inhibit genes and it allows to compare the effects of this inhibition on the fitness of the cells (see for example the [galaxyproject website](https://galaxyproject.github.io/training-material/topics/genome-annotation/tutorials/tnseq/tutorial.html) which explains it in the context of bacteria, but the same principles hold for yeast cells).

- Transposons are small pieces of DNA that can integrate in a random location in the genome.
When the insertion happens at the location of a gene, this gene will be inhibited (i.e. it can still be transcribed, but typically it cannot be translated into a functional protein).

- After a transposon is randomly inserted in the DNA, the growth of the cell is checked.
- If the cell cannot produce offspring, the transposon has likely been inserted in an essential gene or genomic region.
- This is done for many cells at the same time.

All the cells are let to grow after insertion and the cells that have a transposon inserted in an essential part of their DNA (therefore having a low fitness) will, after some time, be outcompeted by the cells with an insertion in a non-essential part of the DNA (cells with a relatively high fitness).
By means of sequencing, the location of the transposon insertion can be checked and related to a specific gene.
Since the cells that have an essential part of the genome blocked do not occur in the population, those cells are not sequenced and hence the location of the transposon in these cells are missing in the sequencing results.
Thus, when the sequencing results of all the cells are aligned to a reference genome, some genomic locations are missing in the sequencing results.
These missing locations corresponds to potentially essential genomic regions.
The genome of all cells (called the library) is saturated when all possible insertion sites have at least one insertion of a transposon.

In that case all regions of the DNA are checked for their essentiality.
```


```{margin} A SATAY FORUM 
If you want to get started just join this [FORUM](https://groups.google.com/forum/#!forum/satayusers). 
And the website for the library generation can be found [HERE](https://sites.google.com/site/satayusers/complete-protocol/library-generation)
```


## Gene Essentiality

Essentiality of genes are defined as its deletion is detrimental to cell in the form that either the cell cannot grow anymore and dies, or the cell cannot give rise to offspring. Essentiality can be grouped in two categories, namely type I and type II [Chen et.al. 2016].
- Type I essential genes are genes, when inhibited, show a loss-of-function that can only be rescued (or masked) when the lost function is recovered by a gain-of-function mechanism.
- Typically these genes are important for some indispensable core function in the cell (e.g. Cdc42 in *S. Cerevisiae* that is type I essential for cell polarity).
- Type II essential genes are the ones that look essential upon inhibition, but the effects of its inhibition can be rescued or masked by the deletion of (an)other gene(s).
- These genes are therefore not actually essential, but when inhibiting the genes some toxic side effects are provoked that are deleterious for the cells.

The idea is that the essentiality of genes (both type I and type II), may change between different genetic backgrounds.
For changes in essentiality four cases are considered:

1. A gene is **essential** in WT and remains **essential** in the
    mutant.

2. A gene is **non-essential** in WT and remains **non-essential** in
    the mutant.

3. A gene is **essential** in WT and becomes **non-essential** in the
    mutant.

4. A gene is **non-essential** in WT and becomes **essential** in the
    mutant.

An example is given in the figure below, where an interaction map is
shown for WT cells and a possible interaction map for a mutant where
both the essentiality and the interactions are changed.

![Example interaction network for a WT cell and a mutant with essential and non-essential genes.](./media/example_interaction_network.png)

Situation 1 and 3 are expected to be the trickiest since those ones are difficult to validate.
To check the synthetic lethality in cells, a double mutation needs to be made where one mutation makes the genetic background and the second mutation should confirm whether the second mutated gene is actually essential or not.
This is typically made by sporulating the two mutants, but deleting a gene that is already essential in wild type prevents the cell from growing or dividing and can therefore not be sporulated with the mutant to create the double deletion.
Therefore, these double mutants are hard or impossible to make.

Another point to be aware of is that some genes might be essential in specific growth conditions (see also the subsection).
For example, cells that are grown in an environment that is rich of a specific nutrient, the gene(s) that are required for the digestion of this nutrient might be essential in this condition.
The SATAY experiments will therefore show that these genes are intolerant for transposon insertions.
However, when the cells are grown in another growth condition where mainly other nutrients are present, the same genes might now not be essential and therefore also be more tolerant to transposon insertions in that region.
It is suggested to compare the results of experiments with cells from the same genetic background grown in different conditions with each other to rule out conditions specific results.

```{margin} The yeast interaction network
The yeast interaction network is already made based on previous research see here >> [thecellmap.org](thecellmap.org).

```

We want to check the essentiality of all genes in different mutants and compare this with both wild type cells and with each other.
The goal is to make an overview of the changes in the essentiality of the genes and the interaction network between the proteins.
With this we aim to eventually be able to predict the synthetic lethality of multiple mutations based on the interaction maps of the individual mutations.

```{margin} Current progress in Machine Learning workflows! 

Currently the idea is to use machine learning that uses the results from the transposon sequencing experiments, the interaction map of genes and possibly GO-terms. Some implementations of Machine learning workflows can be seen [here, in another Jupyter Book :)](https://leilaicruz.github.io/jupyter-book/Machine-learning-related-projects/prot_domains2lethality-PAPER-REPLICATION.html)

```



## Interpreting Transposon Counts & Reads

Once cells have a transposon inserted somewhere in the DNA, the cells are let to grow so they can potentially generate a few new generations.
A cell with a transposon inserted in an essential part of its DNA grows very slowly or might not grow at all (due to its decreased fitness).
Since the sequencing starts only at the location of a transposon insertion (see experimental methods section), it can be concluded that roughly each read from the sequencing corresponds with a transposon insertion (roughly mainly because transposon inserted in essential genes can generate no reads).
Cells with a transposon inserted in an essential genomic region, will not have divided and therefore will not have a contribution to the sequencing reads.
When the reads are aligned to a reference sequence and the number of reads are mapped against the genome, empty regions indicate possible essential genes.
Negative selection can thus be found by looking for empty regions in the reads mapping.
When a transposon is inserted in a non-essential genomic region, these cells can still divide and give rise to offspring and after sequencing the non-essential regions will be represented by relatively many reads.

During processing the genes can be analysed using the number of transposon insertions per gene (or region) or the number of reads per gene.
Reads per gene, instead of transposons per gene, might be a good measure for positive selection since it is more sensitive (bigger difference in number of reads between essential and non-essential genes), but also tends to be nosier.
Transposons per gene is less noisy, but is also less sensitive since the number of transposons inserted in a gene does not change in subsequent generations of a particular cell.
Therefore it is hard to tell the fitness of cells when a transposon is inserted a non-essential region solely based on the number of transposon insertions.

Ideally only the transposons inserted in non-essential genomic regions will have reads (since only these cells can create a colony before sequencing), creating a clear difference between the essential and non-essential genes.
However, sometimes non-essential genes also have few or no transposon insertion sites.
According to [Michel et.al.](https://elifesciences.org/articles/23570) this can have 3 main reasons.

1. During alignment of the reads, the reads that represent repeated DNA sequences are discarded, since there is no unique way of fitting them in the completed sequence.
 (Although the DNA sequence is repeated, the number of transposon counts can differ between the repeated regions)
 Transposons within such repeated sequences are therefore discarded and the corresponding reads not count.
 If this happens at a non-essential gene, it appears that it has no transposons, but this is thus merely an alignment related error in the analysis process.

2. Long stretches of DNA that are without stop codons, called Open Reading Frames (ORF), typically code for proteins.
 Some dubious ORF might overlap with essential proteins, so although these ORF themselves are not essential, the overlapping part is and therefore they do not show any transposons.

3. Some genes are essential only in specific conditions.
 For example, genes that are involved in galactose metabolism are typically not essential, as inhibition of these genes block the cell’s ability to digest galactose, but it can still survive on other nutrition’s.
 In lab conditions however, the cells are typically grown in galactose rich media, and inhibiting the genes for galactose metabolism cause starvation of the cells.

4. A gene might not be essential, but its deletion might cripple the cell so that the growth rate decreases significantly.
 When the cells are grown, the more healthy cells grow much faster and, after some time, occur more frequently in the population than these crippled cells and therefore these cells might not generate many reads or no reads at all.
 In the processing, it might therefore look as if these genes are essential, but in fact they are not.
 The cells just grow very slowly.

The other way around might also occur, where essential genes are partly tolerant to transposons.
This is shown by Michel et.al. to be caused that some regions (that code for specific subdomains of the proteins) of the essential genes are dispensable.
The transposons in these essential genes are clearly located at a specific region in the gene, the one that codes for a non-essential subdomain.
However, this is not always possible, as in some cases deletion of non-essential subdomains of essential genes create unstable, unexpressed or toxin proteins.
The difference in essentiality between subdomains in a single protein only happens in essential genes, not in non-essential genes.
Michel et.al. devised an algorithm to estimate the likelihood $L$ of a gene having an essential subdomain:

$$L = \frac{\text{d }N_{\text{cds}}}{l_{\text{cds}}}$$

where $d$ is the longest interval (in terms of base pairs) between 5 neighbouring transposons in a Coding DNA Sequence (cds) ($\geq 300$ bp), $N_{cds}$ is the total number transposons mapping in the cds ($\geq 20$) transposons) and $l_{cds}$ is the total length of the CDS.
Additionally, it must hold that $0.1 l_{cds} \leq d \leq 0.9 l_{cds}$.

It is expected that only essential genes carry essential subdomains, and indeed what was found by Michel et.al. that the genes with the highest likelihood were mostly also genes previously annotated as essential by other studies.

Because of the reasons mentioned before, not a simple binary conclusion can be made solely based on the amount of transposon insertions or the number of reads.
Instead, a gene with little reads *might* be essential, but to be sure the results from other experiments need to be implemented as well, for example where the cells were grown in a different growth conditions.
Therefore, SATAY analysis only says something about the relative fitness of cells where a specific gene is inhibited in the current growth conditions.

### Methods and File types

Many essential genes in wild type cells are already known and published.
The results from our SATAY experiments we can verify using the known essential genes and check how well they fit with the number of reads in each gene.
Also, we can use a similar approach as was performed by the Kornmann lab [Michel et.al. 2017] using wild type cells. 



### Experimental method
 
```{admonition} Goal
Automatically converts the raw sequencing data to a list of the number of transposons and reads for each gene, and potentially list some of their properties (e.g. essentiality, GO-terms, interactions etc.).

```

```{admonition} Small quality check
Comparison of our results with those obtained by the Kornmann lab might confirm the quality of our experimental and analysis methods.

```
### Experimental Process 
- The process of SATAY starts with inserting a plasmid in the cells that contains a transposase (TPase) and the transposon (MiniDs) flanked on both sides by adenine (ADE).
The transposon has a specific, known, sequence that codes for the transposase that cuts the transposon from the plasmid (or DNA) to (another part of) the DNA.

[Simplified example for the transposon insertion plasmid.](./media/Plasmid_transposon.png)


- The MiniDs transposon is cut loose from the plasmid and randomly inserted in the DNA of the host cell.
- If the transposon is inserted in a gene, the gene can still be transcribed by the ribosomes, but typically cannot be (properly) translated in a functional protein.
- The genomic DNA (with the transposon) is cut in pieces for sequencing using enzymes, for example DpnII.
- This cuts the DNA in many small pieces (e.g. each 75bp long) and it always cuts the transposon in two parts (i.e. digestion of the DNA).
- Each of the two halves of the cut transposon, together with the part of the gene where the transposon is inserted in, is ligated meaning that it is folded in a circle.
- A part of the circle is then the half transposon and the rest of the circle is a part of the gene where the transposon is inserted in.
- Using PCR and primers, this can then be unfolded by cutting the circle at the halved transposon.
- The part of the gene is then between the transposon quarters.
- Since the sequence of the transposon is known, the part of the gene can be extracted.
- This is repeated for the other half of the transposon that includes the other part of the gene.
- When both parts of the gene are known, the sequence from the original gene can be determined.

![Schematic overview of transposon insertion experiments.](./media/satay_experiment_overview.png)

#### Sequence alignment

- To get the order of nucleotides in a genome, shotgun sequencing is used where the genome is cut is small pieces called reads (typically tens to a few hundred basepairs long).
- The reads have overlapping regions that can be used to identify their location with respect to a reference genome and other reads (i.e. mapping of the reads).
- Mapping of the reads result in contigs, which are multiple mapped reads that form continuous assembled parts of the genome (contigs can be the entire target genome itself).
- All contigs should be assembled to form (a large part of) the target genome.

- The sequence assembly problem can be described as: *Given a set of sequences, find the minimal length string containing all members of the set as substrings*.

The reads from the sequencing can be single-end or paired-end, which indicates how the sequencing is performed.
In paired-end sequencing, the reads are sequenced from both directions, making the assembly easier and more reliable, but results in twice as many reads as single-end reads.
The reason of the more reliable results has to do with ambiguous reads that might occur in the single-end sequencing.
Here, a read can be assigned to two different locations on the reference genome (and have the same alignment score).
In these cases, it cannot be determined where the read should actually be aligned (hence its position is ambiguous).
In paired-end sequencing, each DNA fragment has primers on both ends, meaning that the sequencing can start in both the 5’-3’ direction and in the 3’-5’ direction.
Each DNA fragment therefore has two reads both which have a specified length that is shorter than the entire DNA fragment.
This results that a DNA fragment is read on both ends, but the central part will still be unknown (as it is not covered by these two particular reads, but it will be covered by other reads).
Since you know that the two reads belong close together, the alignment of one read can be checked by the alignment of the second read (or paired mate) somewhere in close vicinity on the reference sequence.
This is usually enough for the reads to become unambiguous.

The resulting data from the sequencing is stored in a FASTQ file where all individual reads are stored including a quality score of the sequencing.
The reads are in random order and therefore the first step in the processing is aligning of the reads in the FASTQ files with a reference genome.

Note that the quality of the reads typically decreases near the 3'-end of the reads due to the chemistry processes required for sequencing (this depends on the kind of method used).
For Illumina sequencing, the main reasons are signal decay and dephasing, both causing a relative increase in the background noise.
Dephasing occurs when a DNA fragment is not de-blocked properly.
A DNA fragment is copied many times and all copies incorporate a fluorescent nucleotide that can be imaged to identify the nucleotide.
If there are 1000 copies of the DNA fragment, there are 1000 fluorescent nucleotides that, ideally, are all the same to create a high quality signal.
After imaging, the DNA fragment is de-blocked to allow a new fluorescent nucleotide to bind.
This deblocking might not work for all copies of the DNA fragment.
For example, 100 copies might not be deblocked properly, so for the next nucleotide only 900 copies agree for the next incorporated nucleotide.
For the third round, the 100 copies that were not deblocked properly in the second round, might now be deblocked as well, but now they are lagging behind one nucleotide, meaning that in the coming rounds they have consistently the wrong nucleotide.
As the reads increases in length, more rounds are needed and therefore the chances of dephasing increases causing a decrease in the quality of the reads.
This gives noise in the signal of the new nucleotide and therefore the quality of the signal decreases.
For example, take the next 6bp sequence that is copied 5 times:

1. `GATGTC`
2. `GATGTC`
3. `G ATGT`
4. `GAT GT`
5. `G AT G`

The first two reads are deblocked properly and they give all the right nucleotides.
But the third and fourth have one round that is not deblocked properly (indicated by the empty region between the nucleotides), hence the nucleotide is always lagging one bp after the failed deblocking.
The fifth copy has two failed deblocking situations, hence is lagging two bp.
The first nucleotide is a G for all 5 copies, therefore the quality of this nucleotide is perfect.
But, by the end of the read, only two out of five copies have the correct bp (i.e. a C), therefore the quality is poor for this nucleotide.
(It can either be a C or a T with equal likelyhood or potentially a G, so determining which nucleotide is correct is ambiguous without knowing which reads are lagging, which you don't know).
(See for example [this question on seqanswers](http://seqanswers.com/forums/showthread.php?t=61198) or the paper by [Pfeifer, 2016])

## FASTQ files and FASTA files

The output of sequencing is typically a FASTQ file.
This file contains the sequences of all reads (typically 75 to 100 bp long).
Each read is represented by four lines:

1. Starts with ‘@’ followed by a sequences identifier and can include
    some descriptions.

2. Raw sequence letters representing the nucleotides.

3. Starts with ‘+’ for separating the second and third line.

4. Quality score of the sequence represented by ASCII symbols running
    from ‘\!’ (lowest score) to ‘\~’ (highest score)
    [<http://www.asciitable.com/>]. This is called ASCII-base 33 since
    ‘\!’ has decimal ASCII number 33 and is defined as Q-score 0. This
    typically runs towards ASCII symbol ‘J’ (number 74, Q-score 41). The
    error probability can be calculated based on the Q-score using
    $P_{error} = 10^{-\frac{Q}{10}}$. This means that ‘\!’ has
    an error of ($P_{error} = 100\%$) and ‘J’ an error of
    $P_{error} = 0.008\%$). Typically a Q-score higher than Q=20
    (ASCII symbol ‘5’, ($P_{error} = 1\%$)) is acceptable
    [https://drive5.com/usearch/manual/quality_score.html]. This
    line has the same length as the second line.

A FASTQ file can be used as an input for sequence alignment using a reference genome.
The result of sequence alignment is a SAM file.

Besides FASTQ files, FASTA files are also used for alignment.
These are similar to FASTQ files, but do not include the quality string and therefore FASTA files can be created from FASTQ files by removing line 3 and 4 from each read.
Depending on the sequencing method, FASTA files may be directly given.

#### SAM and BAM files

The FASTQ (or FASTA) files contain all the reads in a random order.
To determine where each belong relative to a reference sequence, the reads
need to be aligned.
After alignment of the reads, the results are typically represented in a Sequencing Alignment Mapping (SAM) file.
For processing purposes this is typically translated to a Binary Alignment Mapping (BAM) file.
This is a compressed, binary version of the SAM file.

The SAM files contain all the short reads from the FASTQ file together with information where and how well the read mapped to a reference
sequence.
The reads are represented by a header and a body.
The header consists of four fields (note that the headers can be different depending which program made the alignment, see for example [cpwardell.com](http://www.cpwardell.com/2015/02/24/understanding-bam-files-part-1/)):

1. @HD lines: version number of SAM specification and how the file is
    sorted (e.g. sorting by genomic coordinates).

2. @SQ: Which reference sequence has been used to align the sequence.
    There will be one line for every chromosome (contig). It also tells
    the sequence name (SN) and length (LN).

3. @RG: Read group line with the tags specifying the origins of the
    input data.

4. @PG: Which programs and commands were used to create the SAM file.

### Readme for the SAM file 

**The body of the SAM file contains the aligned data.**
- Every row contains one read.
- If there is paired data (i.e. a forward reading and a backward reading), then the pair is divided in two rows.
- Every row consists of at least 11 columns:

1. `QNAME` Name of the query sequence.

2. `FLAG` Bitwise flag.
 This consists of twelve binary properties.
 A read typically has multiple flags.
 These flags are all then translated to a decimal number, given by the third column, and these decimal numbers are added up.
 Typical values are 99, 147, 83 or 163.
 To get a proper translation, use [the SAM flag decoder](http://broadinstitute.github.io/picard/explain-flags.html).
 The following flags can be used:

    1. 000000000001 : 1 : read paired

    2. 000000000010 : 2 : read mapped in proper pair

    3. 000000000100 : 4 : read unmapped

    4. 000000001000 : 8 : mate unmapped

    5. 000000010000 : 16 : read reverse strand

    6. 000000100000 : 32 : mate reverse strand

    7. 000001000000 : 64 : first in pair

    8. 000010000000 : 128 : second in pair

    9. 000100000000 : 256 : not primary alignment

    10. 001000000000 : 512 : read fails platform/vendor quality checks

    11. 010000000000 : 1024 : read is PCR or optical duplicate

    12. 100000000000 : 2048 : supplementary alignment

3. `RNAME` Name of the reference contig (chromosome) where the sequence is aligned to (i.e. which chromosome the read is aligned to).

4. `POS` Position of the reference contig that the alignment starts at (given in terms of base pairs).

5. `MAPQ` Mapping quality.
Number indicating the chances that the mapping is wrong, based on phred scaling.
This is logarithmic scaled where 60 is typically the maximum score meaning that the chance of a wrong mapping is the smallest (so a high number is better).
If a value of 255 is shown, that means that the quality is not determined.

6. `CIGAR` Tells how to match the query sequence to the reference sequence using a ‘Compact Idiosyncratic Gapped Alignment Report’ (CIGAR) string.
This contains a sequence of integers and letters.
Possible letters are M (Match), N (Alignment gap), D (Deletion) or I (Insertion).
Thus 76M means that 76 basepairs match the reference sequence (see [JEFworks](https://jef.works/blog/2017/03/28/CIGAR-strings-for-dummies/) for more information).

7. `RNEXT` The name of the reference contig (chromosome) that the other read in the pair (i.e. the next or previous line?) aligns to.
If the two reads in the pair aligns to the same contig, an ‘=’ sign is used.

8. `PNEXT` Position on the contig where the other read in the pair aligns to.
Depending on the size of the DNA fragments in the sequencing library, this is typically a few hundred base pairs away from the current read (i.e. given by the POS column).

9. `TLEN` Total length of the template fragment.
This is the distance from the leftmost base of the first read to the rightmost base pair of the second read in the pair.
This value is assigned a ‘+’ or ‘-‘ to indicate the reading orientation.

10. `SEQ` The DNA sequence of the query sequence.
(Identical to the sequence in the FASTQ file that was aligned to the reference genome).

11. `QUAL` Base quality score of the SEQ.
(Identical to the scores in the FASTQ file).
There are 42 scores, each of which are related to a specific error.
See for example [phred score conversion table](https://drive5.com/usearch/manual/quality_score.html) for a conversion table.


#### Determine essentiality based on transposon counts

Using the number of transposons and reads, it can be determined which genes are potentially essential and which are not.
To check this method, the transposon count for wild type cells are determined.
Currently, genes that are taken as essential are the annotated essentials based on previous research. 

#### Distribution number of insertions and reads compared with essential and non-essential genes

See [this script to visualize the implementation](../statistics_pergene.md)


Ideally, the number of transposon insertions of all essential genes are small and the number of insertions in non-essential genes are is large so that there is a clear distinction can be made.
However, this is not always so clear.
For example, the distribution of transposons in WT cells in the data from Michel et. al. looks like this:

![Reads and transposon density after processing as discussed in Michel et.al. 2017.](./media/Tn_Reads_Pergene_BenoitAnalysis.png)

In this figure, both the reads and the transposon counts are normalized with respect to the length of each gene (hence the graph represents the read density and transposon density).
High transposon counts only occurs for non-essential genes, and therefore when a high transposon count is seen, it can be assigned nonessential with reasonable certainty.
However, when the transposon count is low the there is a significant overlap between the two distributions and therefore there is no certainty whether this gene is essential or not (see also the section about 'Interpreting Transposon Counts & Reads').

The data is also sensitive to postprocessing.
It is expected that the trimming of the sequences is an important step.
The graph below shows the same data as in the previous graph, but with different processing as is done by Michel et. al..
This has a significant influence on the results and as a consequence, no distinction can be made between essential and nonessential genes based on the transposon counts.
Significant attention needs to be given to the postprocessing of the data.

![Reads and transposon density after processing.](./media/Tn_Reads_Pergene_MyAnalysis.png)

#### Profile plot for number of reads

(*See `TransposonRead_Profile_Plot.py`*)

To create a visual overview where the insertions are and how many reads there are for each insertion, a profile plot is created for each chromosome.

![Read profile plot for chromosome XV (note the y-axis is in logarithmic scale).](./media/Read_ProfilePlot_chrxv.png) 

The bars indicate the absolute number of reads for all insertions located in the bars (bar width is 545bp).
The colored background indicate the location of genes, where green are the annotated essential genes and red the non-essential genes.
In general, the essential genes have no or little reads whereas the non-essential genes have many reads.
Note that at location 564476 the ADE2 gene is located that has significant more reads than any other location in the genome, which has to do the way the plasmid is designed (see Michel et.al. 2017).
The examples used here are from a dataset discussed in the paper by Michel et.al. 2017 which used centromeric plasmids where the transposons are excised from.
The transposons tend to reinsert in the chromosome near the chromosomal pericentromeric region causing those regions to have about 20% more insertions compared to other chromosomal regions.

This figure gives a rough overview that allows for identifying how well the data fits the expectation.
Also an alternative version of this plot is made (`TransposonRead_Profile_Compare.py`) that makes a similar plot for two datasets at the same time, allowing the comparison of the two datasets with each other and with the expectation.

![Comparison of the same datasets, but with different processing steps. Shown here is the transposon count for the two files including the absolute difference between the two datasets show in blue. Note also here that some regions has a higher likelihood of bearing transposons compared to the surrounding regions.](./media/Cerevisiae_Michel2017_WT2_Compare_chromIX.png)

### Profile plot number of reads per individual genes

See [this script to visualize the implementation](../gene_reads.md), which is {doc}`../gene_reads.md`

Instead of plotting the number of reads for an entire chromosome, it is also useful to plot the read profile for individual genes to see how the insertion sites and reads are distributed within a gene.
For this a bar plot is made where the number of reads per transposon are determined.
This also shows the distribution of the distances between subsequent transposon insertions for both the gene and the chromosome the gene is located.
It is expected that for essential genes, the median distance between subsequent insertions is larger compared to the distance in the entire chromosome (since important genes have few insertions and large parts will be free of insertions).
For non-essential genes, the distribution for the distance between insertions is expected to follow the distribution of the chromomsome more closely.

![Read per transposon insertion for BEM1. The violin plot gives the distribution for the distance between subsequent insertions for both BEM1 and chromosome II (where BEM1 is located). The small black bars indicate individual transposon insertions.](./media/Read_ProfilePlot_bem1.png)

 ![Read per transposon insertion for the HO-locus. Note that the distribution for the distance between insertions follows the distribution for the chromosome more closely compared to BEM1.](./media/Read_ProfilePlot_ho.png) 

The minimum width of the bars are chosen such that each bar should contain 8 transposon insertions.
The maximum width is set equal to the length such that the probabilty of finding at least one insertion is more than 99% in the whole chromosome.
This is chosen because if now a bar is empty in a gene than this is not a coincidence, but this is an interesting region.

These plots can be used for checking if a gene has transposon free regions which might indicate that this gene is essential.