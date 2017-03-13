Workshop2: Aligning and quantifying reads from RNA-seq
======================================================

Introduction
------------

Today, we are going to process RNA sequencing data into gene expression
values. We will be working with data from [a published study] 
(http://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0097550.s002).
comparing RNA seq from healthy lung samples to lung samples from patients with
idiopathic pulmonary fibrosis (IPF).

We have provided you with FASTQ files from 8 normal and 8 control
samples. For this activity, we will be restricting our analysis to
chromosome 10 to speed up some of the running times.

We will be mapping reads to build hg19 of the human genome. This build
is very similar to GRCh37, which we used last time. The primary difference
is in the naming of the chromosomes: in hg19, the chromosome names include
"chr", whereas in GRCh37 they do not. We chose hg19 because the chromosome
names match those in the gene annotation we are working with and that is
important for the relevant programs to run correctly. (As a side note, a
common source or errors and bugs when working with genomic data is a
mismatch in the chromosome naming schemes between files you are using.) 
There are a few other differences between the genome builds that you can read about
[here](https://wiki.dnanexus.com/Scientific-Notes/human-genome) if you are
interested. 

Goals
-----
At the end of this workshop you should be familiar with the following:
* Gene annotation file formats
* Spliced alignment and gene count generation for RNA-seq data
* Basics of SAM/BAM format
* Visualization of read alignments
* Picard's summary of RNA-seq alignments
* Preliminary inspection of count data in R

Setup
-----
Follow the setup steps outlined in the first workshop (see
[here](https://github.com/zaczap/bios201/blob/master/setup.md) for a refresher).  
Then copy the workshop materials and `cd` into the copied directory:
```
cp -r /afs/ir/class/bios201/workshop2/ .
cd workshop2
```

Understanding gene annotation formats
-------------------------------------

### GenePred format
This format is also referred to as refFlat by Picard tools. This is a file
format used by UCSC (read about the format
[here](http://genome.ucsc.edu/FAQ/FAQformat.html#format9)) to encode
predicted gene positions. Each transcript is described by a single line.

We have retrieved a GenePred file of genes on chromosome 10 from UCSC's
[table browser](https://genome.ucsc.edu/cgi-bin/hgTables) for you.

This is the query we ran:
![screenshot of UCSC table browser](http://web.stanford.edu/class/bios201/table_browser_screenshot.png)

Let's take a look at what this file looks like:
```
less -S annotation/UCSC_table_browser_chr10.txt
(use up/down and side arrows to navigate the file and type q to quit)
```

Since we're working with lung samples, let's check out Surfactant Protein
A2 (SFTPA2), a gene highly expressed in lung.
```
grep "SFTPA2" annotation/UCSC_table_browser_chr10.txt | less -S
```

:question: **How many different transcripts are there for this gene? How
	   many exons per transcript?**  
:question: **Can you figure out how these transcripts relate to each other
	   (which exons are shared, which are different)?**

### GTF/GFF format
Gene Transfer Format (GTF) is an extension of the General Feature Format
(GFF). Take a minute to read about those formats
[here](http://genome.ucsc.edu/FAQ/FAQformat.html#format3). GTF files are
used by many RNA-seq aligners and tools for counting the number of reads
mapping to each gene. In contrast to the GenePred format, each "feature"
occupies a line. There will often be a feature for a gene, a separate
feature for each transcript, and features for each exon, among
others. Therefore, the detailed exon information for each transcript is
spread over multiple lines.

We have provided you with the GTF for the version 19 comprehensive gene annotation available
from [Gencode](https://www.gencodegenes.org/releases/19.html), but
subsetted to chromosome 10.

Let's take a look at it:
```
less -S annotation/gencode.v19.annotation_chr10.gtf
```

Again, let's look at SFTPA2.
```
grep "SFTPA2" annotation/gencode.v19.annotation_chr10.gtf | less -S
```
You should be able to identify the same transcripts as in the GenePred format.

:question: **Do you notice anything different in the order in which the
	   exons are listed compared the GenePred format?**

Spliced Alignment and read counting
-----------------------------------

Understanding SAM/BAM format
----------------------------

Visualizing alignments with IGV
-------------------------------