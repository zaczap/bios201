Workshop2: Aligning and quantifying reads from RNA-seq
======================================================

Introduction
------------

Today, we are going to process RNA sequencing data into gene expression
values. We will be working with data from [a published study comparing RNA
seq from healthy lung samples to lung samples from patients with
idiopathic pulmonary fibrosis
(IPF)](http://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0097550.s002).

We have provided you with trimmed FASTQ files from 8 normal and 8 control
samples.

We will be mapping reads to build hg19 of the human genome. This build
is very similar to GRCh37, which we used last time. The primary difference
is in the naming of the chromosomes: in hg19, the chromosome names include
"chr", whereas in GRCh37 they do not. We chose hg19 because the chromosome
names match the gene annotation we are working with and that is
important for the relevant programs to run correctly. (As a side note, a
common source or errors and bugs when working with genomic data is a
mismatch in the chromosome naming schemes between files you are using.) 
There are a few other minor differences that you can read about
[here](https://wiki.dnanexus.com/Scientific-Notes/human-genome) if you are
interested. 

For this activity, we will be restricting our analysis to
chromosome 10 to speed up some of the running times.
 
