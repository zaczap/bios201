Workshop2: Aligning and quantifying reads from RNA-seq
======================================================

Introduction
------------

Today, we are going to process RNA sequencing data into gene expression
values. We will be working with data from [a published study] 
(http://journals.plos.org/plosone/article/file?type=supplementary&id=info:doi/10.1371/journal.pone.0097550.s002).
comparing RNA seq from healthy lung samples to lung samples from patients with
idiopathic pulmonary fibrosis (IPF).

We have provided you with FASTQ files from 8 normal and 8 IPF
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

Let's first get familiar with a couple file formats used to specify gene
annotations.

### GenePred format
This format is also referred to as refFlat by Picard tools. This is a file
format used by UCSC to encode predicted gene positions. Each transcript is
described by a single line. Read about the format
[here](http://genome.ucsc.edu/FAQ/FAQformat.html#format9).

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
grep "SFTPA2" annotation/UCSC_table_browser_chr10.txt | column -t | less -S
```
(`column -t` helps line up columns in a human-readable way)

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

We will use STAR to align our RNA-seq reads. STAR performs best if you can
provide it with a gene annotation (in GTF format) as it uses the
information about known splice junctions. However, STAR also infers novel
splice junctions based on the data, which it outputs. The tool suggests
that after doing an initial alignment, you provide STAR with the junctions
it has infered and rerun the mapping. This is what it calls 2-pass
mapping.

Before aligning reads with STAR, we need to build the genome index.
If you have known genes in GTF format, you can provide it at this step and
we will.

```
## Do not run this, we did it for you already to save time.
STAR --runThreadN 20 \
     --runMode genomeGenerate \
     --genomeDir STAR_hg19_chr10 \
     --genomeFastaFiles chr10.fa \
     --sjdbGTFfile annotation/gencode.v19.annotation_chr10.gtf \
     --sjdbOverhang 100  # readLength - 1
```

Now we can map the samples.
It takes a few minutes to run each alignment, so we will only have you map
two of the 16 samples and will provide you with the results for the rest.

The generic mapping command we will be using
```
STAR --runThreadN <NumberOfThreads> \
     --genomeDir </path/to/genomeDir> \
     --readFilesIn </path/to/read1> </path/to/read2> \
     --readFilesCommand <UncompressionCommand> \
     --outFileNamePrefix </path/to/output> \
     --outSAMtype <OutputFormat> \
     --sjdbFileChrStartEnd <JunctionFile> \
     --quantMode <QuantificationType>
```
The first three arguments are required. The rest are optional but we will
use them to provide non-default values. The last two we will only use for
the second-pass alignment.

**NOTE:** that there are *lots* of parameters that you can adjust. The default
parameters work well for us, but there is a note on the STARs homepage:
"This release was tested with the default parameters for human and mouse
genomes. Please contact the author for a list of recommended parameters
for much larger or much smaller genomes." Keep that in mind if you want to
use STAR for other organisms.

<!---STAR runs more quickly than some other spliced aligners because it loads
its genome representation into memory. If you are working with a large
genome and need to map multiple samples, you can instruct STAR to load it
once and share that genome between multiple processes.---> 

Run the first-pass alignment for Norm1 and IPF1.
```
STAR --runThreadN 4 \
     --genomeDir STAR_hg19_chr10 \
     --readFilesIn fastq/Norm1_R1.fastq.gz fastq/Norm1_R2.fastq.gz \
     --readFilesCommand gunzip -c \
     --outFileNamePrefix bam_pass1/Norm1_ \
     --outSAMtype BAM Unsorted

STAR --runThreadN 4 \
     --genomeDir STAR_hg19_chr10 \
     --readFilesIn fastq/IPF1_R1.fastq.gz fastq/IPF1_R2.fastq.gz \
     --readFilesCommand gunzip -c \
     --outFileNamePrefix bam_pass1/IPF1_ \
     --outSAMtype BAM Unsorted

```

Then run the second-pass alignment inputting the junction files for **all
16** samples. We'll also get STAR to output the number of reads mapping to
each gene.
```
# This collects the names of all the junction files from the first pass.
# We can then pass this information to the aligner below.
junctions=`ls bam_pass1/*_SJ.out.tab`

STAR --runThreadN 4 \
     --genomeDir STAR_hg19_chr10 \
     --readFilesIn fastq/Norm1_R1.fastq.gz fastq/Norm1_R2.fastq.gz \
     --readFilesCommand gunzip -c \
     --outFileNamePrefix bam_pass2/Norm1_ \
     --outSAMtype BAM Unsorted \
     --sjdbFileChrStartEnd $junctions \
     --quantMode GeneCounts

STAR --runThreadN 4 \
     --genomeDir STAR_hg19_chr10 \
     --readFilesIn fastq/IPF1_R1.fastq.gz fastq/IPF1_R2.fastq.gz \
     --readFilesCommand gunzip -c \
     --outFileNamePrefix bam_pass2/IPF1_ \
     --outSAMtype BAM Unsorted \
     --sjdbFileChrStartEnd $junctions \
     --quantMode GeneCounts

```

Understanding SAM/BAM format
----------------------------

Summaring mapping results
-------------------------

Visualizing alignments with IGV
-------------------------------