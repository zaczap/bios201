# Workshop 1: Variant calling from DNA sequencing data

The exercises in this tutorial require you to create your own commands to run command-line bioinformatics tools.
To get the most from the tutorial, we recommend that you spend some time thinking about how to form these commands, 
and experiment to see what does and doesn’t work. However, if you get stuck or you just need to see the 
fully-formed commands right away, you can instead view [this](https://github.com/zaczap/bios201/tree/master/Workshop1/README_easy.md) version of the tutorial.
The answers to the in-tutorial questions are given [here](https://github.com/zaczap/bios201/tree/master/Workshop1/Answers.md).

_For non-Stanford students_:
The instructions in this tutorial assume you have access to Stanford's `rice` server. Luckily, if you're not at 
Stanford, you can still work through this tutorial on your own Linux machine, using [these files](http://web.stanford.edu/class/bios201/). 
The only tricky part is that you'll have to install the software tools yourself.

## Introduction

In this exercise, we are going to _identify variants_ from DNA sequencing data for a few individuals. A variant (either a base substitution, a short insertion/deletion, or a larger structural variant like a duplication, deletion, translocation, etc.) is defined relative to a reference genome - today we will be working with the [GRCh37 assembly](https://www.ncbi.nlm.nih.gov/grc/human) of the human genome. We will essentially focus our analysis on single nucleotide polymorphisms (SNPs).

We have provided you with raw sequencing data in the FASTQ format for three individuals that were sequenced as part of the 1000 Genomes project ([NA12878](http://www.internationalgenome.org/data-portal/sample/NA12878), [NA12891](http://www.internationalgenome.org/data-portal/sample/NA12891), and [NA12892](http://www.internationalgenome.org/data-portal/sample/NA12892)). 

These samples underwent _paired-end_ sequencing on genomic DNA, generating a pair of FASTQ files (R1 and R2) for each sample. We have subsetted the FASTQ files to reads that originate from the _BRCA1_ locus on chromosome 17 (17:41196312-41277500):

Sample | Read 1 | Read 2
-------|--------|--------
NA12878 | NA12878_R1.fastq | NA12878_R2.fastq
NA12891 | NA12891_R1.qc.trimmed.fastq | NA12891_R2.qc.trimmed.fastq
NA12892 | NA12892_R1.qc.trimmed.fastq | NA12892_R2.qc.trimmed.fastq

We have already performed quality control (QC) on NA12891 and NA12892, but you will be responsible for cleaning up NA12878. Broadly we will be following the [GATK Best Practices pipeline](https://software.broadinstitute.org/gatk/best-practices/bp_3step.php?case=GermShortWGS) for calling variants from DNA sequencing data, though we will be skipping one step due to the small amount of data we are working with (we are skipping variant quality score recalibration, or VQSR).

The general workflow will look like this:

![GATK workflow](https://software.broadinstitute.org/gatk/img/BP_workflow_3.6.png) 

We are going to take raw sequencing reads, align the reads to the reference genome, and then call variants from the aligned reads. If you want to learn more about the main file formats we will be working with, the following pages provide information on [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format), [SAM/BAM](http://genome.sph.umich.edu/wiki/SAM), and [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) formats.

## Goals

After going through this exercise, you will have done the following:

- Run `fastqc` on FASTQ files to assess sequencing quality
- Run `cutadapt` to remove adapters, trim reads, and filter short reads
- Align paired-end reads to your indexed genome using `bwa mem`
- Inspect alignment results with `samtools flagstat`
- Mark PCR duplicates in the aligned reads (Picard tools)
- Perform base quality recalibration on the aligned reads (GATK)
- Call variants in individual samples (GATK)
- Jointly call variants across samples (GATK)
- Compare VCF files (in this case, to a gold-standard reference)

While we have provided all the software for you on our teaching server, you may want to download all of this software on your computer in the future: [fastqc](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [cutadapt](http://cutadapt.readthedocs.io/en/stable/index.html), [bwa](http://bio-bwa.sourceforge.net/), [samtools](http://www.htslib.org/), [Picard tools](https://broadinstitute.github.io/picard/), [GATK](https://software.broadinstitute.org/gatk/), [vcftools](http://www.htslib.org/).

## Setup

First you will have to connect to the server we have set up for you:

	ssh <sunetID>@rice.stanford.edu
	# e.g. ssh zappala@rice.stanford.edu

Once connected, run:

	echo $SHELL

If the response is `tcsh`, run the following command:

	source /afs/ir/class/bios201/setup/setup_tcsh.sh 

If the response is `bash`, run the following command:

	source /afs/ir/class/bios201/setup/setup_bash.sh 
	
(If you're curious what's going on here...
Most command-line tools are run using an _executable_ file, a file that contains the instructions needed for the computer to run a task.
The setup commands above configure your session so that you can run these tools using short commands like `$PICARD` instead of typing
out the full path to the executable file every time, e.g. `/afs/ir/users/e/t/etsang/software/bin/picard.jar`. When you're running your own
analysis, these setup commands aren't necessary; instead, you can just include the full path name.)

Finally, copy the workshop materials and `cd` into it (**NOTE**: copying the materials may take a minute or two):

	cp -r /afs/ir/class/bios201/workshop1/ .
	cd workshop1

## Where we are starting from

Read the first few paragraphs on the Wikipedia page about [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) file format and quality scores.

Now let's take a look at a few FASTQ files.

	head -8 NA12878_R1.fastq
	head -8 NA12878_R2.fastq

:question: **How many reads are displayed using each of these commands?**

:question: **Of the reads you viewed, which has the highest quality? Which has the lowest?**

<!-- The first one in R2 is BAD -->

You can also view the entire FASTQ file if you type `less -S <fastq_file>` and then type PgUp/PgDn to navigate the list. Try it right now.

## Sequencing quality control with FASTQC

When working with sequencing data, you want to make sure your reads do not have systematic biases or adapter contamination, 
and that they are of generally good quality. The tool FASTQC can help assess these artifacts in your sequencing data, and is 
fairly easy to run. First let's look at the required input:

	fastqc --help

Often you will have the ability to run `--help` or `-h` in order to learn what a command does. The relevant arguments
for the `fastqc` command are summarized here:

	fastqc --outdir <output_dir> --format <format> <r1_fastq> <r2_fastq>

	# <output_dir> is the name of a folder where you want the results written
	# <format> is the format of the sequencing reads (optional)
	# <r1_fastq> is the first set of reads
	# <r2_fastq> is the second set of reads 

**NOTE:** We are only going to perform QC on NA12878; the other samples are already QC'ed and trimmed.

Using the format described above, run fastQC on the fastq files from NA12878. `<output_dir>` can be 
whatever you want; it’s the name of a new folder where your results will be placed. `format` should 
be `fastq`. (Hint: if you forget the names of the input files, you can use `ls` to view all files in 
the current working directory.)

This should run pretty quickly. When it's done, let's go look at the results:

	ls fastqc

You can see a `.html` file and a `.zip` file for each FASTQ file. If we download the `.html` files, we can look at them
on our local computer to inspect the results from FASTQC. Normally you would need to download these files to your
machine and open them, but we have made them available at the following URLs:

- [NA12878_R1_fastqc.html](http://web.stanford.edu/class/bios201/NA12878_R1_fastqc.html)
- [NA12878_R2_fastqc.html](http://web.stanford.edu/class/bios201/NA12878_R2_fastqc.html)

The FASTQC files can help you diagnose various issues you may have had - we will just consider adapter 
contamination and trimming low-quality bases today. To understand what all of the different FASTQC results 
mean, [visit this page](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/).

:question: **What happens to the quality of sequenced bases in later cycles?**   
:question: **Do you see some adapter contamination in the data?**   

<!-- The qualities drop off in the latter cycles -->
<!-- Towards the end of the reads yes -->

## Removing adapters, trimming reads, and filtering

From the FASTQC results, you can see a little bit of adapter contamination as well a tell-tale drop off in 
sequencing quality towards the end of the reads. We can remove adapter contamination and trim these 
low-quality bases using the tool `cutadapt`:

	cutadapt --help

The relevant fields for us right now are here:

	cutadapt -q <minimum_quality> --minimum-length <minimum_length> -a <adapter on first read> -A <adapter on second read> -o <r1_trimmed_fastq> -p <r2_trimmed_fastq> <r1_fastq> <r2_fastq>

	# <minimum_quality> is the minimum base quality to allow (trimmed otherwise)
	# <minimum_length> is the minimum length both pairs of reads need to be to be included
	# <adapter on first read> for us: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
	# <adapter on second read> for us: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
	# <r1_trimmed_fastq> the name of the output file for the first read file 
	# <r2_trimmed_fastq> the name of the output file for the second read file 
	# <r1_fastq> is the first set of reads
	# <r2_fastq> is the second set of reads 

	# note: the adapter sequence is the reverse complement of the actual adapter sequence

TODO: Using the above command as a template, with the adapters mentioned above, run `cutadapt` on the NA12878 samples, 
using the following settings:

* Minimum base quality: 10
* Minimum read length: 30

:question: **How many read pairs were there before/after trimming?**

<!-- before: 10,844, after: 10,336 -->

Now we've removed adapters and trimmed low quality bases; we've also filtered out very short reads which would probably not align to the genome anyway.

## Aligning paired-end DNA sequencing data to the human genome

We would like to take our reads (in FASTQ format) and figure out where they originated from in the human genome - this __alignment__ 
is accomplished using a tool called `bwa` (which stands for the Burrows-Wheeler aligner). To learn a little bit about what `bwa` is capable of, run:

	bwa

We first need to __index__ the genome before it can be easily aligned to. Much like
the index of a travel guidebook or a biology textbook, creating an index makes it fast
for us to find where a topic (read) is located in the genome, instead of scanning through
the entire genome to match each read. To index the genome, we run `bwa index`:

	# You don't need to do this today! 
	# We did it for you already.
	# It would take ~5 minutes (longer on a full genome).
	bwa index grch37.fa

**NOTE:** The genome FASTA file we have given you only includes chromosome 17.

Now we can align reads using the following command:

	bwa mem <genome> <r1_fastq> <r2_fastq> | samtools sort --output-fmt BAM > <bam_file>

	# <genome> is the genome we are aligning to - it must be indexed before it can be used!
	# <r1_fastq> the first set of trimmed reads
	# <r2_fastq> the second set of trimmed reads
	# <bam_file> the alignment file where we will save results

Let's do that for NA12878 (normally this can take some time, but we have only given you a small set of reads):

	bwa mem grch37.fa NA12878_R1.qc.trimmed.fastq NA12878_R2.qc.trimmed.fastq \
		| samtools sort --output-fmt BAM > NA12878.bam

This command first aligns the trimmed FASTQ reads to the reference genome. It then
uses the pipe `|` operator to pass the results into `samtools`, which uses the `sort`
command to sort the reads by the order in which they appear in the genome. The resulting
file is a sorted list of reads in a compressed alignment file, called a BAM.

Now form the commands yourself to repeat this process for the next two samples.

We can see various statistics about each of the alignments by running `samtools flagstat`.

:question: **What percentage of reads for NA12878 were successfully mapped?** (Hint: Don't forget to tell `samtools`
which BAM file it should look at.) 
<!-- 93.86% -->

If you're interested in seeing the sorted alignment file, you can do this using `samtools view NA12878.bam | less -S`.
This uses `samtools view` to read the compressed BAM file, and then pipes the result to `less` so you can scroll 
through the file manually.

**NOTE**: Alignment files (SAM/BAM format) are important to understand well - we won't look at them in detail today, 
but in future workshops you will use them a bit more.

## Marking PCR duplicates and performing base quality recalibration

When we call variants, we want to make sure we have good evidence that the alternate allele exists and isn't noise. 
If PCR duplicates are present, we may artificially think that we have a lot of evidence for an alternate allele that 
is really just the original same DNA fragment that has been over-amplified. We also want to take into account the 
quality of sequencing for each base, since a low quality base should intuitively be poor evidence for a variant 
but a high quality base is good evidence.

We use a tool here called Picard which is a piece of software with many useful functions that you can learn about 
[here](https://broadinstitute.github.io/picard/) or by running the help command:

	java -Xmx2g -jar $PICARD

You can also learn more about GATK, the tool we will use to calibrate variants and call variants, by running the help command:

	java -jar $GATK --help
    
This section isn't super exciting, but it's required to call variants, so we've just included commands that you can
copy and paste to achieve each of these tasks.

You can mark duplicates in your BAM files like so:

	java -Xmx2g -jar $PICARD MarkDuplicates \
		INPUT=NA12878.bam \
		OUTPUT=NA12878.markduplicates.bam \
		METRICS_FILE=NA12878.markduplicates.metrics.txt \
		OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
		CREATE_INDEX=true \
		TMP_DIR=/tmp

	java -Xmx2g -jar $PICARD MarkDuplicates \
		INPUT=NA12891.bam \
		OUTPUT=NA12891.markduplicates.bam \
		METRICS_FILE=NA12891.markduplicates.metrics.txt \
		OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
		CREATE_INDEX=true \
		TMP_DIR=/tmp

	java -Xmx2g -jar $PICARD MarkDuplicates \
		INPUT=NA12892.bam \
		OUTPUT=NA12892.markduplicates.bam \
        METRICS_FILE=NA12892.markduplicates.metrics.txt \
		OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
		CREATE_INDEX=true \
		TMP_DIR=/tmp
        
:question: **How many of the reads in `NA12878.bam` were marked as PCR duplicates? What might happen
in our subsequent variant calling if we failed to remove these reads?**

<!-- 95, we might call the wrong variants due to the appearance that there's lots of evidence for them -->

We also need to add read groups. These are somewhat hard to define (see discussions [here](https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups)), 
but are **required** by GATK. To give a very simple explanation, read groups can be used to identify which samples were run on which sequencing runs, which
help GATK to mitigate technical artifacts due to inconsistencies between runs.

In our case, we're just going to mark all reads as coming from the same read group. Again, we can do this using Picard:

	java -jar $PICARD AddOrReplaceReadGroups \
		I=NA12878.markduplicates.bam \
		O=NA12878.markduplicates.rg.bam \
		RGID=1 RGLB=lib1 RGPL=illumina \
		RGPU=unit1 RGSM=NA12878 \
		SORT_ORDER=coordinate CREATE_INDEX=True

	java -jar $PICARD AddOrReplaceReadGroups \
		I=NA12891.markduplicates.bam \
		O=NA12891.markduplicates.rg.bam \
		RGID=1 RGLB=lib1 RGPL=illumina \
		RGPU=unit1 RGSM=NA12891 \
		SORT_ORDER=coordinate CREATE_INDEX=True

	java -jar $PICARD AddOrReplaceReadGroups \
		I=NA12892.markduplicates.bam \
		O=NA12892.markduplicates.rg.bam \
		RGID=1 RGLB=lib1 RGPL=illumina \
		RGPU=unit1 RGSM=NA12892 \
		SORT_ORDER=coordinate CREATE_INDEX=True

Finally, we perform base quality recalibration [BQSR](https://gatkforums.broadinstitute.org/gatk/discussion/44/base-quality-score-recalibration-bqsr). 
From the authors of GATK:

> In a nutshell, [BQSR] is a data pre-processing step that detects systematic errors made by the sequencer when it estimates the quality score of each base call.

<!--java -jar $PICARD CreateSequenceDictionary R=grch37.fa O=grch37.dict # not necessary?-->
	
	java -jar $GATK -T BaseRecalibrator -R grch37.fa \
		-I NA12878.markduplicates.rg.bam \
		-knownSites knownSites.vcf -o NA12878.recal_data.table

	java -jar $GATK -T BaseRecalibrator -R grch37.fa \
		-I NA12891.markduplicates.rg.bam \
		-knownSites knownSites.vcf -o NA12891.recal_data.table

	java -jar $GATK -T BaseRecalibrator -R grch37.fa \
		-I NA12892.markduplicates.rg.bam \
		-knownSites knownSites.vcf -o NA12892.recal_data.table

We then just print out a recalibrated BAM like so:

	java -jar $GATK -T PrintReads -R grch37.fa \
		-I NA12878.markduplicates.rg.bam \
		-BQSR NA12878.recal_data.table -o NA12878.markduplicates.rg.bqsr.bam

	java -jar $GATK -T PrintReads -R grch37.fa \
		-I NA12891.markduplicates.rg.bam \
		-BQSR NA12891.recal_data.table -o NA12891.markduplicates.rg.bqsr.bam

	java -jar $GATK -T PrintReads -R grch37.fa \
		-I NA12892.markduplicates.rg.bam \
		-BQSR NA12892.recal_data.table -o NA12892.markduplicates.rg.bqsr.bam

:question: **Now that you've recalibrated your alignment files, which file contains the calibrated
alignment for NA12891?**

<!-- NA12891.markduplicates.rg.bqsr.bam -->
        
## Calling variants in individual samples

We are going to perform variant calling in two stages: **(1)** in individual samples and **(2)** jointly across samples (n = 3). 
To perform variant calling in each of the samples, we will use the HaplotypeCaller tool in GATK:

	java -Xmx2g -jar $GATK -R grch37.fa -T HaplotypeCaller \
		-I NA12878.markduplicates.rg.bqsr.bam \
		--emitRefConfidence GVCF -o NA12878.g.vcf

	java -Xmx2g -jar $GATK -R grch37.fa -T HaplotypeCaller \
		-I NA12891.markduplicates.rg.bqsr.bam \
		--emitRefConfidence GVCF -o NA12891.g.vcf

	java -Xmx2g -jar $GATK -R grch37.fa -T HaplotypeCaller \
		-I NA12892.markduplicates.rg.bqsr.bam \
		--emitRefConfidence GVCF -o NA12892.g.vcf

**NOTE**: It takes a few minutes to run each of these.

## Joint calling across samples

Calling variants in a single individual is _not_ always ideal, however - we have more confidence in observing a particular variant 
at a given site if we see it in other individuals (because we have greater confidence that it is a real site rather than some 
artifact or a sequencing error). Here we call variants in all samples together:

	java -Xmx2g -jar $GATK -T GenotypeGVCFs -R grch37.fa \
		--variant NA12878.g.vcf \
		--variant NA12891.g.vcf  \
		--variant NA12892.g.vcf \
		-o raw_variants.vcf

After you call variants, take a look at the [VCF](https://en.wikipedia.org/wiki/Variant_Call_Format) file that was produced:

	less -S raw_variants.vcf

:question: **How many variants are in the joint VCF file?** (**Hint**: Run vcftools: `vcftools --vcf raw_variants.vcf` and look at the output.) 

<!-- 201 sites -->

## Genotype phasing - what is it?

The DNA sequencing data you have been working with is pretty interesting because the three individuals are related: NA12878 is the daughter of NA12891 and NA12892. Because of this, we can try to _infer_ which parent particular variants come from (this is called _phasing_). Phasing is useful because it helps us link multiple variants on the same haplotype together (otherwise we don't know what alleles are inherited together).

Here we will be _phasing by transmission_ with GATK:

	java -jar $GATK -T PhaseByTransmission -R grch37.fa \
		-V raw_variants.vcf \
		-ped family.ped -o phased_variants.vcf

We now have a bunch of variants that are phased in the trio...so what next? Typically we need to refine our variant calls using _variant quality score recalibration (VQSR)_, but that is impossible for us to with this data because we don't have enough sites to work with. [VQSR](http://gatkforums.broadinstitute.org/gatk/discussion/39/variant-quality-score-recalibration-vqsr) is a really interesting topic to learn more about.

## Filtering our low-quality variants

Given that we are skipping VQSR, there are still some things we can do to clean up our variant calls - for instance, we can remove variants where there weren't many supporting reads, or the reads that mapped to that locus had generally poor mapping qualities, etc. 

	java -jar $GATK -T VariantFiltration -R grch37.fa \
		-V phased_variants.vcf \
		--filterExpression "DP < 10 || QD < 2.0 || FS > 60.0 || MQ < 40.0" \
		--filterName "BIOS201_FILTER" -o flagged_snps.vcf 

:question: **What do the above filters do? Which would we change if we wanted to require greater read depth?**
(**Hint**: Look in the header section of the VCF file.)

<!-- 
DP filters out sites with low depth; 
QD filters out variants with low confidence; 
FS filters out variants with extreme strand bias; 
MQ filters out sites where the supporting reads had generally poor mapping qualities 
-->

In the previous step, we produced a VCF of variants that we would like to remove from our original VCF. We can remove them like so:

	java -jar $GATK -R grch37.fa -T SelectVariants \
		-V flagged_snps.vcf \
		-o filtered_snps.vcf \
		-selectType SNP -ef \
		--restrictAllelesTo BIALLELIC

Look over this file in order to answer the following questions:

:question: **How many variants are in the filtered VCF file?** (Hint: Use `vcftools` with the `--vcf` parameter, as described above.)   
:question: **For the variant at 17:41204377, what are the ref allele and alt alleles?**   
:question: **For the variant at 17:41204377, how many of the individuals were heterozygous?**   
:question: **For the variant at 17:41204377, how many reads supported the ref allele for NA12878? the alt allele?**   

<!-- 134 variants -->
<!-- A>G -->
<!-- 2 -->  
<!-- 28, 24 -->

## Compare variants to platinum genomes

The last thing we are going to do is compare our results to a [gold standard](https://www.illumina.com/platinumgenomes.html) - these three individuals have been extensively sequenced, and their genotypes in the region we have been investigating are well known.

One way to see how well we are doing is to compare the original VCF and the filtered VCF to the gold standard. We are interested in how many calls are _discordant_ between the two:

	vcftools --vcf raw_variants.vcf \
		--diff platinum_trio.vcf \
		--diff-indv-discordance  --out raw_to_platinum_comparison

	vcftools --vcf filtered_snps.vcf \
		--diff platinum_trio.vcf \
		--diff-indv-discordance  --out filtered_to_platinum_comparison

	less raw_to_platinum_comparison.diff.indv 
	less filtered_to_platinum_comparison.diff.indv 

:question: **How many variants were called for NA12892 in the raw VCF?**   
:question: **How many of those were wrong (compared to the platinum genomes)?**    
:question: **How did the results change after filtering?**   

<!-- 196 -->
<!-- 4 -->
<!-- There was no measured discordance after filtering -->

Congratulations! You've aligned a genome, called variants, and filtered low-quality calls.

## Bonus activity

(Optional, on your own time):

If you want more practice calling variants in DNA sequencing data, download the 1000 Genomes
FASTQ files for your favorite samples, publicly available at ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/.
Descriptions of the ancestry for each sample are given [here](http://www.internationalgenome.org/data-portal/sample).
You can download the files easily using the command

`wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/<person_id>/sequence_read/<filename>`

(Better yet, try this with your own FASTQ files from whatever sequencing project you're 
interested in).

NOTE: This may take many hours to run, since the files supplied for this tutorial were 
intentionally shortened to improve running times.
