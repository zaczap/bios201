# Workshop 1: Variant calling from DNA sequencing data

## Introduction

We've provided you with raw sequencing data (FASTQ format) for 3 individuals that were sequenced as part of the 1000 Genomes project: NA12878, NA12891, and NA12892. Each of these samples has two FASTQ files because the sequencing was done with paired-end reads (R1 and R2), and we've subsetted the FASTQ files to reads that pretty much exclusively from the _BRCA1_ locus on chromosome 17 (17:41196312-41277500 in GRCh37).

- NA12878: NA12878_R1.fastq, NA12878_R2.fastq
- NA12891: NA12891_R1.qc.trimmed.fastq, NA12891_R2.qc.trimmed.fastq
- NA12892: NA12892_R1.qc.trimmed.fastq, NA12892_R2.qc.trimmed.fastq

Our goal is to take these reads, align them to the human reference (GRCh37), call variants following the GATK Best Practices pipeline as best as possible, and inspect variation in the _BRCA1_ locus in these individuals.

In the following workshop, you will learn to:

- Run `fastqc` on FASTQ files to assess sequencing quality
- Run `cutadapt` to remove adapters, trim reads, and filter short reads
- Index a reference genome with `bwa`
- Align paired-end reads to your indexed genome using `bwa mem`
- Inspect alignment results with `samtools flagstat`
- Mark PCR duplicates in the aligned reads (Picard tools)
- Perform base quality recalibration on the aligned reads (GATK)
- Call variants in individual samples (GATK)
- Jointly call variants across samples (GATK)

## Setup

First you will have to connect to the server we have set up for you:

	ssh <sunetID>@corn.stanford.edu

Once connected, run:

	echo $SHELL

If the response is `tcsh`, download and run the following file:

	wget http://www.stanford.edu/~zappala/bios201/setup_tcsh.sh
	source setup_tcsh.sh

If the response is `bash`, download and run the following file:

	wget http://www.stanford.edu/~zappala/bios201/setup_bash.sh
	source setup_bash.sh

Finally, create a directory to work in and `cd` into it:

	mkdir workshop1
	cd workshop1
	wget http://www.stanford.edu/~zappala/bios201/workshop1.gz
	gunzip workshop1

## Sequencing quality control with FASTQC

When working with sequencing data, you want to make sure your reads do not have systemic biases, adapter contamination, and are of generally good quality. The tool FASTQC can help assess these artifacts in your sequencing files, and is fairly easy to run:

	fastqc --outdir <output_dir> --format <format> --adapters <adapter_file> <r1_fastq> <r2_fastq>

	# <output_dir> is the name of a folder where you want the results written
	# <format> is the format of the sequencing reads (optional)
	# <adapter_file> is the name of a tab-delimited file that lists the adapters for each sample
	# <r1_fastq> is the first set of reads
	# <r2_fastq> is the second set of reads 

**NOTE:** We need to create __<output_dir>__ before we can supply it to `fastqc`.

**NOTE:** We are only going to perform QC on NA12878; the other samples are already QC'ed and trimmed.

Let's take a quick look at the `adapters.txt` file:

	# Inspect adapters.txt
	cat adapters.txt

And now let's run FASTQC:

	fastqc --outdir fastqc --format fastq --adapters adapters.txt NA12878_R1.fastq NA12878_R2.fastq

This should run pretty quickly. When it's done, let's go look at the results:

	ls fastqc

You can see a `.html` file and a `.zip` file for each FASTQ file - if we download the `.html` files, we can look at them on our local computer to inspect the results from FASTQC. Normally you would need to download these files to your machine and open them, but we have made them available at the following URLs:

- [NA12878_R1.fastq][http://www.stanford.edu/~zappala/bios201/NA12878_R1.html]
- [NA12878_R2.fastq][http://www.stanford.edu/~zappala/bios201/NA12878_R2.html]

## Removing adapters, trimming reads, and filtering

From the FASTQC results, you can see a little bit of adapter contamination as well a tell-tale drop off in sequencing quality towards the end of the reads. We can remove adapter contamination and trim these low-quality bases using the tool `cutadapt`:

	cutadapt -q <minimum_quality> --minimum-length <minimum_length>
		-a <adapter on first read>
		-A <adapter on second read>
		-o <r1_trimmed_fastq> -p <r2_trimmed_fastq> <r1_fastq> <r2_fastq>

	# <minimum_quality> is the minimum base quality to allow (trimmed otherwise)
	# <minimum_length> is the minimum length both pairs of reads need to be to be included
	# <adapter on first read> for us: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
	# <adapter on second read> for us: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
	# <r1_trimmed_fastq> the name of the output file for the first read file 
	# <r2_trimmed_fastq> the name of the output file for the second read file 
	# <r1_fastq> is the first set of reads
	# <r2_fastq> is the second set of reads 

So let's run the actual command:

	cutadapt -q 10 --minimum-length 30 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o NA12878_R1.qc.trimmed.fastq -p NA12878_R2.qc.trimmed.fastq NA12878_R1.fastq NA12878_R2.fastq

Now we've removed adapters and trimmed low quality bases; we've also filtered out very short reads which would probably not align to the genome anyway.

## Aligning paired-end DNA sequencing data to the human genome

We would like to take our reads (in FASTQ format) and figure out where they originated from in the human genome - this __alignment__ is accomplished using a tool called `bwa` (which stands for the Burrows-Wheeler aligner). We first need to __index__ the genome before it can be easily aligned to:

	bwa index grch37.fa

**NOTE:** The genome FASTA file we have given you only includes chromosome 17.

Now we can align reads using the following command:

	bwa mem <genome> <r1_fastq> <r2_fastq> | samtools sort --output-fmt BAM > <bam_file>

	# <genome> is the genome we are aligning to - it must be indexed before it can be used!
	# <r1_fastq> the first set of trimmed reads
	# <r2_fastq> the second set of trimmed reads
	# <bam_file> the alignment file where we will save results

Let's do that for each of our samples (normally this can take some time, but we have only given you a small set of reads):

	bwa mem grch37.fa NA12878_R1.qc.trimmed.fastq NA12878_R2.qc.trimmed.fastq | samtools sort --output-fmt BAM > NA12878.bam

	bwa mem grch37.fa NA12891_R1.qc.trimmed.fastq NA12891_R2.qc.trimmed.fastq | samtools sort --output-fmt BAM > NA12891.bam

	bwa mem grch37.fa NA12892_R1.qc.trimmed.fastq NA12892_R2.qc.trimmed.fastq | samtools sort --output-fmt BAM > NA12892.bam

We can see how well the alignments were by running the following commands:

	samtools flagstat NA12878.bam
	samtools flagstat NA12891.bam
	samtools flagstat NA12892.bam

## Marking PCR duplicates and performing base quality recalibration

First we mark duplicates:

	java -Xmx2g -jar $PICARD MarkDuplicates INPUT=NA12878.bam OUTPUT=NA12878.markduplicates.bam METRICS_FILE=NA12878.markduplicates.metrics.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 CREATE_INDEX=true TMP_DIR=/tmp

	java -Xmx2g -jar $PICARD MarkDuplicates INPUT=NA12891.bam OUTPUT=NA12891.markduplicates.bam METRICS_FILE=NA12891.markduplicates.metrics.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 CREATE_INDEX=true TMP_DIR=/tmp

	java -Xmx2g -jar $PICARD MarkDuplicates INPUT=NA12892.bam OUTPUT=NA12892.markduplicates.bam METRICS_FILE=NA12892.markduplicates.metrics.txt OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 CREATE_INDEX=true TMP_DIR=/tmp

Then we add read groups, because those were missing but are required by GATK:

	java -jar $PICARD AddOrReplaceReadGroups I=NA12878.markduplicates.bam O=NA12878.markduplicates.rg.bam RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=NA12878

	java -jar $PICARD AddOrReplaceReadGroups I=NA12891.markduplicates.bam O=NA12891.markduplicates.rg.bam RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=NA12891

	java -jar $PICARD AddOrReplaceReadGroups I=NA12892.markduplicates.bam O=NA12892.markduplicates.rg.bam RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=NA12892

Then we index these BAMs:

	samtools index NA12878.markduplicates.rg.bam
	samtools index NA12891.markduplicates.rg.bam
	samtools index NA12892.markduplicates.rg.bam

And we perform base quality recalibration:

	java -jar $PICARD CreateSequenceDictionary R=grch37.fa O=grch37.dict
	
	java -jar $GATK -T BaseRecalibrator -R grch37.fa -I NA12878.markduplicates.rg.bam -knownSites knownSites.vcf -o NA12878.recal_data.table

	java -jar $GATK -T BaseRecalibrator -R grch37.fa -I NA12891.markduplicates.rg.bam -knownSites knownSites.vcf -o NA12891.recal_data.table

	java -jar $GATK -T BaseRecalibrator -R grch37.fa -I NA12892.markduplicates.rg.bam -knownSites knownSites.vcf -o NA12892.recal_data.table

	java -jar $GATK -T PrintReads -R grch37.fa -I NA12878.markduplicates.rg.bam -BQSR NA12878.recal_data.table -o NA12878.markduplicates.rg.bqsr.bam

	java -jar $GATK -T PrintReads -R grch37.fa -I NA12891.markduplicates.rg.bam -BQSR NA12891.recal_data.table -o NA12891.markduplicates.rg.bqsr.bam

	java -jar $GATK -T PrintReads -R grch37.fa -I NA12892.markduplicates.rg.bam -BQSR NA12892.recal_data.table -o NA12892.markduplicates.rg.bqsr.bam

## Calling variants in individual samples

	java -Xmx2g -jar $GATK -R grch37.fa -T HaplotypeCaller -I NA12878.markduplicates.rg.bqsr.bam --emitRefConfidence GVCF -o NA12878.g.vcf

	java -Xmx2g -jar $GATK -R grch37.fa -T HaplotypeCaller -I NA12891.markduplicates.rg.bqsr.bam --emitRefConfidence GVCF -o NA12891.g.vcf

	java -Xmx2g -jar $GATK -R grch37.fa -T HaplotypeCaller -I NA12892.markduplicates.rg.bqsr.bam --emitRefConfidence GVCF -o NA12892.g.vcf

## Joint calling across samples

	time java -Xmx2g -jar $GATK -T GenotypeGVCFs -R grch37.fa --variant NA12878.g.vcf --variant NA12891.g.vcf  --variant NA12892.g.vcf -o raw_variants.vcf

## Phasing within the trio

	java -jar $GATK -T PhaseByTransmission -R grch37.fa -V raw_variants.vcf -ped family.ped -o phased_variants.vcf

## Filter variants

	java -jar $GATK -T VariantFiltration -R grch37.fa -V phased_variants.vcf --filterExpression "DP < 10 || QD < 2.0 || FS > 60.0 || MQ < 40.0" --filterName "BIOS201_FILTER" -o flagged_snps.vcf 

	java -jar $GATK -R grch37.fa -T SelectVariants -V flagged_snps.vcf -o filtered_snps.vcf -o -selectType SNP -ef --restrictAllelesTo BIALLELIC

## Compare variants to platinum genomes

	vcftools --vcf filtered_snps.vcf --diff platinum_trio.vcf --diff-indv-discordance  --out compare

	more compare.diff.indv 
<!-- INDV    N_COMMON_CALLED N_DISCORD       DISCORDANCE
NA12878 134     0       0
NA12891 134     0       0
NA12892 134     0       0 -->
