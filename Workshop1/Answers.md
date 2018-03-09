# Workshop 1: Variant calling from DNA sequencing data

## Answers

**How many reads are displayed using each of these commands?**

>Each of the `head` command displays 2 reads.

**Of the reads you viewed, which has the highest quality? Which has the lowest?**

>The first read in the R2 file is clearly a *bad* read, based on the quality score and
the fact that most bases are called as N's, meaning the base is unknown.
The second read in that file is quite good, though.

**What happens to the quality of sequenced bases in later cycles?**   

>The base qualities drop off in the latter cycles

**Do you see some adapter contamination in the data?**   

>A little bit, but doesn't fail FastQC thresholds

**How many read pairs were there before/after trimming?**

>Before trimming: 10,844
>After trimming: 10,336

**What percentage of reads for NA12878 were successfully mapped?**

>93.86% of the reads

**How many of the reads in `NA12878.bam` were marked as PCR duplicates?**
>95

** What might happen in our subsequent variant calling if we failed to remove these reads?**
>We might call non-existent variants due to the appearance that there's strong of evidence for them, when
in fact we're just seeing many PCR copies of the same starting fragment.

**Now that you've recalibrated your alignment files, which file contains the calibrated
alignment for NA12891?**
>NA12891.markduplicates.rg.bqsr.bam

**How many variants are in the joint VCF file?** 

>201 variants

**What do the above filters do?**

>DP filters out sites with low depth; 
>QD filters out variants with low confidence; 
>FS filters out variants with extreme strand bias; 
>MQ filters out sites where the supporting reads had generally poor mapping qualities 

**Which would we change if we wanted to require greater read depth?**
>We would increase the DP filter.

**How many variants are in the filtered VCF file?**  

>134 variants

**For the variant at 17:41204377, what are the ref allele and alt alleles?**  

>A>G

**For the variant at 17:41204377, how many of the individuals were heterozygous?**   

>2 heterozygotes

**For the variant at 17:41204377, how many reads supported the ref allele for NA12878? the alt allele?**   

>28 reads supported the ref allele
>24 reads supported the alt allele

**How many variants were called for NA12892 in the raw VCF?**   

>196 variants

**How many of those were wrong (compared to the platinum genomes)?**  

>4 mistakes

**How did the results change after filtering?**   

>After filtering there was no discordance