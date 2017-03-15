# Workshop 1: Variant calling from DNA sequencing data

## Answers

**Of the 4 reads you just looked at, did you see any obvious problems?**

>The first read in the R2 file is clearly a *bad* read 

**What happens to the quality of sequenced bases in later cycles?**   

>The base qualities drop off in the latter cycles

**Do you see some adapter contamination in the data?**   

>A little bit, but doesn't fail FastQC thresholds

**How many read pairs were there before/after trimming?**

>Before trimming: 10,844
>After trimming: 10,336

**What percentage of reads for NA12878 were successfully mapped?**

>93.86% of the reads

**How many variants are in the joint VCF file?** 

>201 variants

**How many sites were not phased?** 

>10 variants

**What do the above filters do?**

>DP filters out sites with low depth; 
>QD filters out variants with low confidence; 
>FS filters out variants with extreme strand bias; 
>MQ filters out sites where the supporting reads had generally poor mapping qualities 

**How many variants are in the filtered VCF file?**  

>134 variants

**For the variant at 17:41204377, what are the ref allele and alt alleles?**  

>chr17:g.41204377A>G

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