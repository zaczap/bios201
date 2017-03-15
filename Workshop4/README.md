# Epigenomics Analysis workshop: Peak Calling, Quality Control, and Motif Finding

In today's workshop we will be analyzing some ATAC-seq data as an example of epigenomics data analysis. This analysis is fairly similar to what one might do with ChIP-seq or DNAse-seq data.  

We will be working with ATAC-seq data from two hematopoeitic cell types -- Natural Killer cells and 

We have already aligned and filtered the reads -- we will be starting this analysis with bam files. Additionally, only reads mapping to chr4 have been retained, in the interest of speeding up all the analyses.   

# Setup

Follow the setup steps outlined in the first workshop (see
[here](https://github.com/zaczap/bios201/blob/master/setup.md) for a refresher).  
Then copy the workshop materials and `cd` into the copied directory:
```
cp -r /afs/ir/class/bios201/workshop4/ .
cd workshop4
```

# Peak calling

The first step in our analysis will be to call peaks for each sample.  Peaks are areas of the genome where we have a pileup of signal -- in ATAC-seq these represent "accessible" regions of the genome.  

We will use the tool MACS2 to call peaks. With MACS2 you have the option of calling "narrow" or "broad" peaks. Generally, for TF ChIP-seq, "narrow" peaks are appropriate while for histone modification ChIP-seq "broad" peaks are appropriate. For ATAC-seq, peaks can vary in size and depending on what you want do with the peaks it may make sense to call either broad or narrow peaks.  For this tutorial, we will call narrow peaks. In addition we will use the `--call-summits` option to call a "summit" in each peak-- this reprsents where the peak has the greatest intensity.  

To learn about the various options

Here we will use the following options:

| option | description |
| ------ | ------ |
| --treatment | name of the file |
| --name | name to use as prefix of output |
| --format | what is file format | 
| --call-summits | also output summit positions |
| --keep-dup | count duplicates? |
| -p |p value cutoff |
| -g | genome size -- normally we would use hg19 but we have subset the genome to only include chr4 | 
| --nomodel | don't compute model for fragment size |

```
macs2 callpeak --treatment Donor2596-NK.chr4.bam --name Donor2596-NK --format BAMPE --nomodel --call-summits --nolambda --keep-dup all -p -0.01 -g 1.7e8
```
