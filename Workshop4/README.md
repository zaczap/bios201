# Workshop 4: Epigenomics data -- Visualization, Peak Calling, Quality Control, and Motif Finding

In today's workshop we will be analyzing some ATAC-seq data as an example of epigenomics data analysis. This analysis is fairly similar to what one might do with ChIP-seq or DNAse-seq data.  

We will be working with ATAC-seq data from two hematopoeitic cell types -- Natural Killer (NK) cells and Hematopoietic Stem Cells (HSC). For each cell type, we will use data from two different patients. This data is just a small subset of the data from this [paper](http://www.nature.com/ng/journal/v48/n10/abs/ng.3646.html).

We have already aligned and filtered the reads -- we will be starting this analysis with bam files. Additionally, only reads mapping to chr4 have been retained, in the interest of speeding up all the analyses.   

# Setup

Follow the setup steps outlined in the first workshop.
```
ssh <sunetID>@corn.stanford.edu
# e.g. ssh zappala@corn.stanford.edu
```

Once connected, run:
```
echo $SHELL
```

If the response is `tcsh`, run the following command:
```
source /afs/ir/class/bios201/setup/setup_tcsh.sh 
```

If the response is `bash`, run the following command:
```
source /afs/ir/class/bios201/setup/setup_bash.sh 
```

Then copy the workshop materials and `cd` into the copied directory:
```
cp -r /afs/ir/class/bios201/workshop4/ .
cd workshop4
```

# Making coverage track

For visualizing ATAC-seq data, it can be helpful to make a "coverage" track showing how many fragments map to a particular region of the genome. We will be using the `bamCoverage` tool from the deepTools suite to make coverage tracks. To learn about all the options, check out the [bamCoverage documenation](https://deeptools.readthedocs.io/en/latest/content/tools/bamCoverage.html).

We will be using the following set of options:
```
# Example format for call
bamCoverage --bam <input.bam> -o <output> --binSize 100 --normalizeUsingRPKM --extendReads
```

To do this for all four of our bam files:
```
bamCoverage --bam Donor2596-NK.chr4.bam -o Donor2596-NK.chr4.bw --binSize 100 --normalizeUsingRPKM --extendReads

bamCoverage --bam Donor5483-NK.chr4.bam -o Donor5483-NK.chr4.bw --binSize 100 --normalizeUsingRPKM --extendReads

bamCoverage --bam Donor4983-HSC.chr4.bam -o Donor4983-HSC.chr4.bw --binSize 100 --normalizeUsingRPKM --extendReads

bamCoverage --bam Donor7256-HSC.chr4.bam -o Donor7256-HSC.chr4.bw --binSize 100 --normalizeUsingRPKM --extendReads

```

# Using the UCSC Genome Browser to visualize the data

We are going to use the UCSC Genome Browser to visualize our ATAC-seq data. First, navigate to the [UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgGateway).  For genome build, select "Feb. 2009 GCRCh37/hg19".  For region, type in `chr4:105,079,645-105,815,457`.

![UCSC input](images/ucsc_entry.png)

Then click GO.  This will take you to a page that should look like this:

![UCSC input](images/ucsc_basic.png)

Take a minute to scroll down the page and see some of the different kinds of data that are available to display.  Try adding some additional information and changing the display settings for some of the default information.

Now we're going to add some of our own data.  Click on the 'add custom tracks' gray button that is in the middle below the tracks.  

![UCSC input](images/add_custom_track.png)

We are going to give a URL pointing to our data.  To do that, we first have to copy the data to somewhere it is publically available. Copy all your bigwig files to your WWW folder:

``` 
cp *.bw ~/WWW/
```

Now if you go to http://web.stanford.edu/~sunet/, fillling in sunet with your own sunet id, you will see a list of files.  Now enter http://web.stanford.edu/~sunet/Donor2596-NK.chr4.bw (again substituting sunet with your own id) into the `Paste URLs or data` box.  Then click `submit`. You will get taken to another screen, where you should click the 'go` button next to view in `Genome Browser`. You should now see some data!

Explore different options for displaying the data, selected from the dropdown box int he custom tracks section:

![UCSC input](images/track_options.png)

Use the same procedure to add tracks for the other three samples.  You can add them all at once by pasting each link into the `Paste URLs or data` box.

:question: **1.1 


# Peak calling

The first step in our analysis will be to call peaks for each sample.  Peaks are areas of the genome where we have a pileup of signal -- in ATAC-seq these represent "accessible" regions of the genome.  

We will use the tool MACS2 to call peaks. With MACS2 you have the option of calling "narrow" or "broad" peaks. Generally, for TF ChIP-seq, "narrow" peaks are appropriate while for histone modification ChIP-seq "broad" peaks are appropriate. For ATAC-seq, peaks can vary in size and depending on what you want do with the peaks it may make sense to call either broad or narrow peaks.  For this tutorial, we will call narrow peaks. In addition we will use the `--call-summits` option to call a "summit" in each peak-- this represents where the peak has the greatest intensity.  

Here we will use the following MACS2 options:

| option | description | 
| ------ | ------ |
| --treatment | name of the file(s) | 
| --name | name to use as prefix of output | 
| --format | what is file format | 
| --call-summits | also output summit positions |
| --keep-dup | count duplicates? | 
| -p |p value cutoff |
| -g | genome size -- normally we would use hg19 but we have subset the genome to only include chr4 | 
| --nomodel | don't compute model for fragment size |

To learn more about these options and the other available options, read throught the [MACS2 documentation](https://github.com/taoliu/MACS)

```
macs2 callpeak --treatment Donor2596-NK.chr4.bam --name Donor2596-NK --format BAMPE \
  --nomodel --call-summits --nolambda --keep-dup all -p -0.01 -g 1.7e8
```
