Workshop 2: Aligning and quantifying reads from RNA-seq
=======================================================

Answers
-------
**1.1** Four transcripts with 6,5,6,2 exons, respectively.

**1.2** The gene model looks like (with exons numbered in transcription order):

![gene model SFTPA2](http://web.stanford.edu/class/bios201/SFTPA2.png)

Some of the transcripts have slightly different starting and end
  position for some exons. The first and third transcripts include all 6
  exons. The second transcript is missing exon 2. The fourth transcript
  only includes exons 3 and 4, again with different starting and end positions.

**1.3** The exons are listed in reverse order. This is because the
  transcript is on the reverse strand and, in this GTF file, the exons are
  listed in the order in which they are transcribed.

---------

**2.1** 94.01%

---------

**3.1** This read mapped to chr10:104849444 and its mate mapped to chr10:104849571.

**3.2** The flag is 99, which means that the read is paired, is the first
read in the pair, was mapped in a proper pair (i.e. its mate mapped to the
opposite strand), and its mate mapped to the reverse strand (so this read
mapped to the forward strand). Also, since neither the "not primary alignment"
nor the "supplementary alignment" bits were set, we can infer that this is the
primary alignment.

**3.3** The two reads overlap and span the same exon. The intron is 334
  bases long (both CIGAR strings include N). See next question for a visualization.

---------

**4.1** ![screenshot of IGV](http://web.stanford.edu/class/bios201/screenshot.igv1.png)

**4.2** The vast majority of the reads are pink, which means that the
  first read is mapping to the forward strand. The gene is on the reverse
  strand (look at the arrows in the gene model at the bottom of the
  window). Therefore, we have a stranded library where the strand of the
  second read matches the strand of the gene.

**4.3** The curved lines represent splicing events and the numbers on
  those lines indicate the number of reads supporting the splicing of the
  given exons together. Reads support a given splicing event by spanning
  an exon-exon junction.

---------

**5.1** 92.4% (53.9% to coding, 38.5% to UTRs)

**5.2** You might not have RNA. Alternatively, if the percentage of reads
  mapping to coding bases is moderate, you might have RNA with DNA
  contamination.

---------

**6.1** Norm7

**6.2** In Norm7, these two genes are barely expressed compared to the
  other samples. This suggests thatt Norm7 is not actually a lung sample.