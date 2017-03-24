# Workshop 4 answers

:question: What is the most common fragment length?

The mode of the distribution is at approximately 60 bp.  

:question: What does the peak at around 200 bp likely represent?

Reads spanning nucleosomes.

:question: What is the max coverage for each sample in the chr4:105,079,645-105,815,457 window?



:question: Try adding the insertion track as well -- how does it compare to the coverage track for the same sample?

Insertion track is much less smooth! Shows per base insertion points at base pair resolution rather than fragment coverage in 100 bp windows.

:question: What is the approximate width of the peak at the TSS?

About 500 base pairs -- of course depends somewhat about how you define width.    

:question: What kinds of files were output by MACS2? What is the difference between the different formats?

narrowPeak, xls, and summits bed.  The summits bed file has one line per summit. The narrowPeak file has each peak, with one column indicating the summit position. The xls file is similar to narrowPeak but includes some header information.

:question: How many peaks are there for each sample? [hint: Remember the wc command...]

```
wc -l *.narrowPeak
```

   1961 Donor2596-NK_peaks.narrowPeak
   9652 Donor4983-HSC_peaks.narrowPeak
  15013 Donor5483-NK_peaks.narrowPeak
  10381 Donor7256-HSC_peaks.narrowPeak
  
 Note however that some of these are the same peak, but with different summits. To get just the unique peaks, you could use the bedtools merge command.   

:question: How many peaks are there in the window: chr4:105,410,415-105,419,498? How many summits?

Three peaks and 6 summits.  A peak can have multiple called summits.  

:question: How many peaks are in the merged peak files?

```
bedtools merge -i all_peaks.bed | wc -l
```

22075 peaks

:question: What is the the most enriched known motif in the NK peaks? In the HSC peaks?

In NK:  ETS
In HSC: CTCF

