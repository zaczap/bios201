### FASTQ

The fundamental next-generation sequencing data format is FASTQ - this is the raw data produced from a sequencing run on most common platforms. 

Each read in a FASTQ file gets 4 lines - for example, the following 4 lines describe one read:

	@HWI-D00360:6:H81VLADXX:1:1101:1245:2105 1:N:0:CGATGT
	TTTTTTTCTGAGGCAAGTCCCACTCTCTTGCCCAGGCTGGAGTGCAGTAGTGTGACCTCGGCTCACTGCAACCTCCGTCCCCCAGGTTCAAGTGATTCTCCTGCCTTANCCTCCCAAGTNNCTGNGNTCACAGGNNNCCACCATCATG
	+
	@@@DDDDDD?DBFFF?GCFGEHIE@GGIIEGGIICBFEC<CGG;=B<FCDCDFE@DG:DA<EEEEE;(;>CCBCCB:9<BBBB9?7::ACCCD>@DDDCCC@@@BCC@#++<@BBCC<C##++8#+#+++8<C1###+++8@B@CAC:

The first line, which begins with a `@` character, is a _header_ line that describes where the read came from: what machine, flow cell, etc. What goes here depends on what sequencing machine was used, but it is typically all of the information necessary to figure out where the read came from and to uniquely identify it amongst all other reads.

The second line, which is arguably the most interesting line, contains the bases that were sequenced - it looks like DNA.

The third line, which begins with `+` character, is often blank, and often just repeats the information in the header line otherwise.

Finally, the fourth line contains the _quality scores_ for each of the bases that were sequenced. The encoding scheme is a little confusing, however. You can read all about it [here](https://en.wikipedia.org/wiki/Phred_quality_score) and [here](http://www.drive5.com/usearch/manual/quality_score.html), but basically, each character has a numeric value attached to it that is translated into an per-base accuracy. For example:

	!	encodes a Phred quality score of 0 which ~= 0% accuracy
	"	encodes a Phred quality score of 1 which ~= 21% inaccuracy
	#	encodes a Phred quality score of 2 which ~= 37% inaccuracy
	...
	A	encodes a Phred quality score of 32 which ~= 99.937% inaccuracy
	B	encodes a Phred quality score of 33 which ~= 99.95% inaccuracy
	C	encodes a Phred quality score of 34 which ~= 99.96% inaccuracy
	D	encodes a Phred quality score of 35 which ~= 99.968% inaccuracy
	...

The most important part is that you understand each base is sequenced with varying quality - intutitively you probably want to throw out (or trim) bad bases before continuing your experiment!