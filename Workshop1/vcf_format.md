### VCF

Variant call files contain information on variants observed (relative to a reference genome) in one or more samples. An short example is shown here:

	#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA12878 NA12891 NA12892
	17      41196408        .       G       A       1349.16 PASS    AC=2;AF=0.333;AN=6;BaseQRankSum=-6.450e-01;ClippingRankSum=0.00;DP=121;ExcessHet=3.9794;FS=0.881;MLEAC=2;MLEAF=0.333;MQ=60.00;MQRankSum=0.00;QD=16.06;ReadPosRankSum=-1.180e+00;SOR=0.909  GT:AD:DP:GQ:PL:TP       1|0:17,24:41:99:623,0,430:96    0|0:35,0:35:96:0,96,1440:96     1|0:18,25:43:99:758,0,510:96
