### SAM/BAM

When we align reads to a reference genome, we generated a _sequence alignment_ which we store in either plaintext (SAM) or binary (BAM) - these two file types are the same except one is human readable and larger on disk than the other other. An example line is given here:
	
	HWI-D00360:6:H81VLADXX:1:1101:13347:2151        2177    17      188091  5       73H59M16H       =       74713352        74525262        AAAAATATAAAACTTAGCTGGGCATGGTAGCACACACCTGTAGTCTCAGCTACTCAGGA     CEEEEBCCDEEC>::C:@>>??AC32144::@4?<@5<>>+:4+(4(>((4>>34@49A     NM:i:4  MD:Z:6T21G5T14A9        AS:i:39 XS:i:34 SA:Z:17,4264587,+,101M47S,5,7;  XA:Z:17,+73713822,102S34M12S,0;17,-19146648,12S43M93S,2;
