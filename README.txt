The following is instructions for the offspotter remote utility, not the guide generation.

INSTRUCTIONS:
-------------
Please see https://cm.jefferson.edu/downloads/off-spotter-help/ for more details. Contact us if
you have further questions.

Summary:
--------
You can use this utility to submit remote requests on Off-Spotter (cm.jefferson.edu/Off-Spotter).
You should use this if you want to submit a sequence that's larger than 1500 bp or more than 20
gRNAs. It accepts one input file but each input file can contain several sequences. This utility
will submit the gRNAs in batches of 5. Please don't run more than one such requests at a time.


Requirements:
-------------
1. Java 1.6 or greater must be used (http://www.java.com/en/download/index.jsp to download if
you don't currently have it installed).  
	- issuing the command 'java -version' should tell you your current Java version


How to use:
-----------

1.  Unzip/extract the included zip file into a folder

2.  Open the extracted 'Parameters.properties' file and change desired parameters as needed.
See below and for more details see help page given on top.
	- input: Should be the name of the input file. The input file must be located in the
	same directory as the .class file. This file can have either multiple sequences or
	multiple gRNAs. The gRNAs should be one per line, and the PAM sequences should be omitted. 
	- genome: The available options are hg19 (or GRCh37) for the old human genome, hg38
	(or GRCh38) for the new human genome, mm10 (or GRCm38) for the new mouse genome, and w303
	for yeast.
	- mismatches: The number of different bases between the gRNA and the actual genomic hits.
	- annotation: Whether the results will be annotated using ENSEMBL annotation. Can be on
	or off.
	- PAM: Available choices are NGG, NAG, NNNNACA, and NNGRRT (where R will be automatically
	restricted to A or T).
	- seed: Select the seed schema of your choice. Use - to indicate that the position can
	have mismatches and x to indicate it can not. Default is ---------------xxxxx.
	- output: The utility will create a folder with that name. That folder will hold all the
	results. This option allows you to run Off-Spotter remote from the same folder multiple
	times for different parameters and keep all results easily. Note that, if the folder name
	already exists the folder will be overwritten.
	- Notes:
		- blank lines will be ignored
		- If the input file contains multiple sequences then you will get results for each
		sequence individually.
		- there is no limit to the number of gRNAs or the size of sequences you can have
		in the input file

3.  Run the Off-Spotter_remote utility by issuing the command 'java OffSpotter_remote'
	- You must run this command from the directory you extracted the files to in #1 (OffSpotter_remote.class must be in your current working directory)
	- Note: running the program will fetch results from the server in multiple transactions (5 gRNAs at a time).  This is to allow sharing of our server resources with other people running the same tool.

4.  OffSpotter_remote will create an output folder with the name of your choice for the results. If you are giving sequences as input to Off-Spotter, in this folder you'll find a summary file and another folder containing the results per sequence. If your input file contains gRNAs only, in this folder you'll find a second folder with the results.



Considerations:
Please, be kind to your colleagues! Everyone who uses the Off-Spotter server would like to enjoy reasonable response times. You can help achieve this by submitting not more than a few tens of requests at a time.


Thank you for using Off-Spotter!
