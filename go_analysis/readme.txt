two steps are needed to generate a GO plot

1. go classification: assign genes into GO categories.
You can use one of the following perl scripts : 
	go_chenwh_ver5_lite.pl <- preferred for simpler tasks
	go_chenwh_ver5.pl
	

at least four parameters are needed to run either of the scripts, they are:
	-i : gene to GO file; see 'example_gene2go_file.txt' for example
	-g : GO ontology file in .obo format; download from here : http://geneontology.org/page/download-ontology
	-p : output file for plot
	-o : not important 
	
	an optional parameter -level specifies the GO level with default value of 2. the larger the number, the more groups you will have.
	
2. plot: generate plot in svg format
The script for this step is : graph.pl:
	perl graph.ph infile output.svg
	
here 'infile' is the output file specified by -p.


