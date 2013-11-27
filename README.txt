WiggleTools

Copyright Daniel Zerbino, 2013.
zerbino@ebi.ac.uk

This library parses wiggle files and executes various operations on them streaming through lazy evaluators.

Installation:
	- Install Jim Kent's source code
	- In Jim's git repo, checkout branch bigWigCat, re-compile the library and the wigToBigWig and bigWigCat utilities (You can revert to the master branch after that)
	- Define the environment variable KENT_SRC to point to the /path/to/kent/src/ directory
	- Install the GNU scientific library (GSL): http://www.gnu.org/software/gsl/
	- In this directory type 'make'
	- In the same directory, type 'make test'
	- The binary file wiggletools should be produced in the ./bin directory. 
	- The underlying library and header files are in the ./lib and ./inc directories respectively.
	- A PDF manual should appear in doc/

Running the executable:
 
Inputs:
	The program takes in Wig, BigWig, BedGraph, Bed, BigBed and Bam files, which are distinguished thanks to their suffix (.wig, .bw, .bg, .bed, .bb, .bam respectively).
	Note that wiggletools assumes that every bam file has an index .bai file next to it.

Outputs:
	The program outputs a wiggle file in stdout unless the output is squashed with the 'do' command;

Command line:
	wiggletools --help
	wiggletools program

Program grammar:
	program 	= (iterator) | do (iterator) | (statistic) | (extraction)
	statistic 	= AUC (iterator) | mean (iterator) | variance (iterator) | pearson (iterator) (iterator) | isZero (iterator)
	extraction 	= profile (iterator) (iterator) | profiles (iterator) (iterator)
	                  | apply (out_filename) (statistic) (bed_file) (iterator)
	iterator 	= (filename) | (unary_operator) (iterator) | (binary_operator) (iterator) (iterator) | (reducer) (multiplex) | (setComparison) (multiplex) (multiplex)
	unary_operator 	= unit | stdout | write (filename.wig) | smooth (int)
	binary_operator = diff | ratio
	multiplex 	= (filename_list) | map (unary_operator) (multiplex)
	reducer 	= cat | sum | product | mean | var | stddev | median | min | max
	setComparison 	= ttest | wilcoxon
	filename_list 	= (filename) : | (filename) (filename_list)
	filename 	= *.wig | *.bw | *.bed | *.bb | *.bg | *.bam
