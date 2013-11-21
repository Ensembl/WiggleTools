WiggleTools

Copyright Daniel Zerbino, 2013.
zerbino@ebi.ac.uk

This library parses wiggle files and executes various operations on them streaming through lazy evaluators.

Installation:
	- Install Jim Kent's source code
	- Define the environment variable KENT_SRC to point to the /path/to/kent/src/ directory
	- Install the GNU scientific library (GSL): http://www.gnu.org/software/gsl/
	- In this directory type 'make'
	- The binary file wiggletools should be produced in the ./bin directory. The underlying library and header files are in the ./lib and ./inc directories respectively.
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
	program = do (command) | (command)                ## 'do' is short for 'just do it and don't print it out'
	command = (statistic) (iterator) | (iterator)   
	statistic = AUC | mean | variance | pearson    ## A statistic is simply a number computed across the dataset
	iterator = (filename) | (unary_operator) (iterator) | (binary_operator) (iterator) (iterator) | (reducer) (multiplex)
	unary_operator = unit | stdout | write (filename.wig) | smooth (int)
	binary_operator = diff
	multiplex = (filename_list) | map (unary_operator) (multiplex)
	reducer = cat | add | product | mean | var | stddev | median | min | max
	filename_list = (filename) ; | (filename) (filename_list)
	filename = *.wig | *.bw | *.bed | *.bb | *.bg | *.bam
