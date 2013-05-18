WiggleTools

Copyright Daniel Zerbino, 2013.
zerbino@ebi.ac.uk

This library parses wiggle files and executes various operations on them streaming through lazy evaluators.

Installation:
	- Install Jim Kent's source code
	- Define the environment variable KENT_SRC to point to the /path/to/kent/src/ directory
	- In this directory type 'make'
	- The binary file wiggletools should be produced in the same directory. You can add it to your PATH, or move it onto your PATH

Running the executable:
 
Inputs:
	The program takes in Wig and BigWig files, which are distinguished thanks to their suffix (.wig or .bw respectively).

Outputs:
	The program outputs a bedGraph flat file in stdout.

Parameters:
	// Unary operators
	wiggletools unit file
	wiggletools abs file
	wiggletools exp file
	wiggletools log file

	// Operators between a signal and a scalar
	wiggletools scale file factor
	wiggletools pow file exponent
	wiggletools exp file radix
	wiggletools log file base

	// Reduction operators
	wiggletools add file1 file2 ... 
	wiggletools mult file1 file2 ...
	wiggletools min file1 file2 ...
	wiggletools max file1 file2 ...
	wiggletools mean file1 file2 ...
	wiggletools var file1 file2 ...
	wiggletools stddev file1 file2 ...
	wiggletools median file1 file2 ...
		
	// Calculations
	wiggletools AUC file
	wiggletools pearson file1 file2

	// Other
	wiggletools --help
