WiggleTools

Copyright Daniel Zerbino, 2013.
zerbino@ebi.ac.uk

This library parses wiggle files and executes various operations on them streaming through lazy evaluators. Currently more a skeleton of a library, just type make for compilation.

Installation:
	- Install Jim Kent's source code
	- Define the environment variable KENT_SRC to point to the /path/to/kent/src/ directory
	- In this directory type 'make'
	- The binary file wiggletools should be produced in the same directory. You can add it to your PATH, or move it onto your PATH

Running:

Inputs:
	The program takes in Wig, BedGraph and BigWig files, which are distinguished thanks to their suffix (.wig, .bg or .bw respectively). A dash '-' signifies inputing a flat file through stdin.

Outputs:
	The program outputs a bedGraph flat file in stdout.

Parameters:
	// Unary operators
	wiggletools unit file
	wiggletools exp file
	wiggletools log file

	// Operators between a signal a scalar
	wiggletools scale file factor
	wiggletools pow file exponent
	wiggletools exp file radix
	wiggletools log file base

	// Binary operators between two signal files
	wiggletools add file1 file2
	wiggletools mult file1 file2

