WiggleTools 1.0
===============

Author: [Daniel Zerbino](mailto:zerbino@ebi.ac.uk)
Copyright holder: EMBL-EBI (Apache 2 License)

The WiggleTools package allows genomewide data files to be manipulated as numerical functions, equipped with all the standard functional analysis operators (sum, product, product by a scalar, comparators), and derived statistics (mean, median, variance, stddev, t-test, Wilcoxon's rank sum test, etc).

Installation
------------

WiggleTools requires two dependencies, the Kent and GSL (GNU scientific) libraries:

Installing the Kent library

First download the code:

```
git clone git://genome-source.cse.ucsc.edu/kent.git
setenv KENT_SRC $PWD/kent/src
cd $KENT_SRC
git checkout bigWigCat
```

Follow the instructions to compile the library, as well as the utilities *wigToBigWig* and *bigWigCat*.

You can checkout the master branch if you wish after creating these utilities. 

Obtaining WiggleTools
---------------------

If you didnt download WiggleTools yet:

```
git clone https://github.com/dzerbino/WiggleTools.git
```

Installing WiggleTools
----------------------

Once you installed the two previous libraries and set the KENT\_SRC environment variable, and downloaded the WiggleTools, you can compile the WiggleTools library:

```
cd wiggletools
make
```

The make process produces a number of outputs:

* A statically linked library in lib/
* A header for that library in inc/
* Various executables in bin/
* Documentation (this file) in doc/

There is not installation routine, meaning that you should copy the relevant files onto your path, library path, etc. Note that the executable does not require the libraries to be available.

Just to check, you can launch the tests:

```
make test
```

Basics
------

The WiggleTools library, and the derived program, are centered around the use of iterators. An iterator is a function which produces a sequence of values. The cool thing is that iterators can be built off other iterators, offering many combinations. 

The wiggletools executable is run by giving it a string which describes an iterator function, which it executes, printing the output into stdout.

```
wiggletools <program>
```

If you need a refresher:

```
wiggletools --help
```

Input files
-----------

By default, the executable recognizes the file format from the suffix of the file name:

* Wiggle files

```
wiggletools test/fixedStep.wig 
```

* BigWig files

```
wiggletools test/fixedStep.bw 
```

* BedGraph files

```
wiggletools test/bedfile.bg 
```

* Bed files

```
wiggletools test/overlapping.bed 
```

* BigBed files

```
wiggletools test/overlapping.bb 
```
Operators
---------

However, iterators can be constructed from other iterators, allowing arbitrarily complex constructs to be built. We call these iterators operators. In all the examples below, the iterators are built off simple file readers (for simplicity), but you are free to replace the inputs with other iterators.

1 Unary operators

The following operators are the most straightforward, because they only read data from a single other iterator.

* abs

Returns the absolute value of an iterators output:

```
wiggletools abs test/fixedStep.bw 
```

* log

Returns the natural log of an iterators output:

```
wiggletools log test/fixedStep.bw 
```

* scale

Returns an iterators output multiplied by a scalar (i.e. decimal number):

```
wiggletools scale 10 test/fixedStep.bw 
```

* gt

Returns 1 if the iterator is strictly greater than a given cutoff, 0 otherwise, and merges contiguous positions with the same output value into blocks:

```
wiggletools gt 5 test/fixedStep.bw 
```

This is useful to define regions in the *apply* function (see below).

* unit

Returns 1 if the operator is non-zero, 0 otherwise, and merges contiguous positions with the same output value into blocks:

```
wiggletools unit test/fixedStep.bw 
```

This is useful to define regions in the *apply* function (see below).

* isZero

Does not print anything, just exits with return value 1 (i.e. error) if it encounters a non-zero value:

```
wiggletools isZero test/fixedStep.bw 
```

* seek

Outputs only the points of an iterator within a given genomic region:

```
wiggletools seek chr1 2 8 test/fixedStep.bw 
```

2 Binary operators

The following operators read data from exactly two iterators, allowing comparisons:

* diff

Returns the difference between two iterators outputs:

```
wiggletools diff test/fixedStep.bw test/variableStep.bw 
```

* ratio

Returns the output of the first iterator divided by the output of the second (divisions by 0 are squashed, and no result is given for those bases):

```
wiggletools ratio test/fixedStep.bw test/variableStep.bw 
```

3 Multiplexed iterators

However, sometimes you want to compute statistics across many iterators. In this case, the function is followed by an arbitrary list of iterators, separated by spaces. The list is terminated by a colon (:) separated by spaces from other words. At the very end of a command string, the semi-colon can be omitted (see example in the example for *sum*)

* sum

The sum function sums all the listed iterators. The two following commands are equivalent:

```
wiggletools sum test/fixedStep.bw test/variable Step.bw :
wiggletools sum test/fixedStep.bw test/variable Step.bw
```

However, the semi-colon can be necessary for the program string to be unambiguous, e.g.:

```
wiggletools diff sum test/fixedStep.bw test/variable Step.bw \
            : test/fixedStep
```

* mult

Multiplies the subsequent list of iterators:

```
wiggletools mult test/fixedStep.bw test/variable Step.bw 
```

* mean

Computes the mean of the subsequent list of iterators at each position:

```
wiggletools mean test/fixedStep.bw test/variable Step.bw 
```

* median

Computes the median of the subsequent list of iterators at each position:

```
wiggletools median test/fixedStep.bw test/variable Step.bw 
```

* variance

Computes the variance of the subsequent list of iterators at each position:

```
wiggletools variance test/fixedStep.bw test/variable Step.bw 
```

* stddev

Computes the standard error of the subsequent list of iterators at each position:

```
wiggletools stddev test/fixedStep.bw test/variable Step.bw 
```

* min

Computes the minimum of the subsequent list of iterators at each position:

```
wiggletools min test/fixedStep.bw test/variable Step.bw 
```

* max

Computes the maximum of the subsequent list of iterators at each position:

```
wiggletools max test/fixedStep.bw test/variable Step.bw 
```

4 Comparing sets of sets

* Welch's t-test

Computes the two-tailed p-value of Welch's t-test comparing to sets of numbers, each assumed to have a normal distribution:

```
wiggletools ttest test/fixedStep.bw test/variableStep.bw test/fixedStep.wig \
            : test/fixedStep.wig test/variableStep.bw test/fixedStep.wig
```

* Wilcoxon's sum rank test

Non-parametric equivalent of the above:

```
wiggletools wilcoxon test/fixedStep.bw test/variableStep.bw test/fixedStep.wig \
            : test/fixedStep.wig test/variableStep.bw test/fixedStep.wig
```

5 Mapping a unary function to an iterator list:

If you wish to apply the same function to a list of iterators without typing redundant keywords, you can use the *map* function, which applies said operator to each element of the list:

```
wiggletools sum map ln test/fixedStep.bw test/variableStep.bw
```

Writing into files
------------------

Stdout is great and all, but sometimes you want to specify an output file on the command line without the use of pipes. This is done with the *write* function. It writes the output of an iterator into a wiggle file, and simultaneously returns the same output:

```
wiggletools write copy.wiggle test/fixedStep.wig 
```

The write instruction is itself an iterator, such that you can store data in a file, yet keep it in memory for more computation. For example the following computes the mean of two files, stores the result in a file, and also compares that result to a third file:

```
wiggletools diff test/fixedStep.bw \
write sum.wig mean test/fixedStep.bw test/variable Step.bw 
```

For convenience, if a command starts with a write instruction, the standard output is squashed. Otherwise, if you want to silence standard out, use the *do* command, which simply runs an iterator and returns nothing: 

```
wiggletools do test/fixedStep.wig 
```

If you wish to write into standard output, simply use the dash - symbol.

```
wiggletools write - test/fixedStep.wig 
```

If you wish to have your output in BedGraph format (takes more space but easier to parse line-by-line), use the write\_bg command:

```
wiggletools write_bg - test/fixedStep.wig 
```

Note that BedGraphs and the BedGraph sections within wiggle files are 0-based, whereas the `normal' wiggle lines have 1-based coordinates.

Statistics
----------

Sometimes, you just want a statistic across the genome. The following functions do not return a sequence of numbers, just a single number. All of these outputs are directed to an user defined output file (in this case results.txt) but you can put `-' for standard output:


* AUC

Computes the area under the curve (AUC) of the an iterator:

```
wiggletools AUC results.txt test/fixedStep.bw test/variable Step.bw 
```

* variance

Computes the variance of an iterator across all of its points:

```
wiggletools variance results.txt test/fixedStep.bw 
```

* pearson

Computes the Pearson correlation between two iterators across all their points:

```
wiggletools pearson results.txt test/fixedStep.bw test/fixedStep.bw 
```

* Apply

The apply function reads the regions from one iterator, then computes a given statistic on another iterator across those regions. It ignores regions with value 0.

```
wiggletools apply mean unit test/variableStep.bw test/fixedStep.bw
```

* Apply and Paste

This is a convenience wrapper around the above function: it reads the regions directly from a Bed file, then prints out each line of the file, with the resulting statistic appended at the end of the line. This is useful to keep identifiers and other metadata contained in the same file as the results:

```
wiggletools apply_paste output_file.txt mean test/overlapping.bed test/fixedStep.bw
```

Profiles
--------

To generate a fixed width summary of an iterator across a collection of regions, you can request the profiles function. This will print out the profiles, one for each region:

```
../bin/wiggletools profiles 3 test/overlapping.bed test/fixedStep.wig
```

If you just want a single profile, which sums up the results of all those profiles, you simply do:

```
../bin/wiggletools profile 3 test/overlapping.bed test/fixedStep.wig
```

Parallel processing
-------------------

To aid in running Wiggletools efficiently, a script, *parallelWiggletools.py* was designed to automate the batching of multiple jobs and the merging of their output. At the moment, this scripts requires an LSF job queueing system.

To run this script, you must provide first with a tab-delimited file that specifies the names and legnths of all the chromosomes in your genome, see test/chrom\_sizes for an example.

You then specify a Wiggletools command, note how the write function now points to a BigWig file:

```
parallelWiggletools.py test/chrom_sizes 'write copy.bw test/fixedStep.bw'
```

Because these are asynchronous jobs, they generate a bunch of files as input, stdout and stderr. If these files are annoying to you, you can change the DUMP\_DIR variable in the parallelWiggleTools script, to another directory which is visible to all the nodes in the LSF farm.
