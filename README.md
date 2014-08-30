WiggleTools 1.0
===============

Author: [Daniel Zerbino](mailto:zerbino@ebi.ac.uk)

Copyright holder: EMBL-EBI (Apache 2 License)

The WiggleTools package allows genomewide data files to be manipulated as numerical functions, equipped with all the standard functional analysis operators (sum, product, product by a scalar, comparators), and derived statistics (mean, median, variance, stddev, t-test, Wilcoxon's rank sum test, etc).

Easy Installation (Experimental)
--------------------------------

WiggleTools requires quite a few dependencies, and this can be a hassle to install.

To speed up the process you can use the easy installation script:

```
easy_install.sh
```
It will test for the presence of pre-existing installations, prompt you before each download then install as appropriate. If you already installed Tabix or the Kent source code, be sure to set the $TABIX_SRC and $KENT_SRC environment variables to avoid a redundant installation. Because C libraries are being installed, root permissions are required. 

If you are on Mac or Ubuntu, the script will use HomeBrew or apt-get. If you have a different package installer, it should be pretty easy to install hooks for it in the script. 

Note that this script is quite experimental, and your system is different from mine, so please be nice to it and let me know if it can be corrected in any way.

Installation
------------

WiggleTools requires three dependencies, the Kent and GSL (GNU scientific) libraries:

Installing the Kent library

First download the code:

```
git archive --format=zip -9 --remote=git://genome-source.cse.ucsc.edu/kent.git beta src/userApps > userApps.zip
unzip -d userApps -j userApps.zip
rm userApps.zip
cd userApps
make fetchSource
make
setenv KENT_SRC $PWD/kent/src
# or, if you use bash...
export KENT_SRC=$PWD/kent/src
```
Ensure that your path points to the userApps/bin directory.

Installing the Tabix library

```
git clone https://github.com/samtools/tabix.git
setenv TABIX_SRC $PWD/tabix
# or, if you use bash...
export TABIX_SRC=$PWD/tabix
cd tabix
make
```

Obtaining WiggleTools
---------------------

If you didnt download WiggleTools yet:

```
git clone https://github.com/Ensembl/WiggleTools.git
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

There is not installation routine, meaning that you should copy the relevant files onto your path, library path, etc. Note that the executable does not require the libraries to be available.

If the system complains that it cannot find -lssl or -lcrypto then you need to install the [libssl runtime and development packages](http://www.openssl.org/)

If the system cannot find 'gsl/gsl_cdf.h' then you need to install the [GNU scientific library](http://www.gnu.org/software/gsl/)

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

* Bam files

Requires a .bai index file in the same directory

```
wiggletools test/bam.bam
```

* VCF files

```
wiggletools test/vcf.vcf
```

* BCF files

Requires a .tbi index file in the same directory

```
wiggletools test/bcf.bcf
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

* ln

Returns the natural log of an iterators output:

```
wiggletools ln test/fixedStep.bw 
```

* log

Returns the logarithm in an arbitrary base of an iterators output:

```
wiggletools log 10 test/fixedStep.bw 
```

* scale

Returns an iterator's output multiplied by a scalar (i.e. decimal number):

```
wiggletools scale 10 test/fixedStep.bw 
```

* offset

Returns an iterator's output added to a scalar (i.e. decimal number):

```
wiggletools offset 10 test/fixedStep.bw 
```


* gt

Returns 1 if the iterator is strictly greater than a given cutoff, 0 otherwise, and merges contiguous positions with the same output value into blocks:

```
wiggletools gt 5 test/fixedStep.bw 
```

This is useful to define regions in the *apply* function, or to compute information content (see below).

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

* overlaps

Returns the output of the second iterator that overlaps regions of the first.

```
wiggletools overlaps test/fixedStep.bw test/variableStep.bw 
```

3 Multiplexed iterators

However, sometimes you want to compute statistics across many iterators. In this case, the function is followed by an arbitrary list of iterators, separated by spaces. The list is terminated by a colon (:) separated by spaces from other words. At the very end of a command string, the colon can be omitted (see example in the example for *sum*)

* sum

The sum function sums all the listed iterators. The two following commands are equivalent:

```
wiggletools sum test/fixedStep.bw test/variableStep.bw :
wiggletools sum test/fixedStep.bw test/variableStep.bw
```

However, the colon can be necessary for the program string to be unambiguous, e.g.:

```
wiggletools diff sum test/fixedStep.bw test/variableStep.bw \
            : test/fixedStep
```

* mult

Multiplies the subsequent list of iterators:

```
wiggletools mult test/fixedStep.bw test/variableStep.bw 
```

* mean

Computes the mean of the subsequent list of iterators at each position:

```
wiggletools mean test/fixedStep.bw test/variableStep.bw 
```

* median

Computes the median of the subsequent list of iterators at each position:

```
wiggletools median test/fixedStep.bw test/variableStep.bw 
```

* variance

Computes the variance of the subsequent list of iterators at each position:

```
wiggletools variance test/fixedStep.bw test/variableStep.bw 
```

* stddev

Computes the standard error of the subsequent list of iterators at each position:

```
wiggletools stddev test/fixedStep.bw test/variableStep.bw 
```

* entropy

Computes the Shannon entropy of the subsequent list of iterators at each position, separating 0 from non-0 values. This is probably most useful with the gt (greater than) filter:

```
wiggletools entropy gt 5 test/fixedStep.bw test/overlapping.bb
```


* CV

Computes the coefficient of variation ( = standard deviation / mean) of the subsequent list of iterators at each position:

```
wiggletools CV test/fixedStep.bw test/variableStep.bw 
```

* min

Computes the minimum of the subsequent list of iterators at each position:

```
wiggletools min test/fixedStep.bw test/variableStep.bw 
```

* max

Computes the maximum of the subsequent list of iterators at each position:

```
wiggletools max test/fixedStep.bw test/variableStep.bw 
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
wiggletools sum scale -1 test/fixedStep.bw test/variableStep.bw
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
write sum.wig mean test/fixedStep.bw test/variableStep.bw 
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

Sometimes, you just want a statistic across the genome. The following functions do not return a sequence of numbers, just a single number. Some of these integrating functions have "I" appended to them to distinguish them from the iterators with related (yet different) functions. Note the print statement at the beginning of the command.


* AUC

Computes the area under the curve (AUC) of the an iterator:

```
wiggletools print - AUC test/fixedStep.bw test/variableStep.bw 
```

* meanI

Computes the mean of an iterator across all of its points:

```
wiggletools print - meanI test/fixedStep.bw 
```

* varI

Computes the variance of an iterator across all of its points:

```
wiggletools print - varI test/fixedStep.bw 
```

* stddevI

Computes the standard deviation of an iterator across all of its points:

```
wiggletools print - stddevI test/fixedStep.bw 
```

* CVI 

Computes the coefficient of variation of an iterator across all of its points:

```
wiggletools print - CVI test/fixedStep.bw 
```

* maxI

Computes the maximum of an iterator across all of its points:

```
wiggletools print - maxI test/fixedStep.bw 
```

* minI

Computes the minimum of an iterator across all of its points:

```
wiggletools print - minI test/fixedStep.bw 
```

* pearson

Computes the Pearson correlation between two iterators across all their points:

```
wiggletools print - pearson test/fixedStep.bw test/fixedStep.bw 
```

* Chaining statistics:

All the above functions are actually iterators that transmit the same data as they are given, e.g.:

```
wiggletools test/fixedStep.bw 
wiggletools AUC test/fixedStep.bw 
```

This allows you to plug multiple statistics in a dandelion chain off the same iterator. Note how the print statement simply concatenates the results of the operators as they are reads from left to right:

```
wiggletools print - meanI varI minI maxI test/fixedStep.bw 
```

* Apply

The apply function reads the regions from one iterator, then computes a given statistic on another iterator across those regions. It ignores regions with value 0.

```
wiggletools apply meanI unit test/variableStep.bw test/fixedStep.bw
```

* Apply and Paste

This is a convenience wrapper around the above function: it reads the regions directly from a Bed file, then prints out each line of the file, with the resulting statistic appended at the end of the line. This is useful to keep identifiers and other metadata contained in the same file as the results:

```
wiggletools apply_paste output_file.txt meanI test/overlapping.bed test/fixedStep.bw
```

Profiles
--------

To generate a fixed width summary of an iterator across a collection of regions, you can request the profiles function. This will print out the profiles, one for each region:

```
wiggletools profiles results.txt 3 test/overlapping.bed test/fixedStep.wig
```

If you just want a single profile, which sums up the results of all those profiles, you simply do:

```
wiggletools profile results.txt 3 test/overlapping.bed test/fixedStep.wig
```

As above, the output file name can be replaced by a dash (-) to print to standard output.

Histograms
----------

To generate a histogram of values across the iterator, simply use the *histogram* command. The number of bins must be pre-defined:

```
wiggletools histogram results.txt 10 test/fixedStep.wig
```

The format of the output is hopefully rather self explanatory: each line starts with the lower bound of a bin, and the value for that bin. The last line contains the upper bound of the last bin. 

The algorithm used to compute these histograms is approximate: it adapts the width of the bins to the data received, and requires very little memory or computation. However, the values of the bins is not quite exact, as some points might be counted in a neighbouring bin to the one they should belong to. Normally, over a large datasets, these approximations should roughly even out. 


Parallel processing
-------------------

To aid in running Wiggletools efficiently, a script, *parallelWiggletools.py* was designed to automate the batching of multiple jobs and the merging of their output. At the moment, this scripts requires an LSF job queueing system.

To run this script, you must provide first with a tab-delimited file that specifies the names and legnths of all the chromosomes in your genome, see test/chrom\_sizes for an example.

You then specify a Wiggletools command, note how the write function now points to a BigWig file:

```
parallelWiggletools.py test/chrom_sizes 'write copy.bw test/fixedStep.bw'
```

Because these are asynchronous jobs, they generate a bunch of files as input, stdout and stderr. If these files are annoying to you, you can change the DUMP\_DIR variable in the parallelWiggleTools script, to another directory which is visible to all the nodes in the LSF farm.

Default Values
--------------

A basic underlying question is how to deal with missing values. In some cases, no value in a BigWig file implicitly means 0, typically when working with coverage statistics or peaks. However, sometimes, you want positions with no values to be disregarded. 

* To deal with this, every iterator has a default value. By default, any file being read has a default value of 0. The default value of composed iterators is computed from their inputs. For example, if B is equal to A multiplied by 10, then the default value of B is 10 times that of A. The default value of an iterator can be directly set with the *default* keyword:

```
wiggletools sum test/fixedStep.wig test/variableStep.wig
wiggletools sum test/fixedStep.wig default 10 test/variableStep.wig
```

* When a set of iterators A1, A2 ... is composed by a n-ary iterator M, M will skip the regions which are skipped by all the input iterators. However, in the presence of more than one input iterators that do not perfectly overlap, there will be regions which are covered by say A1, but not A2. Two behaviours are defined: if M is *strict*, it skips those regions, else it replaces missing values with the corresponding default values. 

By default, n-ary iterators are not strict, but they can be made so with the *strict* keyword after the n-ary function name:

```
wiggletools sum test/fixedStep.wig test/variableStep.wig
wiggletools sum strict test/fixedStep.wig test/variableStep.wig
```

* However, when integrating statistics across the genome, or regions of the genome, missing values are discarded by default, because it is not known which regions need to be filled in. 

The *fillIn* operator allows you to define which regions should be filled in with the default value. It takes in two iterators, and behaves like a default multiplexer. Wherever possible, it takes the value of the second iterator, and otherwise takes the default_value of the second iterator. Note that this incurs a slight processing cost, so only use this operator at the last step of your computations, right before computing a statistic. 

```
wiggletools test/variableStep.wig
wiggletools meanI - test/variableStep.wig
wiggletools fillIn test/fixedStep.wig test/variableStep.wig
wiggletools meanI - fillIn test/fixedStep.wig test/variableStep.wig
```

For convenience, the fillIn keyword can be used in the apply commands, as is:

```
wiggletools apply meanI unit test/variableStep.bw test/fixedStep.bw
wiggletools apply meanI fillIn unit test/variableStep.bw test/fixedStep.bw
```
