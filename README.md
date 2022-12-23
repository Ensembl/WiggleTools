[![GitHub license](https://img.shields.io/github/license/Ensembl/WiggleTools)](https://github.com/Ensembl/WiggleTools/blob/master/LICENSE)
[![GitHub stars](https://img.shields.io/github/stars/Ensembl/WiggleTools)](https://github.com/Ensembl/WiggleTools/stargazers)
[![GitHub forks](https://img.shields.io/github/forks/Ensembl/WiggleTools)](https://github.com/Ensembl/WiggleTools/network)
[![GitHub issues](https://img.shields.io/github/issues/Ensembl/WiggleTools)](https://github.com/Ensembl/WiggleTools/issues)

# WiggleTools 1.2

Author: [Daniel Zerbino](mailto:zerbino@ebi.ac.uk)

Copyright holder: [EMBL-European Bioinformatics Institute](https://www.ebi.ac.uk/) (Apache 2 License)

The WiggleTools package allows genomewide data files to be manipulated as numerical functions, equipped with all the standard functional analysis operators (sum, product, product by a scalar, comparators), and derived statistics (mean, median, variance, stddev, t-test, Wilcoxon's rank sum test, etc).

## Conda Installation

Install [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/), then run:

```
conda install -c bioconda wiggletools
```

## Brew Installation

Install [Homebrew](https://brew.sh/), then run:

```
brew install brewsci/bio/wiggletools
```

## Docker Installation

Pull the latest image from Dockerhub:
```
docker pull ensemblorg/wiggletools:latest
```

Run the resulting wiggletools executable, bind-mounting the current working directory into the container:
```
docker container run --rm --mount type=bind,source="$(pwd)",target=/mnt ensemblorg/wiggletools  [...arguments...]
```

## Guix Installation

Install [GNU Guix](https://guix.gnu.org), then run:

```
guix pull
guix install wiggletools
```

## Build from source

### Pre-requisites

WiggleTools requires three main dependencies: [LibBigWig](https://github.com/dpryan79/libBigWig), [HTSLib](https://github.com/samtools/htslib) and [GSL (GNU scientific)](https://www.gnu.org/software/gsl/) libraries. They themselves require [zlib](http://www.zlib.net/) [bzip2](https://sourceware.org/bzip2/) and [libcurl](https://curl.haxx.se/download.html).

#### Installing LibBigWig

```
git clone https://github.com/dpryan79/libBigWig.git
cd libBigWig
make install
```

#### Installing the htslib library

```
git clone --recurse-submodules https://github.com/samtools/htslib.git
cd htslib 
make install
```

#### Installing the GSL library
```
wget ftp://www.mirrorservice.org/sites/ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz 
tar -xvzpf gsl-latest.tar.gz
cd gsl*
./configure
make
make install
```

### Installing WiggleTools

If you didn't download WiggleTools yet:

```
git clone https://github.com/Ensembl/WiggleTools.git
```

Once you installed the previous libraries and downloaded WiggleTools, you can compile the WiggleTools library:

```
cd WiggleTools
make
```

The make process produces a number of outputs:

* A statically linked library in lib/
* A header for that library in inc/
* Various executables in bin/

There is no installation routine, meaning that you should copy the relevant files onto your path, library path, etc. Note that the executable does not require the libraries to be available.

If the system cannot find 'gsl/gsl_cdf.h' then you need to install the [GNU scientific library](http://www.gnu.org/software/gsl/)

Just to check, you can launch the tests (requires Python):

```
make test
```

## Basics

The WiggleTools library, and the derived program, are centered around the use of iterators. An iterator is a function which produces a sequence of values. The cool thing is that iterators can be built off other iterators, offering many combinations. 

The wiggletools executable is run by giving it a string which describes an iterator function, which it executes, printing the output into stdout.

```
wiggletools <program>
```

If you need a refresher:

```
wiggletools --help
```

If you are an intensive user, you may find that processing many files may break limits on commandline commands, especially if shelling out from a scripting language. You may copy the program into a text file, then execute it:

```
wiggletools run program.txt
```

## Input files

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

* Cram files

Requires a .bai index file in the same directory

```
wiggletools test/cram.cram
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

## Streaming data

You can stream data into WiggleTools, e.g.:

```
cat test/fixedStep.wig | wiggletools -
```

The input data is assumed to be in Wig or BedGraph format, but can also be in Sam format:

```
samtools view test/bam.bam | wiggletools sam -
```


## Operators


However, iterators can be constructed from other iterators, allowing arbitrarily complex constructs to be built. We call these iterators operators. In all the examples below, the iterators are built off simple file readers (for simplicity), but you are free to replace the inputs with other iterators.

### 1 Unary operators

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

Returns contiguous boolean regions where the iterator is strictly greater than a given cutoff:

```
wiggletools gt 5 test/fixedStep.bw 
```

This is useful to define regions in the *apply* function, or to compute information content (see below).

* lt

Returns contiguous boolean regions where the iterator is strictly less than a given cutoff:

```
wiggletools lt 5 test/fixedStep.bw 
```

This is useful to define regions in the *apply* function, or to compute information content (see below).

* gte

Returns contiguous boolean regions where the iterator is greater than or equal to a given cutoff:

```
wiggletools gte 5 test/fixedStep.bw 
```

This is useful to define regions in the *apply* function, or to compute information content (see below).

* lte

Returns contiguous boolean regions where the iterator is less than or equal to a given cutoff:

```
wiggletools lte 5 test/fixedStep.bw 
```

This is useful to define regions in the *apply* function, or to compute information content (see below).

* unit

Returns 1 if the operator is non-zero, 0 otherwise, and merges contiguous positions with the same output value into blocks:

```
wiggletools unit test/fixedStep.bw 
```

This is useful to define regions in the *apply* function (see below).

* coverage

Returns a coverage plot of overlapping regions, typically read from a bed file:  

```
wiggletools coverage test/overlapping.bed
```

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

* bin

Sums results into fixed-size bins

```
wiggletools bin 2 test/fixedStep.bw 
```

* toInt

Casts the iterator's output to an `int`, effectively rounding any floating point values toward zero.

```
wiggletools toInt test/fixedStep.bw
```

* floor

Returns the floor of a iterator's output. Note that floor rounds the output toward negative infinity.

```
wiggletools floor test/fixedStep.bw
```

* shiftPos

Returns the iterator given with start and end positions shifted downwards by a specified value. Note the given value must be non-negative, as default behavior is to shift coordinates toward zero.

```
wiggletools shiftPos 10 test/fixedStep.bw
```

### 2 Binary operators

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

* trim

Same as above but trims the regions to the overlapping portions:

```
wiggletools trim test/fixedStep.bw test/variableStep.bw
```

* trimFill

Same as trim, but fills in trimmed regions with the default value of the second iterator.

```
wiggletools trimFill test/fixedStep.bw test/overlapping_coverage.wig
```

* nearest

Returns the regions of the second iterator and their distance to the nearest region in the first iterator.

```
wiggletools nearest test/fixedStep.bw test/variableStep.bw 
```

### 3 Multiplexed iterators

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

* var

Computes the variance of the subsequent list of iterators at each position:

```
wiggletools var test/fixedStep.bw test/variableStep.bw 
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

### 4 Comparing sets of sets

* Welch's t-test

Computes the two-tailed p-value of Welch's t-test comparing two sets of numbers, each assumed to have a normal distribution:

```
wiggletools ttest test/fixedStep.bw test/variableStep.bw test/fixedStep.wig \
            : test/fixedStep.wig test/variableStep.bw test/fixedStep.wig
```

* F-test

Computes the p-value of the F-test comparing sets of numbers, each assumed to have a normal distribution:

```
wiggletools ftest test/fixedStep.bw test/variableStep.bw test/fixedStep.wig \
            : test/fixedStep.wig test/variableStep.bw test/fixedStep.wig
```

* Wilcoxon's sum rank test

Non-parametric equivalent of the above:

```
wiggletools wilcoxon test/fixedStep.bw test/variableStep.bw test/fixedStep.wig \
            : test/fixedStep.wig test/variableStep.bw test/fixedStep.wig
```

### 5 Mapping a unary function to an iterator list:

If you wish to apply the same function to a list of iterators without typing redundant keywords, you can use the *map* function, which applies said operator to each element of the list:

```
wiggletools sum map ln test/fixedStep.bw test/variableStep.bw
wiggletools sum scale -1 test/fixedStep.bw test/variableStep.bw
```

## Writing into files

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

The BedGraph output respects the layout of the input iterator, however if you want consecutive regions with the same value to be collapsed together, you can use the "compress" keyword:

```
wiggletools write_bg - compress test/fixedStep.wig 
```

## Writing multidimensional wiggles into files

Sometimes, for your own reasons, you may want to print out multiple wiggles side by side. This can be done with the designated *mwrite* and *mwrite_bg* operators

**Warning!!** Multidimensional wiggles are not part of the BigWig/BedGraph specs, and will probably spark an error with the Kent apps. These are solely designed for your own usage.

```
wiggletools mwrite_bg - test/overlapping.bed test/fixedStep.bw
```

## Statistics

Sometimes, you just want a statistic across the genome. The following functions do not return a sequence of numbers, just a single number. Some of these integrating functions have "I" appended to them to distinguish them from the iterators with related (yet different) functions. 


* AUC

Computes the area under the curve (AUC) of an iterator:

```
wiggletools AUC test/fixedStep.bw
```

* meanI

Computes the mean of an iterator across all of its points:

```
wiggletools meanI test/fixedStep.bw 
```

* varI

Computes the variance of an iterator across all of its points:

```
wiggletools varI test/fixedStep.bw 
```

* stddevI

Computes the standard deviation of an iterator across all of its points:

```
wiggletools stddevI test/fixedStep.bw 
```

* CVI 

Computes the coefficient of variation of an iterator across all of its points:

```
wiggletools CVI test/fixedStep.bw 
```

* maxI

Computes the maximum of an iterator across all of its points:

```
wiggletools maxI test/fixedStep.bw 
```

* minI

Computes the minimum of an iterator across all of its points:

```
wiggletools minI test/fixedStep.bw 
```

* pearson

Computes the Pearson correlation between two iterators across all their points:

```
wiggletools pearson test/fixedStep.bw test/fixedStep.bw
```

* energy

Computes the energy density at a given wavelength:

```
wiggletools energy 10 test/fixedStep.bw
```

## Chaining statistics

All the above functions are actually iterators that transmit the same data as they are given, e.g.:

```
wiggletools test/fixedStep.bw 
wiggletools scale 1 AUC test/fixedStep.bw 
```

This allows you to plug multiple statistics in a dandelion chain off the same iterator. Note how results of the operators are concatenated as they are read from left to right:

```
wiggletools meanI varI minI maxI test/fixedStep.bw 
```

If you want to save the output of a statistic into a file, you can use the print statement:

```
wiggletools print output.txt AUC test/fixedStep.bw
```

As with other write functions, if a command starts with a print statement, the standard output is squashed.

Apply
-----

The *apply* function reads the regions from one iterator, then computes a given statistic on another iterator across those regions. You can chain the statistics as above. Because of this feature, the *apply* operator returns a multiplexer (i.e. a multidimensional wiggle), hence the *mwrite* operator before it: 

```
wiggletools mwrite_bg - apply meanI stddevI unit test/variableStep.bw test/fixedStep.bw
```

For convenience, if the *apply* operator is used in a context which expects a standard unidimensional wiggle, it is transformed into one. In particular, if you computed multiple statistics in parallel as above, only the first is retained. 

* Apply and Paste

This is a convenience wrapper around the above function: it reads the regions directly from a Bed file, then prints out each line of the file, with the resulting statistic appended at the end of the line. This is useful to keep identifiers and other metadata contained in the same file as the results. Note that the *mwrite* operator is unnecessary in this case:

```
wiggletools apply_paste output_file.txt meanI test/overlapping.bed test/fixedStep.bw
```

## Profiles

To generate a fixed width summary of an iterator across a collection of regions, you can request the profiles function. This will print out the profiles, one for each region:

```
wiggletools profiles results.txt 3 test/overlapping.bed test/fixedStep.wig
```

If you just want a single profile, which sums up the results of all those profiles, you simply do:

```
wiggletools profile results.txt 3 test/overlapping.bed test/fixedStep.wig
```

As above, the output file name can be replaced by a dash (-) to print to standard output.

## Histograms

To generate a histogram of values across the iterator, simply use the *histogram* command. The number of bins must be pre-defined:

```
wiggletools histogram results.txt 10 test/fixedStep.wig
```

A histogram can hold multiple distributions:

```
wiggletools histogram results.txt 10 test/fixedStep.wig test/variableStep.wig
```

The format of the output is hopefully rather self explanatory: each line starts with the midpoint value of a bin, and the values for that bin, tabbed-delimited.

The algorithm used to compute these histograms is approximate: it adapts the width of the bins to the data received, and requires very little memory or computation. However, the values of the bins is not quite exact, as some points might be counted in a neighbouring bin to the one they should belong to. Normally, over a large datasets, these approximations should roughly even out. 


## Parallel processing

To aid in running Wiggletools efficiently, a script, *parallelWiggletools.py* was designed to automate the batching of multiple jobs and the merging of their output. At the moment, this scripts requires an LSF job queueing system.

To run this script, you must provide first with a tab-delimited file that specifies the names and legnths of all the chromosomes in your genome, see test/chrom\_sizes for an example.

You then specify a Wiggletools command, note how the write function now points to a BigWig file:

```
parallelWiggletools.py test/chrom_sizes 'write copy.bw test/fixedStep.bw'
```

Because these are asynchronous jobs, they generate a bunch of files as input, stdout and stderr. If these files are annoying to you, you can change the DUMP\_DIR variable in the parallelWiggleTools script, to another directory which is visible to all the nodes in the LSF farm.

## Default Values

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

## Creating your own functions

If you're looking for guidance on creating your own iterators, you can have a look at:
* Iterator -> Float: src/statistics.c, MeanIntegrator
* Iterator -> Iterator: src/unaryOps.c, ScaleWiggleIterator
* Set of iterators -> Iterator: src/reducers.c, MaxReduction
* Set of Set of Iterators -> Iterator: src/setComparisons.c, TTestReduction

More info on the basic objects at:
* wiggleIterator.h: Simple iterator.
* multiplexer.h: Set of synchronised iterators.
* multiSet.h: Set of set of synchronised iterators.

You may need some help hooking your new functions to the parser, we can help you out.

## Citing WiggleTools

[Zerbino DR, Johnson N, Juettemann T, Wilder SP and Flicek PR: **WiggleTools: parallel processing of large collections of genome-wide datasets for visualization and statistical analysis.** *Bioinformatics* 2014 **30**:1008-1009.](http://bioinformatics.oxfordjournals.org/content/30/7/1008)
