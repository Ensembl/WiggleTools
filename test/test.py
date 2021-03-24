import sys
import os
import shutil
import subprocess

def test(cmd):
	print('Testing: %s' % cmd)
	return subprocess.call(cmd, shell = True)

def testOutput(cmd):
	print('Testing: %s' % cmd)
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell = True)
	assert p.wait() == 0
	out, err = p.communicate()
	return out

if os.path.exists('tmp'):
	shutil.rmtree('tmp')
os.mkdir('tmp')

# Negative control
assert test('../bin/wiggletools do isZero diff fixedStep.bw variableStep.wig') == 1

# Test trailing token check:
assert test('../bin/wiggletools diff fixedStep.bw variableStep.wig fixedStep.wig') == 1

# Positive control
assert test('../bin/wiggletools do isZero diff fixedStep.bw fixedStep.wig') == 0

# Testing ratios and offset
assert test('../bin/wiggletools do isZero offset -1 ratio variableStep.bw variableStep.wig') == 0

# Testing BAM & BedGraph
assert test('../bin/wiggletools do isZero diff bam.bam pileup.bg') == 0

# Testing BAM & CRAM
assert test('../bin/wiggletools do isZero diff bam.bam cram.cram') == 0

# Testing BAM & SAM
assert test('../bin/wiggletools do isZero diff bam.bam sam.sam') == 0

# Testing fast BAM and SAM
assert test('../bin/wiggletools do isZero diff read_count bam.bam read_count sam.sam') == 0

# Testing BAM & SAM
assert test('cat sam.sam | ../bin/wiggletools do isZero diff bam.bam sam -') == 0

# Testing Bed and BigBed
assert test('../bin/wiggletools do isZero diff overlapping.bed overlapping.bb') == 0

# Testing Wig and BigWig
assert test('../bin/wiggletools do isZero diff variableStep.bw variableStep.wig') == 0

# Testing VCF and BCF
assert test('../bin/wiggletools do isZero diff vcf.vcf bcf.bcf') == 0

# Testing BAM & BedGraph
assert test('../bin/wiggletools do isZero seek GL000200.1 1 1000 diff bam.bam pileup.bg') == 0

# Testing BAM & CRAM
assert test('../bin/wiggletools do isZero seek GL000200.1 1 1000 diff bam.bam cram.cram') == 0

# Testing SAM & BedGraph
assert test('../bin/wiggletools do isZero seek GL000200.1 1 1000 diff sam.sam pileup.bg') == 0

# Testing Bed and BigBed
assert test('../bin/wiggletools do isZero seek chr1 2 6 diff overlapping.bed overlapping.bb') == 0

# Testing Wig and BigWig
assert test('../bin/wiggletools do isZero seek chr1 2 6 diff variableStep.bw variableStep.wig') == 0

# Testing VCF and BCF
assert test('../bin/wiggletools do isZero seek chr1 2 6 diff vcf.vcf bcf.bcf') == 0

# Testing sum, scale and multiplexers
assert test('../bin/wiggletools do isZero diff sum fixedStep.bw fixedStep.bw : scale 2 fixedStep.bw') == 0

# Testing open-ended lists
assert test('../bin/wiggletools do isZero diff sum fixedStep.bw fixedStep.bw : sum fixedStep.bw fixedStep.bw ') == 0

# Testing map
assert test('../bin/wiggletools do isZero diff ln fixedStep.bw sum map ln fixedStep.bw ') == 0

# Testing log and exponential
assert test('../bin/wiggletools do isZero diff ln exp fixedStep.bw fixedStep.wig') == 0

# Testing power and multiplication
assert test('../bin/wiggletools do isZero diff pow 2 fixedStep.bw mult fixedStep.wig fixedStep.wig') == 0

# Testing smoothing
# TODO : Find better test
# assert test('../bin/wiggletools do isZero diff smooth 2 fixedStep.wig fixedStep.wig') == 0

# Testing filters
assert float(testOutput('../bin/wiggletools AUC lt 4 fixedStep.wig')) == 4
assert float(testOutput('../bin/wiggletools AUC lte 4 fixedStep.wig')) == 5
assert float(testOutput('../bin/wiggletools AUC gte 4 fixedStep.wig')) == 6
assert float(testOutput('../bin/wiggletools AUC gt 4 fixedStep.wig')) == 5

# Testing apply
assert test('../bin/wiggletools apply_paste tmp/regional_means.txt meanI overlapping.bed fixedStep.wig') == 0

# Testing pearson
assert test('../bin/wiggletools print tmp/pearson.txt pearson fixedStep.wig variableStep.wig') == 0

# Testing profiles
assert test('../bin/wiggletools profiles tmp/profiles.txt 3 overlapping.bed fixedStep.wig') == 0

# Testing profile
assert test('../bin/wiggletools profile tmp/profile.txt 3 overlapping.bed fixedStep.wig') == 0

# Test overlap
assert test('../bin/wiggletools do isZero diff fixedStep.wig overlaps fixedStep.wig fixedStep.wig') == 0

# Test nearest #1
assert test('../bin/wiggletools write_bg tmp/nearest_overlapping.bg nearest variableStep.wig overlapping.bed') == 0

# Test nearest #1
assert test('../bin/wiggletools write_bg tmp/nearest_fixedStep.bg nearest variableStep.wig fixedStep.bw') == 0

# Test min
assert float(testOutput('../bin/wiggletools print - minI fixedStep.wig')) == 0

# Test max
assert float(testOutput('../bin/wiggletools print - maxI fixedStep.wig')) == 9

# Test coverage
assert test('../bin/wiggletools do isZero diff overlapping_coverage.wig coverage overlapping.bed') == 0

#Test trim
assert test('../bin/wiggletools do isZero diff trim overlapping.bed variableStep.wig mult overlapping.bed variableStep.wig') == 0

#Test floor
assert test('../bin/wiggletools do isZero diff floor fixedStep.wig floor fixedStep.wig') == 0

#Test toInt
assert test('../bin/wiggletools do isZero diff toInt fixedStep.wig toInt variableStep.wig') == 1

# Test program file
assert test('../bin/wiggletools run program.txt') == 0

assert test('diff tmp expected') == 0


shutil.rmtree('tmp')

print('All tests OK')
