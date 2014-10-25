import sys
import os
import shutil
import subprocess

def test(cmd):
	print 'Testing: %s' % cmd
	return subprocess.call(cmd, shell = True)

def testOutput(cmd):
	print 'Testing: %s' % cmd
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell = True)
	assert p.wait() == 0
	out, err = p.communicate()
	return out

if not os.path.exists('tmp'):
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

# Testing Bed and BigBed
assert test('../bin/wiggletools do isZero diff overlapping.bed overlapping.bb') == 0

# Testing Wig and BigWig
assert test('../bin/wiggletools do isZero diff variableStep.bw variableStep.wig') == 0

# Testing VCF and BCF
assert test('../bin/wiggletools do isZero diff vcf.vcf bcf.bcf') == 0

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
assert test('../bin/wiggletools do isZero diff smooth 1 fixedStep.wig fixedStep.wig') == 0

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

# Test min
assert float(testOutput('../bin/wiggletools print - minI fixedStep.wig')) == 0

# Test max
assert float(testOutput('../bin/wiggletools print - maxI fixedStep.wig')) == 9

assert test('diff tmp expected') == 0
shutil.rmtree('tmp')

print 'All tests OK'
