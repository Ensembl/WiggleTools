import sys
import os
import shutil
import subprocess

def test(cmd):
	print 'Testing: %s' % cmd
	return subprocess.call(cmd, shell = True)

if not os.path.exists('tmp'):
	os.mkdir('tmp')

# Negative control
assert test('../bin/wiggletools isZero diff fixedStep.bw variableStep.wig') == 1

# Positive control
assert test('../bin/wiggletools isZero diff fixedStep.bw fixedStep.wig') == 0

# Testing ratios
assert test('../bin/wiggletools isZero log ratio fixedStep.bw fixedStep.wig') == 0

# Testing BAM & BedGraph 
assert test('../bin/wiggletools isZero diff bam.bam pileup.bg') == 0

# Testing Bed and BigBed
assert test('../bin/wiggletools isZero diff overlapping.bed overlapping.bb') == 0

# Testing Wig and BigWig
assert test('../bin/wiggletools isZero diff variableStep.bw variableStep.wig') == 0

# Testing sum and scale
assert test('../bin/wiggletools isZero diff sum fixedStep.bw fixedStep.bw : scale 2 fixedStep.bw fixedStep.bw ') == 0

# Testing log and exponential
assert test('../bin/wiggletools isZero diff ln exp fixedStep.bw fixedStep.wig') == 0

# Testing power and multiplication
assert test('../bin/wiggletools isZero diff pow 2 fixedStep.bw mult fixedStep.wig fixedStep.wig') == 0

# Testing smoothing
assert test('../bin/wiggletools isZero diff smooth 1 fixedStep.wig fixedStep.wig') == 0

# Testing apply
assert test('../bin/wiggletools apply tmp/regional_means.txt mean overlapping.bed fixedStep.wig') == 0

# Testing pearson
assert test('../bin/wiggletools pearson fixedStep.wig variableStep.wig > tmp/pearson.txt') == 0

# Testing profiles
assert test('../bin/wiggletools profiles 3 overlapping.bed fixedStep.wig > tmp/profiles.txt') == 0

# Testing profile
assert test('../bin/wiggletools profile 3 overlapping.bed fixedStep.wig > tmp/profile.txt') == 0

assert test('diff tmp expected') == 0
shutil.rmtree('tmp')

print 'All tests OK'
