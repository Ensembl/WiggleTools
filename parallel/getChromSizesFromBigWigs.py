#!/usr/bin/env python

import sys
import re
import subprocess

regexp = re.compile("\W*(\w*) \d* (\w*)")

chrom_lengths = dict()

for file in sys.argv[1:]:
	cmd = 'bigWigInfo -chroms ' + file
	p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
	for line in p.communicate()[0].split("\n"):
		m = regexp.match(line)
		if m is not None:
			chrom_lengths[m.group(1)] = m.group(2)

for chrom in sorted(chrom_lengths.keys()):
	print "\t".join([chrom, chrom_lengths[chrom]])

