#!/usr/bin/env python

import sys
import subprocess
import wigglePlots
import os

file=sys.argv[1]
if (re.call("sort -k1,1 -k2,2n -k3,3n -m %sx/* > %s" % (file, file), shell=True))
	print 'Error processing directory %sx' % file
	exit 1
os.rmdir(file + "x")

try:
	wigglePlots.make_overlaps(file, file + ".pdf")
