#!/usr/bin/env python

import sys
import subprocess
import wigglePlots
import os
import shutil

file=sys.argv[1]
if (subprocess.call("sort -k1,1 -k2,2n -k3,3n -m %sx/* > %s" % (file, file), shell=True)):
	print 'Error processing directory %sx' % file
	sys.exit(1)
shutil.rmtree(file + "x")
wigglePlots.make_profiles_matrix(file, file + ".pdf")
