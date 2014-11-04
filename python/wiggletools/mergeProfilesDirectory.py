#!/usr/bin/env python

import sys
import subprocess
import wigglePlots
import os
import os.path
import shutil

try:
	file=sys.argv[1]
	if (subprocess.call("sort -k1,1 -k2,2n -k3,3n -m %sx/* > %s" % (file, file), shell=True)):
		print 'Error processing directory %sx' % file
		sys.exit(100)
	shutil.rmtree(file + "x")

	if os.path.getsize(file) > 0:
		wigglePlots.make_profiles_matrix(file, file + ".png", format='png')
	else:
		# Create empty file with .empty suffix
		open(file + ".empty", "w").close()
except:
	sys.exit(100)
