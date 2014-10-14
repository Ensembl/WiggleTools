#!/usr/bin/env python

import sys
import subprocess
import glob
import shutil
import os

directory = sys.argv[1]
bigwigs = glob.glob(os.path.join(directory + "x", '*.bw'))

if len(bigwigs) > 0:
	command = ['catBigWigs', directory] + bigwigs
	p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	return_code = p.wait()

	if return_code != 0:
		out, err = p.communicate()
		print "Failed: %s" % " ".join(command)
		print "Stdout:"
		print out
		print "Stderr:"
		print err
		sys.exit(1)

else:
	# Create empty file with .empty suffix
	open(directory + ".empty", "w").close()

shutil.rmtree(directory + "x")
