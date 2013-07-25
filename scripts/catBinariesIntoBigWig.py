#!/usr/bin/env python

import sys
import subprocess
import os
import tempfile

fd, path = tempfile.mkstemp()

p = subprocess.Popen(["wiggletools", 'write %s cat %s' % (path, " ".join(sys.argv[3:]))])
err = p.wait()
if err != 0:
	print "Concatenation failed!"
	stdout, stderr = p.communicate()
	if stdout is not None:
		sys.stdout.write(stdout)
	if stderr is not None:
		sys.stderr.write(stderr)
	sys.exit(err)

p = subprocess.Popen(["bedGraphToBigWig", path, sys.argv[1], sys.argv[2]])
err = p.wait()
if err != 0:
	print "indexing failed!"
	stdout, stderr = p.communicate()
	if stdout is not None:
		sys.stdout.write(stdout)
	if stderr is not None:
		sys.stderr.write(stderr)
	sys.exit(err)

os.remove(path)
for file in sys.argv[3:]:
	os.remove(file)
