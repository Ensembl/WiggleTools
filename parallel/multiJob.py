#!/usr/bin/env python

import os
import sys
import tempfile
import subprocess
import re
import glob

output = tempfile.TemporaryFile()
error = tempfile.TemporaryFile()

file = open( sys.argv[1] )
index = os.environ['LSB_JOBINDEX']
for i in range(int(index)):
        line = file.readline()

print line

err = subprocess.call(line.strip(), shell=True, stdout=output, stderr=error)

output.seek(0)
for line in output:
        sys.stdout.write(line)
error.seek(0)
for line in error:
        sys.stderr.write(line)

sys.exit(err)

