#!/usr/bin/env python

# Copyright [1999-2017] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import sys
import subprocess
import wigglePlots
import os
import os.path

file=sys.argv[1]
if (re.call("sort -k1,1 -k2,2n -k3,3n -m %sx/* > %s" % (file, file), shell=True)):
	print 'Error processing directory %sx' % file
	sys.exit(100)
os.rmdir(file + "x")

try:
	if os.path.getsize(wiggle_file) > 0:
		wigglePlots.make_overlaps(file, file + ".png")
	else:
		# Create empty file with .empty suffix
		open(file + ".empty", "w").close()
except:
	sys.exit(100)
