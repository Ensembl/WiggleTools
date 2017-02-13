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
import glob
import shutil
import wigglePlots

try:
	target = sys.argv[1]
	firstFile = True
	array = None

	files = glob.glob(target + "x/*")

	if len(files) > 0:
		# Scan through files and add up their profiles
		for file in files:
			for line in open(file):
				items = line.strip().split('\t')
				index = int(items[0])
				value = float(items[1])
				if firstFile and array is None:
					array = [float(value)]
				elif firstFile:
					array.append(float(value))
				else:
					array[index] += float(value)

			firstFile = False

		if array is None:
			# Create empty file with .empty suffix
			open(target + ".empty", "w").close()
		else:
			# Write out results
			file = open(target, 'w')
			for index in range(len(array)):
				file.write("%i\t%f\n" % (index, array[index]))
			file.close()

			# Do a little drawing
			wigglePlots.make_profile_curve(target, target + ".png", format='png')
	else:
		# Create empty file with .empty suffix
		open(target + ".empty", "w").close()

	# Remove unnecessary files
	shutil.rmtree(target + "x")
except:
	sys.exit(100)
