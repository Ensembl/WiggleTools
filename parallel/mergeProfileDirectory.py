#!/usr/bin/env python

import sys
import glob
import shutil

target = sys.argv[1]
firstFile = True
array = None

# Scan through files and add up their profiles
for file in glob.glob(target + "x/*"):
	if not firstFile:
		index = 0
	
	for line in file:
		if firstFile and array is None:
			array = float(line.strip())
		elif firstFile:
			array.append(float(line.strip()))
		else:
			array[index] += float(line.strip())
			index += 1

	firstFile = False

# Write out results
file = open(target, 'w')
for value in array:
	file.write("%f\n", value)
file.close()

# Remove unnecessary files
shutil.rmtree(target + "x")
