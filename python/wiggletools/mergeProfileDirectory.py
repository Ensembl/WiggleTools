#!/usr/bin/env python

import sys
import glob
import shutil
import wigglePlots

target = sys.argv[1]
firstFile = True
array = None

# Scan through files and add up their profiles
for file in glob.glob(target + "x/*"):
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

# Write out results
file = open(target, 'w')
for index in range(len(array)):
	file.write("%i\t%f\n" % (index, array[index]))
file.close()

# Remove unnecessary files
shutil.rmtree(target + "x")

# Do a little drawing
wigglePlots.make_profile_curve(target, target + ".pdf")
