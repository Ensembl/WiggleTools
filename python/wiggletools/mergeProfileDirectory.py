#!/usr/bin/env python
# Copyright 2014 EMBL-EBI                                                                             
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
