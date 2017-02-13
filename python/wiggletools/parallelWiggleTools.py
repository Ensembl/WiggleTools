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
import os.path
import re
import tempfile
import subprocess
import glob
import multiJob

# Directory where job stdin, stdout and stderr are stored
# Must be visible to all LSF nodes
# By default, right on your doorstep ;P
DUMP_DIR = '.'

################################################
## Splitting a wiggletools command into regional jobs
################################################

def create_dirs(command):
	for match in re.finditer(r'(write|write_bg|profile|profiles|AUC|mean|variance|pearson)\s+(\S+)\s', command):
		path = match.group(2) + "x"
		if not os.path.exists(path):
			os.makedirs(path)

def create_new_command(command, chr, start, finish, chrom_sizes_file):
	command = re.sub(r'write\s+(\S+.wig)\s',r'write \1x/%s_%i_%i.wig ' % (chr, start, finish), command)
	# Careful: the following line has to be AFTER the one above, else they overwrite each other
	command = re.sub(r'write\s+(\S+.bw)\s',r'write \1x/%s_%i_%i.wig ' % (chr, start, finish), command)
	command = re.sub(r'write_bg\s+(\S+)\s',r'write \1x/%s_%i_%i.wig ' % (chr, start, finish), command)
	command = re.sub(r'(profile|profiles)\s+(\S+)\s',r'\1 \2x/%s_%i_%i ' % (chr, start, finish), command)
	command = re.sub(r'^(AUC|mean|variance|pearson)\s+(\S+)\s',r'\1 \2x/%s_%i_%i ' % (chr, start, finish), command)
	
	m = re.match(r'(profile|profiles)\s+(\S+)\s(\S+)\s+(.*)', command)
	m2 = re.match(r'(AUC|mean|variance|pearson)\s+(\S+)\s+(.*)', command)
	if m is not None:
		plot = m.group(1)
		output = m.group(2)
		width = m.group(3)
		iterator = m.group(4)
		return " ".join(map(str, ['wiggletoolsIndex.py', chrom_sizes_file, plot, output, width,'seek',chr,start,finish,iterator]))
	elif m2 is not None:
		plot = m.group(1)
		output = m.group(2)
		iterators = m.group(3)
		return " ".join(map(str, ['wiggletoolsIndex.py', chrom_sizes_file, plot, output,'seek',chr,start,finish,iterators]))
	else :
		# Command does not start with extraction function:
		return " ".join(map(str, ['wiggletoolsIndex.py', chrom_sizes_file,'do','seek',chr,start,finish,command]))


def makeMapCommand(command, chrom_sizes_file, chrom_sizes, region_size):
	create_dirs(command)
	return [create_new_command(command, chr, start, min(chrom_sizes[chr], start + int(region_size)), chrom_sizes_file) for chr in sorted(chrom_sizes.keys()) for start in range(1, chrom_sizes[chr], int(region_size))]

def test_makeMapCommand():
	chrom_sizes = dict([("chr1", 20), ("chr2", 30)])
	chrom_sizes_file = 'chrom_sizes.txt'
	cmd = 'add toto.bw tata.bam'
	destination = 'sum.bg'
	print makeMapCommand(cmd, chrom_sizes_file, chrom_sizes, region_size=int(3e7))

################################################
## Merging the results of parallel wiggletools runs 
################################################

def makeReduceCommands(command):
	mergeBigWigCommands = ['mergeBigWigDirectory.py %s' % match.group(1) for match in re.finditer(r'write\s+(\S+.bw)\s', command)]
	mergeProfileCommand = ['mergeProfileDirectory.py %s' % match.group(1) for match in re.finditer(r'profile\s+(\S+)\s', command)]
	mergeProfilesCommand = ['mergeProfilesDirectory.py %s' % match.group(1) for match in re.finditer(r'profiles\s+(\S+)\s', command)]
	mergeWigglesCommand = ['mergeBedLikeDirectory.sh %s' % match.group(1) for match in re.finditer(r'write\s+(\S+.wig)\s', command)]
	mergeBedGraphsCommand = ['mergeBedLikeDirectory.sh %s' % match.group(1) for match in re.finditer(r'write_bg\s+(\S+.bg)\s', command)]
	return mergeBigWigCommands + mergeProfileCommand + mergeProfilesCommand + mergeWigglesCommand + mergeBedGraphsCommand


################################################
## Main function
################################################

def readChromSizes(file):
	chrom_sizes = dict()
	fh = open(file)
	for line in fh:
		items = line.strip().split()
		chrom_sizes[items[0]] = int(items[1])
	fh.close()
	return chrom_sizes

def run(cmds, chrom_file, batch_system='local', tmp='.'):
	for cmd in cmds:
		if re.search('(apply_paste|histogram)', cmd) is not None:
			print "Cannot parallelize the computation of histograms or apply_paste operations"
			sys.exit(1)
	chrom_sizes = readChromSizes(chrom_file)
	mapCommands = sum((makeMapCommand(cmd, chrom_file, chrom_sizes, region_size=3e7) for cmd in cmds), [])
	jobID1, filename1 = multiJob.submit(mapCommands, batch_system=batch_system, working_directory=tmp)
	reduceCommands = sum((makeReduceCommands(cmd) for cmd in cmds), [])
	jobID2, filename2 = multiJob.submit(reduceCommands, dependency = jobID1, batch_system=batch_system, mem=8, working_directory=tmp)
	return jobID2, [filename1, filename2]

def main():
	if len(sys.argv) == 3:
		chrom_file = sys.argv[1]
		cmds = sys.argv[2:]
		run(cmds, chrom_file)
	else:
		print """
parallelWiggletools.py: wrapper script to run wiggletools in parallel on LSF

Usage: parallelWiggletools.py chrom_sizes.txt 'command1' ['command2' [...]]

Where:
* chrom_sizes.txt is a tab-delimited text file with the chromosome names and lengths	
* command* is a valid wiggletools command, without histogram or apply_paste keywords.
		"""

if __name__ == "__main__":
	main()
