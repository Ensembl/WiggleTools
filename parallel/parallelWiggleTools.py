#!/usr/bin/env python

# Copyright 2013 EMBL-EBI
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

# Directory where job stdin, stdout and stderr are stored
# Must be visible to all LSF nodes
# By default, right on your doorstep ;P
DUMP_DIR = '.'

################################################
## Splitting a wiggletools command into regional jobs
################################################

def create_dirs(command):
	for match in re.finditer(r'write\s*(\S*.bw)\s', command):
		path = match.group(1) + "x"
		if not os.path.exists(path):
			os.makedirs(path)

def create_new_command(command, chr, start, finish, chrom_sizes_file):
	command = re.sub(r'write\s*(\S*.bw)\s',r'write \1x/%s_%i_%i.wig ' % (chr, start, finish), command)
	command = re.sub(r'(apply|profile|profiles)\s*(\S*)\s',r'\1 \2x/%s_%i_%i ' % (chr, start, finish), command)
	command = re.sub(r'^(AUC|mean|variance|pearson)\s*(\S*)\s',r'\1 \2x/%s_%i_%i ' % (chr, start, finish), command)
	return " ".join(map(str, ['wiggletoolsIndex.py', chrom_sizes_file, 'do','seek',chr,start,finish,command]))

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

def makeReduceCommand(command):
	mergeBigWigCommands = ['mergeBigWigDirectory.py %s' % match.group(1) for match in re.finditer(r'write\s*(\S*.bw)\s', command)]
	mergeApplyCommand = ['mergeApplyDirectory.sh %s' % match.group(1) for match in re.finditer(r'apply\s*(\S*)\s', command)]
	mergeProfileCommand = ['mergeProfileDirectory.py %s' % match.group(1) for match in re.finditer(r'profile\s*(\S*)\s', command)]
	mergeProfilesCommand = ['mergeProfilesDirectory.sh %s' % match.group(1) for match in re.finditer(r'profiles\s*(\S*)\s', command)]
	return mergeBigWigCommands + mergeApplyCommand + mergeProfileCommand + mergeProfilesCommand

################################################
## LSF MultiJob
################################################

def makeMultiJobCommand(filename, count, dependency=None, mem=4):
	name = os.path.basename(filename)
	bsub_cmd = "bsub -q normal -R'select[mem>%i] rusage[mem=%i]' -M%i -J'%s[1-%s]'" % (1024*mem, 1024*mem, 1024*mem, name, count)
	if dependency is not None:
		bsub_cmd += " -w '%s[*]'" % dependency
	output = "-o %s_%%I.out -e %s_%%I.err" % (filename, filename)
	jobCmd = " ".join([bsub_cmd, output, 'LSFwrapper.sh', "' multiJob.py ", filename, "'"])
	print jobCmd
	return jobCmd

def submitMultiJobToLSF(cmds, dependency=None, mem=4):
	descr, filename = tempfile.mkstemp(dir='.')

	fh = open(filename, 'w')
	fh.write("\n".join(cmds))
	fh.close()

	multi_job_cmd = makeMultiJobCommand(filename, len(cmds), dependency, mem)
	p = subprocess.Popen(multi_job_cmd, shell=True, stdout=subprocess.PIPE)
	err = p.wait()
	if err != 0:
		print "Could not start job:"
		print multi_job_cmd
		sys.exit(err)

	out, err = p.communicate()

	for line in out.split('\n'):
		match = re.match(r'Job <([0-9]*)>', line)
		if match is not None:
			return match.group(1)

	raise RuntimeError


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

def main():
	if len(sys.argv) == 2:
		chrom_file = sys.argv[1]
		chrom_sizes = readChromSizes(chrom_file)
		cmd = sys.argv[2]
		jobID = submitMultiJobToLSF(makeMapCommand(cmd, chrom_file, chrom_sizes, region_size=3e7))
		submitMultiJobToLSF(makeReduceCommand(cmd), dependency = jobID, mem=8)
	else:
		print """
parallelWiggletools.py: wrapper script to run wiggletools in parallel on LSF

Usage: parallelWiggletools.py chrom_sizes.txt ' command '

Where:
chrom_sizes.txt is a tab-delimited text file with the chromosome names and lengths	
command is a valid wiggletools command, between single quotes.
		"""

if __name__ == "__main__":
	main()
