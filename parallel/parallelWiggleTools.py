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

################################################
## Splitting a wiggletools command into regional jobs
################################################

def create_dirs(command):
	for match in re.finditer(r'write\s*(\S*.bw)\s', command):
		path = match.group(1) + "x"
		if not os.path.exists(path):
			os.makedirs(path)

def create_new_command(command, chr, start, finish, chrom_sizes_file):
	substitute_cmd = re.sub(r'write\s*(\S*.bw)\s',r'write \1x/%s_%i_%i.wig ' % (chr, start, finish), command)
	return " ".join(map(str, ['wiggletoolsIndex.py', chrom_sizes_file, "'", 'do','seek',chr,start,finish,substitute_cmd, "'"]))

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
## LSF MultiJob
################################################

def makeMultiJobCommand(filename, count):
	name = os.path.basename(filename)
	bsub_cmd = "bsub -q normal -R'select[mem>4000] rusage[mem=4000]' -M4000 -J%s[1-%s] " % (name, count)
	output = "-o %s_%%I.out -e %s_%%I.err" % (filename, filename)
	jobCmd = " ".join([bsub_cmd, output, 'LSFwrapper.sh', "' multiJob.py ", filename, "'"])
	print jobCmd
	return jobCmd

def submitMultiJobToLSF(cmds):
	descr, filename = tempfile.mkstemp(dir="/lustre/scratch109/ensembl/dz1")

	fh = open(filename, 'w')
	fh.write("\n".join(cmds))
	fh.close()

	p = subprocess.Popen(makeMultiJobCommand(filename, len(cmds)), shell=True)
	err = p.wait()
	if err != 0:
		print "Could not start job %s" % cmds
		print line
		sys.exit(err)

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
	chrom_file = sys.argv[1]
	chrom_sizes = readChromSizes(chrom_file)
	cmd = sys.argv[2]
	submitMultiJobToLSF(makeMapCommand(cmd, chrom_file, chrom_sizes, region_size=3e8))

if __name__ == "__main__":
	main()
