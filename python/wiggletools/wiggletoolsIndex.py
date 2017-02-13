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
import os
import re
import os.path

chrom_lengths_file = sys.argv[1]
wiggletools_args = sys.argv[2:]

def run(args):
	p = subprocess.Popen(" ".join(args), stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	err = p.wait()
	stdout, stderr = p.communicate()
	sys.stdout.write(stdout)
	sys.stderr.write(stderr)
	if err != 0:
		sys.exit(100)

def main():
	cmd = ['wiggletools'] + wiggletools_args
	run(cmd)

	for i in range(1, len(wiggletools_args)):
		if wiggletools_args[i-1] == 'write' and re.search(r'\S*.wig', wiggletools_args[i]) is not None:
			wiggle_file = wiggletools_args[i]
			if os.path.getsize(wiggle_file) > 0:
				bigwig_file = re.sub('.wig$','.bw', wiggle_file)
				run(['wigToBigWig -keepAllChromosomes -fixedSummaries',wiggle_file,chrom_lengths_file,bigwig_file])
			os.remove(wiggle_file)

if __name__ == "__main__":
	main()
