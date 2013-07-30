#!/usr/bin/env python

import sys
import subprocess
import os
import re

chrom_lengths_file = sys.argv[1]
wiggletools_args = sys.argv[2]

def run(args):
	p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	err = p.wait()
	if err != 0:
		stdout, stderr = p.communicate()
		sys.stdout.write(stdout)
		sys.stderr.write(stderr)
		sys.exit(err)

def main():
	run(['wiggletools', wiggletools_args])

	for match in re.finditer('write\W*(\w*.wig)\W', wiggletools_args):
		wiggle_file = match.group(1)
		bigwig_file = re.sub('.wig$','.bw', wiggle_file)
		run(['wigToBigWig',wiggle_file,chrom_lengths_file,bigwig_file])
		os.remove(wiggle_file)

if __name__ == "__main__":
	main()
