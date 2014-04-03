#!/usr/bin/env python

import os
import os.path
import sys
import tempfile
import subprocess
import re
import glob

################################################
## File hygiene
################################################

def clean_temp_file(file):
	if file is None:
		return 

	if os.path.exists(file):
		os.remove(file)

	for file2 in glob.glob("%s_[0-9]*.out" % file) + glob.glob("%s_[0-9]*.err" % file):
		os.remove(file2)

def clean_temp_files(files):
	map(clean_temp_file, files)

################################################
## LSF MultiJob
################################################

def makeCommand(filename, count, dependency=None, mem=4):
	name = os.path.basename(filename)
	bsub_cmd = "bsub -q normal -R'select[mem>%i] rusage[mem=%i]' -M%i -J'%s[1-%s]'" % (1024*mem, 1024*mem, 1024*mem, name, count)
	if dependency is not None:
		bsub_cmd += " -w '%s[*]'" % dependency
	output = "-o %s_%%I.out -e %s_%%I.err" % (filename, filename)
	jobCmd = " ".join([bsub_cmd, output, 'LSFwrapper.sh', "' multiJob.py ", filename, "'"])
	print jobCmd
	return jobCmd

def submit(cmds, dependency=None, mem=4):
	if len(cmds) == 0:
		sys.stderr.write("No commands in list")
		raise RuntimeError
		return None, None

	descr, filename = tempfile.mkstemp(dir='.')

	fh = open(filename, 'w')
	fh.write("\n".join(cmds))
	fh.close()

	multi_job_cmd = makeCommand(filename, len(cmds), dependency, mem)
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
			return match.group(1), filename

	sys.stderr.write("Could not find job id in lsf output: %s" % out)
	raise RuntimeError
	return None, None

################################################
## Worker 
################################################

def main():
	output = tempfile.TemporaryFile()
	error = tempfile.TemporaryFile()

	file = open( sys.argv[1] )
	index = os.environ['LSB_JOBINDEX']
	for i in range(int(index)):
		line = file.readline()

	print line

	err = subprocess.call(line.strip(), shell=True, stdout=output, stderr=error)

	output.seek(0)
	for line in output:
		sys.stdout.write(line)
	error.seek(0)
	for line in error:
		sys.stderr.write(line)

	sys.exit(err)

if __name__ == "__main__":
	main()
