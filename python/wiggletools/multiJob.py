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

def makeCommand(filename, count, batchSystem='LSF', dependency=None, mem=4, working_directory='.'):
	name = os.path.basename(filename)
	if batchSystem == 'LSF':
		bsub_cmd = "bsub -q normal -R'select[mem>%i] rusage[mem=%i]' -M%i -J'%s[1-%s]'" % (1024*mem, 1024*mem, 1024*mem, name, count)
		if dependency is not None:
			bsub_cmd += " -w '%s[*]'" % dependency
		output = "-o %s/%s_%%I.out -e %s/%s_%%I.err" % (working_directory, filename, working_directory, filename)
		jobCmd = " ".join([bsub_cmd, output, 'LSFwrapper.sh', "' multiJob.py ", filename, batchSystem, "'"])
	elif batchSystem == 'SGE':
		bsub_cmd = "qsub -terse -cwd -V -b y -t 1-%s -N %s" % (count, os.path.basename(filename))
		if dependency is not None:
			bsub_cmd += " -hold_jid %s" % dependency
		output = "-o %s -e %s" % (working_directory, working_directory)
		jobCmd = " ".join([bsub_cmd, output, "' multiJob.py ", filename, batchSystem, "'"])
	elif batchSystem == 'local':
		jobCmd = "sh " + filename + ">& " + os.path.join(working_directory, filename + ".oe");
	else:
		raise NameError

	return jobCmd

def submit(cmds, batchSystem="LSF", dependency=None, mem=4, working_directory='.'):
	if len(cmds) == 0:
		sys.stderr.write("No commands in list")
		raise RuntimeError
		return None, None

	descr, filename = tempfile.mkstemp(dir=working_directory)

	fh = open(filename, 'w')
	fh.write("\n".join(cmds))
	fh.close()

	multi_job_cmd = makeCommand(filename, len(cmds), batchSystem, dependency, mem, working_directory)
	p = subprocess.Popen(multi_job_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	ret = p.wait()
	out, err = p.communicate()
	if ret != 0:
		print "Could not start job:"
		print multi_job_cmd
		print "OUT: " + out
		print "ERR: " + err
		print os.environ
		if 'SGE_ROOT' in os.environ:
			print "SGE_ROOT: " + os.environ['SGE_ROOT']
		else:
			print "SGE_ROOT: UNDEF" 
  		print 'USER: ' + os.environ['USERNAME']
		assert False

	if batchSystem == 'LSF':
		for line in out.split('\n'):
			match = re.match(r'Job <([0-9]*)>', line)
			if match is not None:
				return match.group(1), filename
		sys.stderr.write("Could not find job id in lsf output: %s" % out)
		raise RuntimeError
		return None, None
	else:
		return re.split(r'[\n\.]', out)[0], filename


################################################
## Worker 
################################################

def main():
	output = tempfile.TemporaryFile()
	error = tempfile.TemporaryFile()

	file = open( sys.argv[1] )
	batchSystem = sys.argv[2]

	if batchSystem == 'LSF':
		index = os.environ['LSB_JOBINDEX']
	elif batchSystem == 'SGE':
		index = os.environ['SGE_TASK_ID']
	else:
		raise NameError

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
