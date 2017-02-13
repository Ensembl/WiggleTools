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

def makeCommand(filename, count, batch_system='LSF', dependency=None, mem=4, working_directory='.'):
	name = os.path.basename(filename)
	if batch_system == 'LSF':
		bsub_cmd = "bsub -q normal -R'select[mem>%i] rusage[mem=%i]' -M%i -J'%s[1-%s]'" % (1024*mem, 1024*mem, 1024*mem, name, count)
		if dependency is not None:
			bsub_cmd += " -w '%s[*]'" % dependency
		output = "-o %s/%s_%%I.out -e %s/%s_%%I.err" % (working_directory, filename, working_directory, filename)
		jobCmd = " ".join([bsub_cmd, output, 'LSFwrapper.sh', "' multiJob.py ", filename, batch_system, "'"])
	elif batch_system == 'SGE':
		bsub_cmd = "qsub -terse -cwd -V -b y -t 1-%s -N %s" % (count, os.path.basename(filename))
		if dependency is not None:
			bsub_cmd += " -hold_jid %s" % dependency
		output = "-o %s -e %s" % (working_directory, working_directory)
		jobCmd = " ".join([bsub_cmd, output, "' multiJob.py ", filename, batch_system, "'"])
	elif batch_system == 'local':
		jobCmd = "sh " + filename + ">& " + os.path.join(working_directory, filename + ".oe");
	else:
		raise NameError

	return jobCmd

def submit(cmds, batch_system="LSF", dependency=None, mem=4, working_directory='.'):
	if len(cmds) == 0:
		sys.stderr.write("No commands in list")
		raise RuntimeError
		return None, None

	descr, filename = tempfile.mkstemp(dir=working_directory)

	fh = open(filename, 'w')
	fh.write("\n".join(cmds))
	fh.close()

	multi_job_cmd = makeCommand(filename, len(cmds), batch_system, dependency, mem, working_directory)
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

	if batch_system == 'LSF':
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
	batch_system = sys.argv[2]

	if batch_system == 'LSF':
		index = os.environ['LSB_JOBINDEX']
	elif batch_system == 'SGE':
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
