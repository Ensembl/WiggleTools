#!/usr/bin/env python

import sys
import subprocess
import os
import os.path
import json

import wiggletools.wiggleDB
import wiggletools.multiJob 
import wiggletools.wigglePlots

class Struct(object):
        def __init__(self, **entries):
            self.__dict__.update(entries)

def get_options():
	assert len(sys.argv) == 2
	options = Struct(**(json.load(open(sys.argv[-1]))))
	return options, wiggletools.wiggleDB.read_config_file(options.config)

def copy_to_longterm(data, config):
	if 's3_bucket' in config:
		os.environ['AWS_CONFIG_FILE'] = config['aws_config']
		cmd = "aws s3 cp %s s3://%s/%s --acl public-read" % (data, config['s3_bucket'], os.path.basename(data))
		if subprocess.call(cmd, shell=True) != 0:
			print "Failed to copy over results"
			print cmd
			sys.exit(100)

def main():
	try:
		options, config = get_options()
		empty = os.path.exists(options.data + ".empty")

		# Optional graphics
		if options.histogram is not None:
			if subprocess.call("wiggletools "  + options.histogram, shell=True):
				print "Failed to construct histogram"
				sys.exit(1)
			if os.path.getsize(options.data) > 0:
				wiggletools.wigglePlots.make_histogram(options.data, options.labels, options.data + ".png", format='png')
			else:
				empty = True

		if options.apply_paste is not None:
			if subprocess.call("wiggletools "  + options.apply_paste, shell=True):
				print "Failed to construct overlap graph"
				sys.exit(1)
			if os.path.getsize(options.data) > 0:
				wiggletools.wigglePlots.make_overlaps(options.data, options.data + ".png", format='png')
			else:
				empty = True

		# Signing off
		if empty:
			wiggletools.wiggleDB.report_empty_to_user(options, config)
			wiggletools.wiggleDB.mark_job_status(options.db, options.jobID, 'EMPTY')
		else:
			copy_to_longterm(options.data, config)
			if os.path.exists(options.data + ".png"):
				copy_to_longterm(data + ".png", config)
			wiggletools.wiggleDB.report_to_user(options, config)
			wiggletools.wiggleDB.mark_job_status(options.db, options.jobID, 'DONE')

		# Housekeeping
		if options.temps is not None:
			wiggletools.multiJob.clean_temp_files(options.temps)
		os.remove(sys.argv[-1])
	except:
		raise
		sys.exit(100)

if __name__ == "__main__":
	main()
