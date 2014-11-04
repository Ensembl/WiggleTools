#!/usr/bin/env python

import sys
import argparse
import subprocess
import glob
import os
import os.path

import wiggletools.wiggleDB
import wiggletools.multiJob 
import wiggletools.wigglePlots

def get_options():
	parser = argparse.ArgumentParser(description='WiggleDB worker.')
	parser.add_argument('--db', '-d', dest='db_name', help='Database file', required=True)
	parser.add_argument('--jobID','-i',dest='jobID', help='Job ID number', type=int, required=True)
	parser.add_argument('--temp','-t',dest='temps', help='Temporary files to be cleaned up', nargs='*')
	parser.add_argument('--histogram','-p',dest='histogram', help='Draw histogram')
	parser.add_argument('--apply_paste','-v',dest='apply_paste', help='Run Apply_paste command')
	parser.add_argument('--labels','-l',dest='labels', help='Labels for plot categories', nargs='*')
	parser.add_argument('--emails','-e',dest='emails', help='Notification e-mail addresses', nargs='*')
	parser.add_argument('--config','-c',dest='conf', help='Config file')
	return parser.parse_args()

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
		options = get_options()
		config = wiggletools.wiggleDB.read_config_file(options.conf)
		data = wiggletools.wiggleDB.get_job_location(options.db_name, options.jobID)
		empty = os.path.exists(data + ".empty")

		# Optional graphics
		if options.histogram is not None:
			if subprocess.call("wiggletools "  + options.histogram, shell=True):
				print "Failed to construct histogram"
				sys.exit(1)
			if os.path.getsize(data) > 0:
				wiggletools.wigglePlots.make_histogram(data, options.labels, data + ".png", format='png')
			else:
				empty = True

		if options.apply_paste is not None:
			if subprocess.call("wiggletools "  + options.apply_paste, shell=True):
				print "Failed to construct overlap graph"
				sys.exit(1)
			if os.path.getsize(data) > 0:
				wiggletools.wigglePlots.make_overlaps(data, data + ".png", format='png')
			else:
				empty = True

		# Housekeeping
		if options.temps is not None:
			wiggletools.multiJob.clean_temp_files(options.temps)

		# Signing off
		if empty:
			wiggletools.wiggleDB.report_empty_to_user(options.emails, options.jobID, config)
			wiggletools.wiggleDB.mark_job_status(options.db_name, options.jobID, 'EMPTY')
		else:
			copy_to_longterm(data, config)
			if os.path.exists(data + ".png"):
				copy_to_longterm(data + ".png", config)
			wiggletools.wiggleDB.report_to_user(options.emails, options.jobID, data, config)
			wiggletools.wiggleDB.mark_job_status(options.db_name, options.jobID, 'DONE')
	except:
		raise
		sys.exit(100)

if __name__ == "__main__":
	main()
