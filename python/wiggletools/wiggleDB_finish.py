#!/usr/bin/env python

import sys
import argparse
import subprocess
import glob
import os
import os.path

import wiggleDB
import multiJob 
import wigglePlots

def get_options():
	parser = argparse.ArgumentParser(description='WiggleDB worker.')
	parser.add_argument('--db', '-d', dest='db_name', help='Database file', required=True)
	parser.add_argument('--jobID','-i',dest='jobID', help='Job ID number', type=int, required=True)
	parser.add_argument('--temp','-t',dest='temps', help='Temporary files to be cleaned up', nargs='*', required=True)
	parser.add_argument('--histogram','-p',dest='histogram', help='Draw histogram')
	parser.add_argument('--apply_paste','-v',dest='apply_paste', help='Run Apply_paste command')
	parser.add_argument('--labels','-l',dest='labels', help='Labels for plot categories', nargs='*')
	parser.add_argument('--email','-e',dest='emails', help='Notification e-mail addresses', nargs='*')
	parser.add_argument('--s3',dest='s3', help='S3 bucket')
	return parser.parse_args()

def main():
	options = get_options()
	data = wiggleDB.get_job_location(options.db_name, options.jobID)
	if options.temps is not None:
		multiJob.clean_temp_files(options.temps)
	if options.histogram is not None:
		if subprocess.call("wiggletools "  + options.histogram, shell=True):
			print "Failed to construct histogram"
			sys.exit(1)
		wigglePlots.make_histogram(data, options.labels, data + ".pdf")
	if options.apply_paste is not None:
		if subprocess.call("wiggletools "  + options.apply_paste, shell=True):
			print "Failed to construct overlap graph"
			sys.exit(1)
		wigglePlots.make_overlaps(data, data + ".pdf")
	if options.s3 is not None:
		cmd = "aws s3 cp %s s3://%s/%s --acl public-read" % (data, options.s3, os.path.basename(data))
		if subprocess.call(cmd, shell=True):
			print "Failed to copy over results"
			print cmd
			sys.exit(1)
	wiggleDB.mark_job_as_done(options.db_name, options.jobID)

if __name__ == "__main__":
	main()
