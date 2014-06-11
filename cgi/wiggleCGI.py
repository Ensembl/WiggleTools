#!/usr/bin/env python

# SQLite3 file location
database_location = "/data/wiggletools/datasets.sqlite3"
# Empty directory to catch CGI log files
logdir = '/data/wiggletools/log/'
# Directory to catch results files
working_directory = '/data/tmp/'
# S3 references
s3_bucket = 'wiggletools-data'
s3_region = 'eu-west-1'

import sys
import cgi
import cgitb
import json
import sqlite3
import re
import wiggletools.wiggleDB

cgitb.enable(logdir=logdir)

class WiggleDBOptions(object):
	def __init__(self):
		self.conn = None
		self.assembly = None
		self.wa  = None
		self.working_directory = None
		self.s3 = None
		self.wb = None
		self.a = None
		self.b = None

def main():
	print "Content-Type: application/json"
	print

	try:
		form = cgi.FieldStorage()
		conn = sqlite3.connect(database_location)
		cursor = conn.cursor()
		if "result" in form:
			status, info = wiggletools.wiggleDB.query_result(cursor, form["result"].value)
			if status == "DONE":
				base_url = 's3-%s.amazonaws.com/%s/' % (s3_region, s3_bucket)
				url = re.sub(working_directory, base_url, info)		
				print json.dumps({'status':status, 'url':url})
			else:
				print json.dumps({'status':status})

		elif "count" in form:
			assembly = form['assembly'].value
			params = dict((re.sub("^._", "", X), form.getlist(X)) for X in form if X != "count")
			count = len(wiggletools.wiggleDB.get_dataset_locations(cursor, params, assembly))
			print json.dumps({'query':params,'count':count})

		elif 'annotations' in form:
			assembly = form['assembly'].value
			print json.dumps({"annotations": wiggletools.wiggleDB.get_annotations(cursor, assembly)})

		elif 'wa' in form:
			options = WiggleDBOptions()
			options.conn = conn
			options.assembly = form['assembly'].value
			options.wa = form['wa'].value
			options.working_directory = working_directory
			options.s3 = s3_bucket

			if 'wb' in form:
				options.wb = form['wb'].value
			else:
				options.wb = None

			if 'merge' in form:
				options.fun_merge = form['merge'].value
			else:
				options.fun_merge = None

			options.a = dict((X[2:], form.getlist(X)) for X in form if X[:2] == "A_")
			options.b = dict((X[2:], form.getlist(X)) for X in form if X[:2] == "B_")
			
			print json.dumps({'ID':1})
			#jobID = wiggletools.wiggleDB.request_compute(cursor, options)
			#print json.dumps({'ID':jobID, 'options':options})

		else:
			print json.dumps("No params, no output")

		conn.commit()
		conn.close()
        except:
                print json.dumps("ERROR")
		raise

if __name__ == "__main__":
	main()
