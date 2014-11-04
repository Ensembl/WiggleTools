#!/usr/bin/env python

import sys
import cgi
import cgitb
import json
import sqlite3
import re
import wiggletools.wiggleDB

DEBUG = False
CONFIG_FILE = '/data/wiggletools/wiggletools.conf'

config = wiggletools.wiggleDB.read_config_file(CONFIG_FILE)
cgitb.enable(logdir=config['logdir'])

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
		self.dry_run = DEBUG
		self.remember = False
		self.db = config['database_location']
		self.config = CONFIG_FILE
		self.emails = None
		

def report_result(result):
	base_url = 'http://s3-%s.amazonaws.com/%s/' % (config['s3_region'], config['s3_bucket'])
	url = re.sub(config['working_directory'], base_url, result['location'])		
	if result['location'][-3:] == ".bw" or result['location'][-3:] == ".bb":
		ensembl = 'http://%s/%s/Location/View?g=%s;contigviewbottom=url:%s' % (config['ensembl_server'], config['ensembl_species'], config['ensembl_gene'], url)
	else:
		ensembl = url + ".png"
	print json.dumps({'status':result['status'], 'url':url, 'view':ensembl})

def main():
	print "Content-Type: application/json"
	print

	try:
		form = cgi.FieldStorage()
		conn = sqlite3.connect(config['database_location'])
		cursor = conn.cursor()
		if "result" in form:
			result = wiggletools.wiggleDB.query_result(cursor, form["result"].value, config['batch_system'])
			if result['status'] == "DONE":
				report_result(result)
			else:
				print json.dumps({'status':result['status']})

		elif "count" in form:
			assembly = form['assembly'].value
			params = dict((re.sub("^._", "", X), form.getlist(X)) for X in form if X != "count")
			count = len(wiggletools.wiggleDB.get_dataset_locations(cursor, params, assembly))
			print json.dumps({'query':params,'count':count})

		elif 'annotations' in form:
			assembly = form['assembly'].value
			print json.dumps({"annotations": [X[1] for X in wiggletools.wiggleDB.get_annotations(cursor, assembly)]})

		elif 'wa' in form:
			options = WiggleDBOptions()
			options.conn = conn
			options.assembly = form['assembly'].value
			options.wa = form['wa'].value
			options.working_directory = config['working_directory']
			options.s3 = config['s3_bucket']
			if 'email' in form:
				options.emails = form.getlist('email')

			if 'wb' in form:
				options.wb = form['wb'].value
			else:
				options.wb = None

			if 'w' in form:
				options.fun_merge = form['w'].value
			else:
				options.fun_merge = None

			options.a = dict((X[2:], form.getlist(X)) for X in form if X[:2] == "A_")
			options.b = dict((X[2:], form.getlist(X)) for X in form if X[:2] == "B_")
			if len(options.b.keys()) == 0:
				options.b = None

			if options.a['type'] == 'regions' and options.b['type'] == 'signal':
				tmp = options.b
				options.b = options.a
				options.a = tmp
			
			result = wiggletools.wiggleDB.request_compute(cursor, options, config['batch_system'])
			if result['status'] == 'DONE':
				report_result(result)
			else:
				print json.dumps(result)

		else:
			print json.dumps("No params, no output")

		conn.commit()
		conn.close()
        except:
                print json.dumps("ERROR")
		raise

if __name__ == "__main__":
	main()
