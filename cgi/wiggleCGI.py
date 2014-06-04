#!/usr/bin/env python

import cgi
import cgitb
import json
import sqlite3
import wiggleDB
import re

database_location = "/Users/dzerbino/Desktop/WiggleCGI/test.db"
logdir = '/Users/dzerbino/Desktop/WiggleCGI/'
root_url = 'localhost'
root_dir = r'^'

cgitb.enable(logdir=logdir)

def main():
	print "Content-Type: application/json"
	print

	try:
		form = cgi.FieldStorage()
		conn = sqlite3.connect(database_location)
		cursor = conn.cursor()
		if "result" in form:
			status, info = wiggleDB.query_result(cursor, form["result"].value)
			if status == "DONE":
				url = re.sub(root_dir, root_url, info)		
				print json.dumps({'status':status, 'url':url})
			else:
				print json.dumps({'status':status})

		elif "count" in form:
			assembly = form['assembly'].value
			params = dict((re.sub("^._", "", X), form.getlist(X)) for X in form if X != "count")
			count = len(wiggleDB.get_dataset_locations(cursor, params, assembly))
			print json.dumps({'query':params,'count':count})

		elif 'a' in form:
			options.conn = conn
			options.assembly = form['assembly'].value
			options.wa = form['wa'].value

			if 'wb' in form:
				options.wb = form['wb'].value
			else:
				options.wb = None

			if 'merge' in form:
				options.fun_merge = form['merge'].value
			else:
				options.fun_merge = None

			options.a = dict((X[2:], form[X].value) for X in form if X[:2] == "A_")
			options.b = dict((X[2:], form[X].value) for X in form if X[:2] == "B_")
			
			jobID = wiggleDB.request_compute(cursor, options)
			print json.dumps({'ID':jobID, 'options':options})

		else:
			print json.dumps("No params, no output")

		conn.commit()
		conn.close()
        except:
                print json.dumps("ERROR")

if __name__ == "__main__":
	main()
