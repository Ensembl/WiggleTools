#!/usr/bin/env python

import cgi
import cgitb
import json
import sqlite3
import wiggleDB

cgitb.enable(logdir='.')

def main():
	form=cgi.FieldStorage()
	conn = sqlite3.connect("wiggletools.db")
	cursor = conn.cursor()

	if "result" in form:
		status, info = wiggleDB.query_result(cursor, form["result"])
		if status == "DONE":
			print json.dumps({'status':status, 'url':info})
		else:
			print json.dumps({'status':status})
	else:
		jobID = wiggleDB.request_compute(cursor, form)
		print json.dumps({'ID':jobID})

	conn.commit()
	conn.close()

if __name__ == "__main__":
	main()
