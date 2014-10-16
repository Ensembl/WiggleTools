#!/usr/bin/env python

import sys
import re
import argparse
import sqlite3
import subprocess
import tempfile
import os
import os.path
import json

import parallelWiggleTools
import multiJob

verbose = False

###########################################
## Configuration file
###########################################

def read_config_file(filename):
	return dict(line.strip().split('\t') for line in open(filename) if line[0] != '#' and len(line) > 1)


###########################################
## Command line interface
###########################################

def normalise_spaces(string):
	if string is None:
		return None
	else:
		return re.sub("^\W+", "", re.sub("\W+", " ", string))

def get_options():
	parser = argparse.ArgumentParser(description='WiggleDB backend.')
	parser.add_argument('--db', '-d', dest='db', help='Database file',required=True)
	parser.add_argument('-a',dest='a',help='A set of SQL constraints',nargs='*')
	parser.add_argument('-wa',dest='wa',help='WiggleTools command for A')
	parser.add_argument('-b',dest='b',help='A second set of SQL constraints',nargs='*')
	parser.add_argument('-wb',dest='wb',help='WiggleTools command for B')
	parser.add_argument('--wiggletools','-w',dest='fun_merge',help='Wiggletools command')
	parser.add_argument('--load','-l',dest='load',help='Datasets to load in database')
	parser.add_argument('--load_assembly','-la',dest='load_assembly',help='Assembly name and path to file with chromosome lengths',nargs=2)
	parser.add_argument('--assembly','-y',dest='assembly',help='File with chromosome lengths')
	parser.add_argument('--clean',dest='clean',help='Delete cached datasets older than X days', type=int)
	parser.add_argument('--dump_cache',dest='dump_cache',help='Dump cache info', action='store_true')
	parser.add_argument('--datasets',dest='datasets',help='Print dataset info', action='store_true')
	parser.add_argument('--clear_cache',dest='clear_cache',help='Reset cache info', action='store_true')
	parser.add_argument('--remember',dest='remember',help='Preserve dataset from garbage collection', action='store_true')
	parser.add_argument('--dry-run',dest='dry_run',help='Do not run the command, print wiggletools command', action='store_true')
	parser.add_argument('--result','-r',dest='result',help='Return status or end result of job', type=int)
	parser.add_argument('--attributes','-t',dest='attributes',help='Print JSON hash of attributes and values', action='store_true')
	parser.add_argument('--verbose','-v',dest='verbose',help='Turn on status output',action='store_true')
	parser.add_argument('--config','-c',dest='config',help='Configuration file')
	parser.add_argument('--annotations','-n',dest='annotations',help='Print list of annotation names', action='store_true')
	parser.add_argument('--jobs','-j',dest='jobs',help='Print list of jobs', action='store_true')

	options = parser.parse_args()
	if options.load is not None:
		assert not os.path.exists(options.db), "Cannot overwrite pre-existing database %s" % options.db
	if all(X is None for X in [options.load, options.clean, options.result, options.load_assembly, options.datasets]) and not options.dump_cache and  not options.clear_cache and not options.attributes and not options.annotations:
		assert options.a is not None, 'No dataset selection to run on'
		assert options.wa is not None, 'No dataset transformation to run on'
		assert options.assembly is not None, 'No assembly name specified'
		if options.b is not None:
			assert options.fun_merge is not None, 'No action command (load,clean,compute) specified'

	options.wa = normalise_spaces(options.wa)	
	options.wb = normalise_spaces(options.wb)	
	options.fun_merge = normalise_spaces(options.fun_merge)

	if options.config is not None:
		config = read_config_file(options.config)
	else:
		config = None

	global verbose
	verbose = options.verbose
	return options, config

###########################################
## Creating a database
###########################################

def create_database(cursor, filename):
	if verbose:
		print 'Creating database'
	create_assembly_table(cursor)
	create_cache(cursor)
	create_dataset_table(cursor, filename)
	create_job_table(cursor)

def create_assembly_table(cursor):
	cursor.execute('''
	CREATE TABLE
	assemblies
	(
	name varchar(255),
	location varchar(1000)
	)
	''')

def create_job_table(cursor):
	cursor.execute('''
	CREATE TABLE
	jobs
	(
	job_id INTEGER PRIMARY KEY AUTOINCREMENT,
	lsf_id int,
	lsf_id2 int,
	temp varchar(1000),
	status varchar(255),
	email varchar(255)
	)
	''')

def create_cache(cursor):
	cursor.execute('''
	CREATE TABLE
	cache
	(
	job_id int,
	query varchar(10000),
	location varchar(1000),
	remember bit,
	primary_loc bit,
	last_query datetime
	)
	''')

def create_dataset_table(cursor, filename):
	file = open(filename)
	items = file.readline().strip().split('\t')
	assert items[:5] == list(('location','name','type','annotation','assembly')), "Badly formed dataset table, please ensure the first five columns refer to location, name, type, annotation and assembly"
	header = '''
			CREATE TABLE IF NOT EXISTS 
			datasets 
			(
			location varchar(1000),
 			name varchar(100), 
			type varchar(100), 
			annotation bit, 
			assembly varchar(100),
		 '''
	cursor.execute('\n'.join([header] + [",\n".join(['%s varchar(255)' % X for X in items[5:]])] + [')']))

	for line in file:
		cursor.execute('INSERT INTO datasets VALUES (%s)' % ",".join("'%s'" % X for X in line.strip().split('\t')))
	file.close()

###########################################
## Loading assembly info
###########################################

def load_assembly(cursor, assembly_name, chrom_sizes):
	if verbose:
		print 'Loading path to assembly chromosome length %s for %s' % (chrom_sizes, assembly_name)
	cursor.execute('INSERT INTO assemblies VALUES(\'%s\',\'%s\')' % (assembly_name, chrom_sizes))

###########################################
## Garbage cleaning 
###########################################

def clean_database(cursor, days):
	for location in cursor.execute('SELECT location FROM cache WHERE julianday(\'now\') - julianday(last_query) > %i AND remember = 0' % days).fetchall():
		if verbose:
			print 'Removing %s' % location[0]
		os.remove(location[0])
	cursor.execute('DELETE FROM cache WHERE julianday(\'now\') - julianday(last_query) > %i AND remember = 0' % days)

	for temp in cursor.execute('SELECT temp FROM jobs WHERE status="DONE" OR status="EMPTY"').fetchall():
		if verbose:
			print 'Removing %s and derived files' % temp[0]
		multiJob.clean_temp_file(temp[0])

###########################################
## Search datasets
###########################################

def get_dataset_attributes_2(cursor):
	return [X[1] for X in cursor.execute('PRAGMA table_info(datasets)').fetchall()]

def get_dataset_attributes(cursor):
	return list(set(get_dataset_attributes_2(cursor)) - set(["annotation","name","assembly","location","type"]))

def get_attribute_values_2(cursor, attribute):
	return [X[0] for X in cursor.execute('SELECT DISTINCT %s FROM datasets' % (attribute)).fetchall()]

def get_attribute_values(cursor):
	return dict((attribute, get_attribute_values_2(cursor, attribute)) for attribute in get_dataset_attributes(cursor))

def get_annotations(cursor, assembly):
	return cursor.execute('SELECT * FROM datasets WHERE assembly=? AND annotation', (assembly,)).fetchall()

def get_jobs(cursor):
	res = cursor.execute('SELECT * FROM jobs').fetchall()
	desc = cursor.description
	return [[X[0] for X in desc]] + res

def get_datasets(cursor):
	res = cursor.execute('SELECT * FROM datasets').fetchall()
	desc = cursor.description
	return [[X[0] for X in desc]] + res

def attribute_selector(attribute, params):
	return "( %s )" % " OR ".join("%s=:%s_%i" % (attribute,attribute,index) for index in range(len(params[attribute])))

def denormalize_params(params):
	return dict(("%s_%i" % (attribute, index),value) for attribute in params for (index, value) in enumerate(params[attribute]))

def get_dataset_locations(cursor, params, assembly):
	# Quick check that all the keys are purely alphanumeric to avoid MySQL injections
	assert not any(re.match('\W', X) is not None for X in params)
	params['assembly'] = [assembly]
	query = " AND ".join(attribute_selector(X, params) for X in params)
	res = cursor.execute('SELECT location FROM datasets WHERE ' + query, denormalize_params(params)).fetchall()
	if verbose:
		query_txt = " AND ".join("%s=%s" % (X,params[X]) for X in params)
		print 'Query: SELECT location FROM datasets WHERE ' + query_txt
		print 'Found:\n' + "\n".join(X[0] for X in res)
	return sorted(X[0] for X in res)

###########################################
## Search cache
###########################################

def reset_time_stamp(cursor, cmd):
	cursor.execute('UPDATE cache SET last_query= date(\'now\') WHERE query = \'%s\'' % cmd)

def get_precomputed_jobID(cursor, cmd):
	reset_time_stamp(cursor, cmd)
	reports = cursor.execute('SELECT job_id FROM cache WHERE query = \'%s\'' % cmd).fetchall()
	if len(reports) == 0:
		if verbose:
			print 'Did not find prior job for query: %s' % cmd
		return None
	else:
		if verbose:
			print 'Found prior job for query: %s' % cmd
			print reports[0][0]
		return reports[0][0]

def get_job_location_2(cursor, jobID):
	reports = cursor.execute('SELECT location FROM cache WHERE job_id = \'%s\' and primary_loc=1' % jobID).fetchall()
	assert len(reports) == 1
	return reports[0][0]

def get_job_location(db, jobID):
	connection = sqlite3.connect(db)
	cursor = connection.cursor()
	res = get_job_location_2(cursor, jobID)
	connection.close()
	return res

def get_precomputed_location(cursor, cmd):
	reports = cursor.execute('SELECT location FROM jobs NATURAL JOIN cache WHERE status="DONE" OR status="EMPTY" AND query = \'%s\'' % cmd).fetchall()
	if len(reports) > 0:
		reset_time_stamp(cursor, cmd)
		if verbose:
			print 'Found pre-computed file for query: %s' % cmd
			print reports[0]
		return reports[0][0]
	else:
		if verbose:
			print 'Did not find pre-computed file for query: %s' % cmd
		return None

def reuse_or_write_precomputed_location(cursor, cmd, working_directory):
	pre_location = get_precomputed_location(cursor, cmd)
	if pre_location is not None:
		return pre_location, pre_location, False
	else:
		fh, destination = tempfile.mkstemp(suffix='.bw',dir=working_directory)
		return 'write %s %s' % (destination, cmd), destination, True

def launch_compute(cursor, fun_merge, fun_A, data_A, fun_B, data_B, options, normalised_form, batch_system):
	destination = None
	destinationA = None
	destinationB = None
	cmds = None
	histogram = None
	apply_paste = None

	cmd_A = " ".join([fun_A] + data_A + [':'])
	cmd_A2, destinationA, computeA = reuse_or_write_precomputed_location(cursor, cmd_A, options.working_directory)

	if data_B is not None:
		merge_words = fun_merge.split(' ')

		assert fun_merge is not None
		if fun_B is not None:
			cmd_B = " ".join([fun_B] + data_B + [':'])
			cmd_B2, destinationB, computeB = reuse_or_write_precomputed_location(cursor, cmd_B, options.working_directory)
		else:
			cmd_B2 = " ".join(data_B)
			computeB = False

		if merge_words[0] == 'histogram':
			cmds = []
			if computeA:
				cmds = [cmd_A2]
			if computeB:
				cmds.append(cmd_B2)
			width = merge_words[1]

			fh, destination = tempfile.mkstemp(suffix='.txt',dir=options.working_directory)
			if fun_B is not None:
		        	histogram = "histogram %s %s %s mult %s %s" % (destination, width, destinationA, destinationA, destinationB)
			elif data_B is not None:
				histogram = "histogram %s %s %s" % (destination, width, " ".join("mult %s %s" % (destinationA, X) for X in data_B))
		elif merge_words[0] == 'profile':
			fh, destination = tempfile.mkstemp(suffix='.txt',dir=options.working_directory)
			cmds = [" ".join(['profile', destination, merge_words[1], cmd_B2, cmd_A2])]
		elif merge_words[0] == 'profiles':
			fh, destination = tempfile.mkstemp(suffix='.txt',dir=options.working_directory)
			cmds = [" ".join(['profiles', destination, merge_words[1], cmd_B2, cmd_A2])]
		elif merge_words[0] == 'apply_paste':
			cmds = []
			if computeA:
				cmds = [cmd_A2]
			if computeB:
				cmds.append(cmd_B2)
			fh, destination = tempfile.mkstemp(suffix='.txt',dir=options.working_directory)
			assert len(data_B) == 1, "Cannot apply_paste to multiple files %s\n" % " ".join(data_B)
			apply_paste = " ".join(['apply_paste', destination, 'AUC', data_B[0], destinationA])
		else:
			fh, destination = tempfile.mkstemp(suffix='.bw',dir=options.working_directory)
			cmds = [" ".join(['write', destination, fun_merge, cmd_A2, cmd_B2])]
	else:
		computeB = False
		cmds = [cmd_A2]
		destination = destinationA
		destinationA = None

	if verbose:
		print "PARALLEL CMDS : " + "\n".join(cmds)
		if histogram is not None:
			print "HISTOGRAM " + histogram
		if apply_paste is not None:
			print "APPLY PASTE : " + apply_paste

	chrom_sizes = get_chrom_sizes(cursor, options.assembly)
	if len(cmds) > 0 and not options.dry_run:
		lsfID, files = parallelWiggleTools.run(cmds, chrom_sizes, batch_system=batch_system, tmp=options.working_directory)
		cursor.execute('INSERT INTO jobs (lsf_id, status) VALUES (%i, "LAUNCHED")' % int(lsfID))
		jobID = cursor.execute('SELECT LAST_INSERT_ROWID()').fetchall()[0][0]
		cursor.execute('INSERT INTO cache (job_id,primary_loc,query,remember,last_query,location) VALUES (\'%s\',1,\'%s\',\'%i\',date(\'now\'),\'%s\')' % (jobID, normalised_form, int(options.remember), destination))
		if computeA: 
			cursor.execute('INSERT INTO cache (job_id,primary_loc,query,remember,last_query,location) VALUES (\'%s\',0,\'%s\',\'0\',date(\'now\'),\'%s\')' % (jobID, cmd_A, destinationA))
		if computeB:
			cursor.execute('INSERT INTO cache (job_id,primary_loc,query,remember,last_query,location) VALUES (\'%s\',0,\'%s\',\'0\',date(\'now\'),\'%s\')' % (jobID, cmd_B, destinationB))
		options.conn.commit()
	else:
		lsfID = None
		files = ['NO FILES']
		jobID = -1

	finishCmd = 'wiggleDB_finish.py --db %s --jobID %s --temp %s --config %s' % (options.db, jobID, " ".join(files), options.config)
	if histogram is not None:
		if fun_B is not None:
			finishCmd += ' --histogram \'%s\' --labels Overall Regions' % (histogram)
		elif data_B is not None:
			finishCmd += ' --histogram \'%s\' --labels %s' % (histogram, " ".join(".".join(os.path.basename(X).split(".")[:-1])))
		else:
			finishCmd += ' --histogram \'%s\' --labels Overall' % (histogram)
	elif apply_paste is not None:
		finishCmd += ' --apply_paste \'%s\'' % (apply_paste)

	if options.dry_run:
		return ";".join(cmds + [finishCmd])
	else:
		lsfID2, temp = multiJob.submit([finishCmd], batch_system=batch_system, dependency=lsfID, working_directory=options.working_directory)
		cursor.execute('UPDATE jobs SET lsf_id2=\'%s\',temp=\'%s\' WHERE job_id=\'%s\'' % (lsfID2, temp, jobID))
		return jobID

def get_chrom_sizes(cursor, assembly):
	res = cursor.execute('SELECT location FROM assemblies WHERE name = \'%s\'' % (assembly)).fetchall()
	assert len(res) == 1, 'Could not find any match to the query %s' % query
	return res[0][0]

def make_normalised_form(fun_merge, fun_A, data_A, fun_B, data_B):
	cmd_A = " ".join([fun_A] + data_A)
	if data_B is not None:
		if fun_B is not None:
			cmd_B = " ".join([fun_B] + data_B)
		else:
			cmd_B = " ".join(data_B)
		res = "; ".join([fun_merge, cmd_A, cmd_B])
	else:
		res = cmd_A

	if verbose:
		print 'CMD A: ' + cmd_A
		print 'CMD B: ' + str(cmd_B)
		print 'CMD  : ' + res

	return res

def request_compute(cursor, options, batch_system):
	fun_A = options.wa 
	data_A = get_dataset_locations(cursor, options.a, options.assembly)
	assert len(data_A) > 0
	cmd_A = " ".join([fun_A] + data_A)

	if options.b is not None:
		fun_B = options.wb
		data_B = get_dataset_locations(cursor, options.b, options.assembly)
		assert len(data_B) > 0
	else:
		data_B = None
		fun_B = None
		cmd_B = None

	normalised_form = make_normalised_form(options.fun_merge, fun_A, data_A, fun_B, data_B)
	prior_jobID = get_precomputed_jobID(cursor, normalised_form)
	if prior_jobID is not None:
		return query_result(cursor, prior_jobID, batch_system)
	else:
		return {'ID':launch_compute(cursor, options.fun_merge, fun_A, data_A, fun_B, data_B, options, normalised_form, batch_system), 'status':'LAUNCHED'}

def query_result(cursor, jobID, batch_system):
	reports = cursor.execute('SELECT status, lsf_id2 FROM jobs WHERE job_id =?', (jobID,)).fetchall()

	if len(reports) == 0:
		return {'ID':jobID, 'status':'UNKNOWN'}
	else:
		assert len(reports) == 1, 'Found %i status reports for job %s' % (len(reports), jobID)

	status, lsfID = reports[0]
	if status == 'DONE':
		return {'ID':jobID, 'status':'DONE', 'location':get_job_location_2(cursor, jobID)}
	elif status == 'EMPTY':
		return {'ID':jobID, 'status':'EMPTY'}
	elif batch_system == 'LSF':
		p = subprocess.Popen(['bjobs','-noheader',str(lsfID)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		ret = p.wait()
		(stdout, stderr) = p.communicate()
		assert ret == 0, 'Error when polling LSF job %i' % lsfID
		values = []
		for line in stdout.split('\n'):
			items = re.split('\W*', line)
			if len(items) > 2:
				values.append(items[2])
		return {'ID':jobID, 'status':"WAITING", 'return_values':values}
	elif batch_system == 'SGE':
		p = subprocess.Popen(['qstat','-j',str(lsfID)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		ret = p.wait()
		(stdout, stderr) = p.communicate()
		if ret == 0:
			count = 0
			for line in stdout.split('\n'):
				items = re.split('\W*', line)
				if items[0] == 'usage':
					count += 1
			values = '%i RUNNING' % (count)
			return {'ID':jobID, 'status':"WAITING", 'return_values':values}
		else:
			p = subprocess.Popen(['qacct','-j',str(lsfID)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			p.wait()
			(stdout, stderr) = p.communicate()
			assert p.returncode == 0, 'Error when polling SGE job %i' % lsfID
			values = []
			failedTask = False
			for line in stdout.split('\n'):
				items = re.split('\W*', line)
				if items[0] == 'failed' and items[1] != '0':
					values.append(" ".join(items[1:]))
					failedTask = True
				elif items[0] == 'exit_status': 
					if failedTask:
						failedTaskID = False
					else:
						values.append(items[1])
			if any(X != '0' for X in values):
				return {'ID':jobID, 'status':"ERROR", 'return_values':values, 'LSF_ID':lsfID}
			else:
				return {'ID':jobID, 'status':"WAITING", 'return_values':values, 'LSF_ID':lsfID}
			  
	else:
		raise NameError
		return {'ID':jobID, 'status':'CONFIG_ERROR'}

###########################################
## When a job finishes:
###########################################

def mark_job_as_done(db, jobID, empty=False):
	conn = sqlite3.connect(db)
	cursor = conn.cursor()
	if empty:
		cursor.execute('UPDATE jobs SET status = \'EMPTY\' WHERE job_id = \'%s\'' % jobID)
	else:
		cursor.execute('UPDATE jobs SET status = \'DONE\' WHERE job_id = \'%s\'' % jobID)
	conn.commit()
	conn.close()

def send_SMTP(msg, email, config):
	import smtplib
	s = smtplib.SMTP(config['smtp_server'], config['smtp_port'])
        s.ehlo()
        s.starttls()
        s.ehlo()
        s.login(config['user'], config['password'])
	s.sendmail(config['reply_to'], [email], msg.as_string())
	s.quit()

def send_email(text, email, config):
	from email.mime.text import MIMEText
	msg = MIMEText(text)
	msg['Subject'] = '[WiggleTools] Job succeeded'
	msg['From'] = config['reply_to']
	msg['To'] = email
        msg['sendername'] = config['sendername']
	send_SMTP(msg, email, config)

def visible_url(location, config):
	if 's3_bucket' in config:
		base_url = 'http://s3-%s.amazonaws.com/%s/' % (config['s3_region'], config['s3_bucket'])
		return re.sub(config['working_directory'], base_url, location)		
	else:
		return location

def report_to_user(email, jobID, data, config):
	if email is None:
		return
	else:
		url = visible_url(data, config)
		text = "Hello\n\n"
		text += "Your job %i is now finished, please refer to the WiggleTools server for your results\n\n" % jobID
		text += "You can download all the results here:\n\n"
		text += url + "\n\n"
		if url[-3] == '.bw' or url[:-3] == ".bb":
			text += "Or you can view them directly on Ensembl:\n\n"
			text += 'http://%s/%s/Location/View?g=%s;contigviewbottom=url:%s' % (config['ensembl_server'], config['ensembl_species'], config['ensembl_gene'], url)
		else:
			text += "Or you can view the graphic:\n\n"
			text += url + ".pdf"
		text += "\nBest regards,\n\n"
		text += "The WiggleTools team"
		send_email(text, email, config)

def report_empty_to_user(email, jobID, config):
	if email is None:
		return
	else:
		url = visible_url(data, config)
		text = "Hello\n\n"
		text += "Your job %i is now finished, but returned empty results\n\n" % jobID
		text += "\nBest regards,\n\n"
		text += "The WiggleTools team"
		send_email(text, email, config)

###########################################
## Main
###########################################

def main():
	options, config = get_options()
	options.conn = sqlite3.connect(options.db)
	cursor = options.conn.cursor()

	if config is None or 'batch_system' not in config:
		batch_system = 'SGE'
	else:
		batch_system = config['batch_system']

	if options.load is not None:
		create_database(cursor, options.load)
	elif options.load_assembly is not None:
		load_assembly(cursor, options.load_assembly[0], options.load_assembly[1])
	elif options.clean is not None:
		clean_database(cursor, options.clean)
	elif options.result is not None:
		print json.dumps(query_result(cursor, options.result, batch_system))
	elif options.dump_cache:
		for entry in cursor.execute('SELECT * FROM cache').fetchall():
			print entry
	elif options.clear_cache:
		cursor.execute('DROP TABLE cache')
		create_cache(cursor)
		cursor.execute('DROP TABLE jobs')
		create_job_table(cursor)
	elif options.attributes:
		print json.dumps(get_attribute_values(cursor))
	elif options.jobs:
		print "\n".join("\t".join(map(str, X)) for X in get_jobs(cursor))
	elif options.datasets:
		print "\n".join("\t".join(map(str, X)) for X in get_datasets(cursor))
	elif options.annotations:
		print "\n".join("\t".join(map(str, X)) for X in get_annotations(cursor, options.assembly))
	else:
		if options.a is not None:
			options.a = dict((X[0],X[1]) for X in (Y.split('=') for Y in options.a))
		if options.b is not None:
			options.b = dict((X[0],X[1]) for X in (Y.split('=') for Y in options.b))
		print json.dumps(request_compute(cursor, options, batch_system))

	options.conn.commit()
	options.conn.close()

if __name__=='__main__':
	main()
