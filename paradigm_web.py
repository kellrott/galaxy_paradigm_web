#!/usr/bin/env python

import 	os
import sys
import urllib
import json
from argparse import ArgumentParser
import uuid
import time
import requests

UPLOAD_URL = "https://dna.five3genomics.com/api/v1/legacy_upload_file/"
RUN_URL = "https://dna.five3genomics.com/api/v1/paradigm_legacy_run/"
RESULT_URL = "https://dna.five3genomics.com/api/v1/paradigm_legacy_run/"

if __name__ == "__main__":
	parser = ArgumentParser()
	parser.add_argument("-user", dest="user", help="Username", default=None)
	parser.add_argument("-api", dest="apikey", help="API Key", default=None)
	
	parser.add_argument("-exp", dest="exp", help="mRNA matrix", default=None)
	parser.add_argument("-cna", dest="cna", help="Copy number matrix", default=None)
	parser.add_argument("-path", dest="path", help="Pathway File", default=None)

	parser.add_argument("-null_batches", dest="null_batches", help="Null batches", default=None)
	parser.add_argument("-skip_em", dest="skip_em", action="store_true", help="Skip EM", default=False)
	parser.add_argument("-disc", dest="disc", help="Discretization", default=False)

	parser.add_argument("-out", dest="out", help="Output Path", default="paradigm.out")
	parser.add_argument("-resume", dest="resume", help="Resume watching job UUID", default=None)
	args = parser.parse_args()
	
	if args.user is None or args.apikey is None:
		sys.stderr.write("Need Web API login/key\n")
		sys.exit(1)
	
	job_uuid = None
	
	if args.resume is None:
		uploaded_names = []
		
		###
		#Upload the pathway file
		###
		if args.path is None:
			sys.stderr.write("Need to provide pathway file\n")
			sys.exit(1)

		u_name = str(uuid.uuid4())
		r = requests.post(UPLOAD_URL, 
			data={ 'username' : args.user, 'api_key' : args.apikey, 'file' : u_name, 'format' : 'path' }, 
			files={ 'file' : (u_name, open(args.path, "rb"))  })
		meta = json.loads(r.text)
		if "error_message" in meta:
			sys.stderr.write("%s\n" % (meta['error_message']))
			sys.exit(1)
		uploaded_names.append(u_name)

		###
		#Upload the expression matrix
		###		
		if args.exp is not None:
			u_name = str(uuid.uuid4())
			r = requests.post(UPLOAD_URL, 
				data={ 'username' : args.user, 'api_key' : args.apikey, 'file' : u_name, 'format' : 'he' }, 
				files={ 'file' : (u_name, open(args.exp, "rb") ) })
			meta = json.loads(r.text)
			if "error_message" in meta:
				sys.stderr.write("%s\n" % (meta['error_message']))
				sys.exit(1)
			uploaded_names.append(u_name)

		
		###
		#Upload the copy number matrix
		###		
		if args.exp is not None:
			u_name = str(uuid.uuid4())
			r = requests.post(UPLOAD_URL, 
				data={ 'username' : args.user, 'api_key' : args.apikey, 'file' : u_name, 'format' : 'hcnv' }, 
				files={ 'file' : (u_name, open(args.cna, "rb"))  })
			meta = json.loads(r.text)
			if "error_message" in meta:
				sys.stderr.write("%s\n" % (meta['error_message']))
				sys.exit(1)
			uploaded_names.append(u_name)
		
		run_data = {'username' : args.user, 'api_key' : args.apikey, 'name' : 'Galaxy Paradigm Run', 'files' : uploaded_names }
		
		if args.null_batches is not None:
			run_data['null_batches'] = args.null_batches
		if args.skip_em:
			run_data['skip_em'] = 1
		if args.disc is not None:
			run_data['disc'] = args.disc
		
		r = requests.post(RUN_URL,
			data=run_data
		)
		meta = json.loads(r.text)
		if "error_message" in meta:
			sys.stderr.write("%s\n" % (meta['error_message']))
			sys.exit(1)
		print meta
		job_uuid = meta['paradigm_legacy_run']['uuid']
	else:
		job_uuid = args.resume
		
	if job_uuid is not None:
		while True:
			r = requests.get(RESULT_URL + job_uuid + "/", 
				params={'username' : args.user, 'api_key' : args.apikey }
			)
			res = json.loads(r.text)
			if "error_message" in res:
				sys.stderr.write("%s\n" % (meta['error_message']))
				sys.exit(1)
			if 'result' not in res: #debugging something I observed once
				print res
			if res['result'] is not None: break
			print res
			time.sleep(30)
		
		sample_set = {}
		probe_set = {}
		for sample, probe, value in res['result']:
			sample_set[sample] = True
			if probe not in probe_set:
				probe_set[probe] = [ [sample, value] ]
			else:
				probe_set[probe].append( [sample, value] )
		probe_list = sorted(probe_set.keys())
		sample_list = sorted(sample_set.keys())
		sample_order = {}
		for x, y in enumerate(sample_list):
			sample_order[y] = x
		handle = open(args.out, "w")
		
		handle.write("probe\t%s\n" % ("\t".join(sample_list)))
		for cur_probe in probe_list:
			v = [float("nan")] * len(sample_list)
			for sample, value in probe_set[cur_probe]:
				v[ sample_order[sample] ] = value
			handle.write("%s\t%s\n" % (cur_probe.encode('ascii', 'ignore'), "\t".join( map(lambda x : "%f" % x, v) ) ))
		handle.close()
		