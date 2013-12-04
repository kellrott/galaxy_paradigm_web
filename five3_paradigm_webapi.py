#!/usr/bin/env python

#############################################################################
# Copyright 2012-2013
# Kyle Ellrott <kellrott@soe.ucsc.edu>
#       and
# Five3 Genomics, LLC ( http://five3genomics.com )
# All Rights Reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
############################################################################

import sys
import json
from argparse import ArgumentParser
import uuid
import time
import requests

UPLOAD_URL = "https://dna.five3genomics.com/api/v1/paradigm/file/"
RUN_URL = "https://dna.five3genomics.com/api/v1/paradigm/run/"
RESULT_URL = "https://dna.five3genomics.com/api/v1/paradigm/run/"
API_POLL_DELAY = 600  # please don't lower this unless you have a really good reason
VERBOSE = True


def upload_file(username, api_key, filename, format, verbose=False):
    """
    Upload a file, converting filename to a uuid for uniqueness
    """
    if verbose:
        print "Uploading file {0} as type {1}".format(filename, format)
    u_name = str(uuid.uuid4())
    r = requests.post(UPLOAD_URL,
        params={'username': username, 'api_key': api_key},
        data={'file': u_name, 'format': format},
        files={'file': (u_name, open(filename, "rb"))})
    
    try:
        meta = json.loads(r.text)
    except ValueError:
        sys.stderr.write("Error Invalid response: {0}\n\tHTTPcode={1}\n".format(r.text, r.status_code))
        sys.exit(1)
        
    if "error_message" in meta:
        sys.stderr.write("Upload failed for file: {0}\n{1}\n".format(filename, meta['error_message']))
        sys.exit(1)
    return u_name


def monitor_job_uuid(username, api_key, job_uuid, verbose=False):
    """
    Monitor a run - returns False until it's finished and then returns True
    Will error out (sys.exit) if the run fails or encounters another issue
    """
    r = requests.get(RESULT_URL + job_uuid + "/",
        params={'username': username, 'api_key': api_key}
    )
    res = json.loads(r.text)
    if "error_message" in res:
        sys.stderr.write("%s\n" % (res['error_message']))
        sys.exit(1)
    elif res.get("status", "") == "Failed":
        sys.stderr.write("Run failed, please contact support@five3genomics.com for more info.\n")
        sys.exit(1)

    if verbose:
        print res
    if res.get('result', None):
        return True

    return False


def save_job_result(username, api_key, job_uuid, output_type, output_filename):
    r = requests.get(RESULT_URL + job_uuid + "/" + output_type + "/",
        params={'username': username, 'api_key': api_key},
        stream=True
    )
    fout = open(output_filename, 'w')
    for chunk in r.iter_content(chunk_size=1024): 
        if chunk: # filter out keep-alive new chunks
            fout.write(chunk)
            fout.flush()
    fout.write(r.text)
    fout.close()


def main(args):
    job_uuid = None

    if not (args.resume is not None or args.resume_file is not None):
        uploaded_names = []
        u_name = upload_file(args.user, args.api_key, args.path, 'path', VERBOSE)
        uploaded_names.append(u_name)

        ###
        # Upload the expression matrix
        ###
        if args.exp:
            u_name = upload_file(args.user, args.api_key, args.exp, 'he', VERBOSE)
            uploaded_names.append(u_name)

        ###
        # Upload the protein matrix
        ###
        if args.prot:
            u_name = upload_file(args.user, args.api_key, args.prot, 'hprot', VERBOSE)
            uploaded_names.append(u_name)

        ###
        # Upload the copy number matrix
        ###
        if args.cna:
            u_name = upload_file(args.user, args.api_key, args.cna, 'hcnv', VERBOSE)
            uploaded_names.append(u_name)

        ###
        # Upload the protein active matrix
        ###
        if args.prota:
            u_name = upload_file(args.user, args.api_key, args.prota, 'hprota', VERBOSE)
            uploaded_names.append(u_name)

        ###
        # Upload a params file
        ###
        if args.param:
            u_name = upload_file(args.user, args.api_key, args.param, 'param', VERBOSE)
            uploaded_names.append(u_name)

        ###
        # Upload a dogma file
        ###
        if args.dogma:
            u_name = upload_file(args.user, args.api_key, args.dogma, 'dogma', VERBOSE)
            uploaded_names.append(u_name)

        ###
        # Upload a imap file
        ###
        if args.imap:
            u_name = upload_file(args.user, args.api_key, args.imap, 'imap', VERBOSE)
            uploaded_names.append(u_name)

        ###
        # Assemble the run_data object and initiate the run
        ###
        run_data = {'username': args.user, 'api_key': args.api_key,
                    'name': args.name, 'files': uploaded_names,
                    'disc_0': args.disc_low, 'disc_1': args.disc_high}

        if args.posterior:
            run_data['config_top'] = 1
        if args.null_batches:
            run_data['null_batches'] = args.null_batches
        if args.skip_clustering:
            run_data['skip_clustering'] = 1
        if args.skip_em:
            run_data['skip_em'] = 1
        elif args.skip_link_em:
            run_data['config_top_em'] = 1

        r = requests.post(RUN_URL, data=run_data)
        try:
            meta = json.loads(r.text)
        except ValueError:
            sys.stderr.write("Server Message: %s" % (r.text))
            sys.exit(1)
        if "error_message" in meta:
            sys.stderr.write("Run Request failed: %s\n" % (meta['error_message']))
            sys.exit(1)
        if VERBOSE:
            print meta
        job_uuid = meta['run']['uuid']
        if args.submit_only is not None:
            job_finished = True
            handle = open(args.submit_only, "w")
            handle.write(job_uuid)
            handle.close()
        else:
            job_finished = False
        print "Running Job:", job_uuid
    else:
        if args.resume is not None:
            job_uuid = args.resume
        if args.resume_file is not None:
            handle = open(args.resume_file)
            job_uuid = handle.readline().rstrip()
            handle.close()
        print "Resuming Job:", job_uuid
        job_finished = monitor_job_uuid(args.user, args.api_key, job_uuid, VERBOSE)

    while not job_finished:
        time.sleep(API_POLL_DELAY)
        job_finished = monitor_job_uuid(args.user, args.api_key, job_uuid, VERBOSE)

    if args.submit_only is None:
        save_job_result(args.user, args.api_key, job_uuid, 'results', args.out)
        save_job_result(args.user, args.api_key, job_uuid, 'nulls', args.out_nulls)
        save_job_result(args.user, args.api_key, job_uuid, 'params', args.out_params)


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-u", "--user", dest="user", help="Username", required=True)
    parser.add_argument("-a", "--api", dest="api_key", help="API Key", required=True)
    parser.add_argument("-r", "--resume", dest="resume", help="Resume watching job UUID", default=None)
    parser.add_argument("--resume-file", dest="resume_file", help="Resume watching job UUID File", default=None)    
    parser.add_argument("-n", "--name", dest="name", help="Run Name", default='Galaxy Paradigm Run')
    parser.add_argument("--submit-only", dest="submit_only", help="Submit and return job UUID file", default=None)

    file_options = parser.add_argument_group('Files', 'Files provided to the Paradigm Run')
    file_options.add_argument("-e", "--exp", dest="exp", help="mRNA matrix", default=None)
    file_options.add_argument("-c", "--cna", dest="cna", help="Copy number matrix", default=None)
    file_options.add_argument("-p", "--path", dest="path", help="Pathway File", default=None)
    file_options.add_argument("--prot", dest="prot", help="Protein File", default=None)
    file_options.add_argument("--prota", dest="prota", help="Protein active file", default=None)
    file_options.add_argument("--param", dest="param", help="Params File", default=None)
    file_options.add_argument("--dogma", dest="dogma", help="Dogma File", default=None)
    file_options.add_argument("--imap", dest="imap", help="IMAP File", default=None)

    parser.add_argument("--null-batches", dest="null_batches", help="Number of null batches", default=None)
    parser.add_argument("--skip-em", dest="skip_em", action="store_true", help="Skip EM", default=False)
    parser.add_argument("--posterior", dest="posterior", action="store_true", help="Request all posteriors", default=False)
    parser.add_argument("--disc-low", dest="disc_low", help="Discretization lower bound", default=0.3333)
    parser.add_argument("--disc-high", dest="disc_high", help="Discretization upper bound", default=0.6667)
    parser.add_argument("--skip-link_em", dest="skip_link_em", action="store_true", help="Skip Link-Learning EM", default=False)
    parser.add_argument("--skip-clustering", dest="skip_clustering", action="store_true", help="Skip Clustering", default=False)

    parser.add_argument("-o", "--out", dest="out", help="Output filename", default="paradigm.out")
    parser.add_argument("--out-nulls", dest="out_nulls", help="Output filename for Nulls", default="paradigm.nulls.out")
    parser.add_argument("--out-params", dest="out_params", help="Output filename for Params", default="paradigm.params.out")
    args = parser.parse_args()

    if args.user is None or args.api_key is None:
        sys.stderr.write("Need Web API login/key\n")
        sys.exit(1)

    if not (args.resume or args.resume_file) and args.path is None:
        sys.stderr.write("Need to provide pathway file\n")
        sys.exit(1)

    main(args)
