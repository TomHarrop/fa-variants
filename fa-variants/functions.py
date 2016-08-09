#!/usr/bin/python3
# -*- coding: utf-8 -*-

import subprocess
import re
import os
import datetime

#############
# UTILITIES #
#############


def flatten_list(l):
    for x in l:
        if hasattr(x, '__iter__') and not isinstance(x, str):
            for y in flatten_list(x):
                yield y
        else:
            yield x


# touch function for updating ruffus flag files
def touch(fname, mode=0o666, dir_fd=None, **kwargs):
    flags = os.O_CREAT | os.O_APPEND
    with os.fdopen(os.open(fname, flags=flags, mode=mode, dir_fd=dir_fd)) as f:
        os.utime(f.fileno() if os.utime in os.supports_fd else fname,
                 dir_fd=None if os.supports_fd else dir_fd, **kwargs)

############################
# JOB SUBMISSION FUNCTIONS #
############################


def submit_job(job_script, ntasks, cpus_per_task, job_name, extras=[]):
    # type: (str, str, str, str, str, list) -> str
    '''
    Submit the job using salloc hack. When complete return job id and write
    output to file.
    '''
    # call salloc as subprocess
    proc = subprocess.Popen(['salloc', '--ntasks=' + ntasks,
                             '--cpus-per-task=' + cpus_per_task,
                             '--job-name=' + job_name, job_script] +
                            list(extras),
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    # get stdout and stderr
    out, err = proc.communicate()
    # parse stderr (salloc output) for job id
    job_regex = re.compile(b'\d+')
    job_id_bytes = job_regex.search(err).group(0)
    job_id = job_id_bytes.decode("utf-8")
    # write stderr & stdout to log file
    out_file = 'ruffus/' + job_name + '.' + job_id + '.ruffus.out.txt'
    with open(out_file, 'wb') as f:
        f.write(out)
    err_file = 'ruffus/' + job_name + '.' + job_id + '.ruffus.err.txt'
    with open(err_file, 'wb') as f:
        f.write(err)
    # mail output
    if proc.returncode != 0:
        subject = "[Tom@SLURM] Pipeline step " + job_name + " FAILED"
    else:
        subject = "[Tom@SLURM] Pipeline step " + job_name + " finished"
    mail = subprocess.Popen(['mail', '-s', subject, '-A', out_file, '-A',
                             err_file, 'tom'], stdin=subprocess.PIPE)
    mail.communicate()
    # check subprocess exit code
    os.remove(out_file)
    os.remove(err_file)
    assert proc.returncode == 0, ("Job " + job_name +
                                  " failed with non-zero exit code")
    return(job_id)


def print_job_submission(job_name, job_id):
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
    print('[', now, '] : Job ' + job_name + ' run with JobID ' + job_id)

######################
# FUNCTION GENERATOR #
######################


def generate_job_function(
        job_script, job_name, job_type='transform', ntasks=1,
        cpus_per_task=1, extras=False, verbose=False):

    '''Generate a function for a pipeline job step'''

    # job_type determines the number of arguments accepted by the returned
    # function. Ruffus 'transform' and 'merge' functions should  expect two
    # positional arguments, because Ruffus will pass input_files  and
    # output_files positionally. 'originate' functions should expect
    # output_files as a positional argument and no input_files. Additional
    # arguments for the 'extra' parameter will be passed as a list from
    # Ruffus.

    # check job_type
    _allowed_job_types = ['transform', 'merge', 'originate', 'download']
    if job_type not in _allowed_job_types:
        raise ValueError('{job_type} not an allowed job_type')

    # set up the args
    function_args = []
    # if we expect input_files, they go first
    if job_type in ['transform', 'merge']:
        function_args.append('input_files')
    # all job_types have output_files
    function_args.append('output_files')
    # download jobs have logon details
    if job_type == 'download':
        function_args.append('jgi_logon')
        function_args.append('jgi_password')
    # extras go at the end
    if extras:
        function_args.append('extras')

    # define the function
    def job_function(*function_args):

        # standardise submit_args and handle them in bash.
        # provide arguments to submit_job extras argument in the following
        # order:
        # -i: input_files
        # -o: output_files
        # -e: jgi_logon (email)
        # -p: jgi_password
        # extra arguments passed verbatim from Ruffus

        function_args_list = list(function_args)
        submit_args = []

        if verbose:
            print("\nfunction_args: ", function_args)
            print("\nfunction_args_list: ", function_args_list)

        # if we expect input_files, they go first
        if job_type in ['transform', 'merge']:
            input_files = [function_args_list.pop(0)]
            input_files_flat = list(flatten_list(input_files))
            y = ['-i'] * len(input_files_flat)
            new_args = [x for t in
                        zip(y, input_files_flat)
                        for x in t]
            submit_args.append(new_args)

            if verbose:
                print("\ninput_files_flat:", input_files_flat)
                print("\nnew_args: ", new_args)
                print("\nsubmit_args: ", submit_args)

        # output_files required for all job_types
        output_files = [function_args_list.pop(0)]
        output_files_flat = list(flatten_list(output_files))
        y = ['-o'] * len(output_files_flat)
        new_args = [x for t in
                    zip(y, output_files_flat)
                    for x in t]
        submit_args.append(new_args)

        if verbose:
            print("\noutput_files_flat:", output_files_flat)
            print("\nnew_args: ", new_args)
            print("\nsubmit_args: ", submit_args)

        # if we have logon details they go here
        if job_type == 'download':
            submit_args.append('-e')
            submit_args.append(function_args_list.pop(0))
            submit_args.append('-p')
            submit_args.append(function_args_list.pop(0))

            if verbose:
                print("\nsubmit_args: ", submit_args)

        # extras go at the end
        if extras:
            submit_args.append(function_args_list.pop(0))

            if verbose:
                print("\nsubmit_args: ", submit_args)

        # did we use everything?
        if verbose:
            print("\nsubmit_args: ", submit_args)
            print("end_fal: ", function_args_list)

        if len(function_args_list) > 0:
            raise ValueError('unused function_args_list')

        # flatten the list
        submit_args_flat = list(flatten_list(submit_args))

        if verbose:
            print("\nsubmit_args_flat: ", submit_args_flat)

        # submit the job. n.b. the job script has to handle the extras
        # properly, probably by parsing ${1}.. ${n} in bash.
        job_id = submit_job(
            job_script=job_script,
            ntasks=str(ntasks),
            cpus_per_task=str(cpus_per_task),
            job_name=job_name,
            extras=list(submit_args_flat))
        print_job_submission(job_name, job_id)

    return job_function
