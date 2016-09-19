#!/usr/bin/python3
# -*- coding: utf-8 -*-

import subprocess
import re
import os
import datetime
import tempfile

#############
# UTILITIES #
#############


# flatten list cf. `unlist()` in R
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


# convert input and output files to a bash-friendly option
def io_file_to_bash_flag(file_name, file_type, debug=False):
    # dictionary of options:
    # -b: input_bam
    # -c: output_bam
    # -d: output_bai
    # -e: jgi_logon
    # -f: input_fa
    # -g: output_fa
    # -h: output_dict
    # -i: output_fai
    # -j: input_gtf
    # -k: output_gtf
    # -l: input_bed
    # -m: output_bed
    # -p: jgi_password
    # -r: output_pdf
    # -t: input_table
    # -u: output_table
    # -v: input_vcf
    # -w: output_vcf
    # -y: other_input
    # -z: other_output

    # hard code input/output type
    input_flags = {
        '.bam': 'b',
        '.fa': 'f',
        '.fasta': 'f',
        '.gtf': 'j',
        '.bed': 'l',
        '.table': 't',
        '.vcf': 'v',
        '.Rds': 'y'}
    output_flags = {
        '.bam': 'c',
        '.bai': 'd',
        '.fa': 'g',
        '.fasta': 'g',
        '.dict': 'h',
        '.fai': 'i',
        '.gtf': 'k',
        '.bed': 'm',
        '.pdf': 'r',
        '.table': 'u',
        '.vcf': 'w',
        '.Rds': 'z'}

    # get the first extension and deal with .gz files
    file_ext = os.path.splitext(file_name)[1]
    if file_ext == '.gz':
        file_ext = os.path.splitext(os.path.splitext(file_name)[0])[1]

    if debug:
        print('file_type: ', file_type)
        print('file_name: ', file_name)
        print(' file_ext: ', file_ext)

    # find extension and switch in input or output files
    if file_type == 'input':
        if file_ext not in input_flags:
            raise KeyError(('input file extension ' +
                            file_ext +
                            ' not recognised'))
        else:
            return([('-' + input_flags[file_ext]), file_name])
    if file_type == 'output':
        if file_ext not in output_flags:
            raise KeyError(('output file extension ' +
                            file_ext +
                            ' not recognised'))
        else:
            return([('-' + output_flags[file_ext]), file_name])


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
    _allowed_job_types = ['transform', 'originate', 'download']
    if job_type not in _allowed_job_types:
        raise ValueError('{job_type} not an allowed job_type')

    # set up the args
    function_args = []
    # if we expect input_files, they go first
    if job_type == 'transform':
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

        # use io_file_to_bash_flag function to generate arguments to pass to
        # bash script.

        function_args_list = list(function_args)
        submit_args = []

        if verbose:
            print("\nfunction_args: ", function_args)
            print("\nfunction_args_list: ", function_args_list)

        # if we expect input_files, they go first
        if job_type == 'transform':
            input_files = [function_args_list.pop(0)]
            input_files_flat = list(flatten_list(input_files))
            input_args = [io_file_to_bash_flag(x, 'input')
                          for x in input_files_flat]
            input_args_flat = list(flatten_list(input_args))

            submit_args.append(input_args_flat)

            if verbose:
                print("\n input_files_flat:", input_files_flat)
                print("\n  input_args_flat: ", input_args_flat)
                print("\n      submit_args: ", submit_args)

        # output_files required for all job_types
        output_files = [function_args_list.pop(0)]
        output_files_flat = list(flatten_list(output_files))
        output_args = [io_file_to_bash_flag(x, 'output')
                       for x in output_files_flat]
        output_args_flat = list(flatten_list(output_args))

        submit_args.append(output_args_flat)

        if verbose:
            print("\noutput_files_flat:", output_files_flat)
            print("\n output_args_flat: ", output_args_flat)
            print("\n      submit_args: ", submit_args)

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
            print("remaining function_args_list: ", function_args_list)

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


def generate_queue_job_function(job_script, job_name, verbose=False):
    # type: (str, str) -> NoneType
    '''
    Run the HaplotypeCaller shell script, capture and mail output.
    '''

    def job_function(input_files, output_files):

        if verbose:
            print("\n\nGenerating Queue job function\n")

        submit_args = []

        # flatten file lists
        input_files_flat = list(flatten_list([input_files]))
        output_files_flat = list(flatten_list([output_files]))

        if verbose:
            print(" input_files_flat: ", input_files_flat)
            print("output_files_flat: ", output_files_flat)

        # construct argument list
        input_args = [io_file_to_bash_flag(x, 'input')
                      for x in input_files_flat]
        input_args_flat = list(flatten_list(input_args))

        submit_args.append(input_args_flat)

        output_args = [io_file_to_bash_flag(x, 'output')
                       for x in output_files_flat]
        output_args_flat = list(flatten_list(output_args))

        submit_args.append(output_args_flat)

        submit_args_flat = list(flatten_list([submit_args]))

        if verbose:
            print(" submit_args_flat: ", submit_args_flat)

        # run the job and grab output
        proc = subprocess.Popen(
            [job_script] + list(submit_args_flat),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = proc.communicate()

        # write output for emailing
        tmp_out = tempfile.mkstemp(
            prefix=(job_name + '.'), suffix=".out.txt", text=True)[1]
        tmp_err = tempfile.mkstemp(
            prefix=(job_name + '.'), suffix=".err.txt", text=True)[1]
        with open(tmp_out, 'wb') as f:
            f.write(out)
        with open(tmp_err, 'wb') as f:
            f.write(err)

        # mail output
        if proc.returncode != 0:
            subject = "[Tom@SLURM] Pipeline step " + job_name + " FAILED"
        else:
            subject = "[Tom@SLURM] Pipeline step " + job_name + " finished"
        mail = subprocess.Popen(['mail', '-s', subject, '-A', tmp_out, '-A',
                                 tmp_err, 'tom'], stdin=subprocess.PIPE)
        mail.communicate()
        os.remove(tmp_out)
        os.remove(tmp_err)

        # check subprocess exit code
        assert proc.returncode == 0, ("Job " + job_name +
                                      " failed with non-zero exit code")

    return job_function
