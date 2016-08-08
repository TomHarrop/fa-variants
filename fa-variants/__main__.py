#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

############################
# Proforma ruffus pipeline #
############################

import functions
import ruffus
import os
import datetime
import re


def main():

    #########
    # SETUP #
    #########

    # catch jgi logon and password from cli
    parser = ruffus.cmdline.get_argparse(
        description='5 accessions variant calling pipeline.')
    parser.add_argument('--email', '-e',
                        help='Logon email address for JGI',
                        type=str,
                        dest='jgi_logon')
    parser.add_argument('--password', '-p',
                        help='JGI password',
                        type=str,
                        dest='jgi_password')
    options = parser.parse_args()
    jgi_logon = options.jgi_logon
    jgi_password = options.jgi_password

    ##################
    # PIPELINE STEPS #
    ##################

    # initialise pipeline
    main_pipeline = ruffus.Pipeline.pipelines["main"]

    # test originate job
    test_originate_files = ['ruffus/foo.txt', 'ruffus/bar.txt']
    test_originate = main_pipeline.originate(
        name='test_originate',
        task_func=functions.generate_job_function(
            job_script='src/test_originate',
            job_name='test_originate',
            job_type='originate'),
        output=test_originate_files)

    # test download job
    test_download = main_pipeline.originate(
        name='test_download',
        task_func=functions.generate_job_function(
            job_script='src/test_download',
            job_name='test_download',
            job_type='download'),
        output='ruffus/download.txt',
        extras=[jgi_logon, jgi_password])

    # test transform with multiple outputs (e.g. bamfile, FASTA etc)
    test_transform = main_pipeline.transform(
        name="test_transform",
        task_func=functions.generate_job_function(
            job_script='src/test_transform',
            job_name='test_transform',
            job_type='transform'),
        input=test_originate,
        filter=ruffus.suffix(".txt"),
        output=["_transformed.txt", "_transformed.bam"])

    # Transform ONLY the bam files produced by test_transform

    # The filtering here is a bit crazy. `input` has to be an object, not
    # ruffus.output_from(). `replace_inputs` should use `ruffus.inputs()` to
    # match the files, but `filter` has to match the first file produced by
    # the previous step, NOT necessarily the file that will be transformed!
    test_selective_transform = main_pipeline.transform(
        name="test_selective_transform",
        task_func=functions.generate_job_function(
            job_script='src/test_selective_transform',
            job_name='test_selective_transform',
            job_type='transform'),
        input=test_transform,
        replace_inputs=ruffus.inputs(r"\1.bam"),
        filter=ruffus.suffix(".txt"),
        output=".bof")

    test_merge = main_pipeline.merge(
        name='test_merge',
        task_func=functions.generate_job_function(
            job_script='src/test_merge',
            job_name='test_merge',
            job_type='merge'),
        input=test_transform,
        output='ruffus/foobar_merge.txt'
        )

    ###################
    # RUFFUS COMMANDS #
    ###################

    # print the flowchart
    ruffus.pipeline_printout_graph("ruffus/flowchart.pdf", "pdf",
                                   pipeline_name="5 accessions analysis "
                                                 "pipeline")

    # run the pipeline
    ruffus.cmdline.run(options, multithread=8)

if __name__ == "__main__":
    main()
