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

    # bamfiles
    raw_files = [x.path for x in os.scandir('data/bam') if
                 x.name.endswith('.bam') and x.is_file]
    mapped_raw = main_pipeline.originate(
        name='mapped_raw',
        task_func=os.path.isfile,
        output=raw_files)

    # mark duplicates with picard
    deduped = main_pipeline.transform(
        name='dedupe',
        task_func=functions.generate_job_function(
            job_script='src/sh/mark_duplicates_and_sort',
            job_name='dedupe',
            job_type='transform',
            cpus_per_task=2),
        input=mapped_raw,
        filter=ruffus.regex(r"data/bam/(.*).Aligned.out.bam"),
        output=[r"output/mark_duplicates_and_sort/\1.deduped.bam"])

    ###################
    # RUFFUS COMMANDS #
    ###################

    # print the flowchart
    ruffus.pipeline_printout_graph(
        "ruffus/flowchart.pdf", "pdf",
        pipeline_name="5 accessions variant calling pipeline")

    # run the pipeline
    ruffus.cmdline.run(options, multithread=8)

if __name__ == "__main__":
    main()
