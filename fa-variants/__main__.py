#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

############################
# Proforma ruffus pipeline #
############################

import functions
import ruffus
import os


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

    # subset the files while the pipeline is in development
    raw_files_subset = [x for x in raw_files if
                        'G1' in x or 'G4' in x or 'J1' in x or 'J4' in x]

    # check that the files exist
    mapped_raw = main_pipeline.originate(
        name='mapped_raw',
        task_func=os.path.isfile,
        output=raw_files_subset)

    # genome fasta
    ref_fa = main_pipeline.originate(
        name='ref_fa',
        task_func=functions.generate_job_function(
            job_script='src/sh/download_genome',
            job_name='ref_fa',
            job_type='download'),
        output='data/genome/Osativa_323_v7.0.fa',
        extras=[jgi_logon, jgi_password])

    # indexes
    fa_idx = main_pipeline.transform(
        name='fa_idx',
        task_func=functions.generate_job_function(
            job_script='src/sh/fa_idx',
            job_name='fa_idx',
            job_type='transform',
            cpus_per_task=6),
        input=ref_fa,
        filter=ruffus.suffix(".fa"),
        output=['.dict', '.fa.fai'])

    # annotation
    annot = main_pipeline.originate(
        name='annot',
        task_func=functions.generate_job_function(
            job_script='src/sh/download_genome',
            job_name='annot',
            job_type='download'),
        output=('data/genome/'
                'Osativa_323_v7.0.gene_exons.gffread.rRNAremoved.gtf'),
        extras=[jgi_logon, jgi_password])

    # mark duplicates with picard
    deduped = main_pipeline.subdivide(
        name='dedupe',
        task_func=functions.generate_job_function(
            job_script='src/sh/mark_duplicates_and_sort',
            job_name='dedupe',
            job_type='transform',
            cpus_per_task=2),
        input=mapped_raw,
        filter=ruffus.regex(r"data/bam/(.*).Aligned.out.bam"),
        output=(r"output/mark_duplicates_and_sort/\1.deduped.bam",
                r"output/mark_duplicates_and_sort/\1.deduped.bai"))

    # Split'N'Trim and reassign mapping qualities
    split_and_trimmed = main_pipeline.transform(
        name='split_trim',
        task_func=functions.generate_job_function(
            job_script='src/sh/split_trim',
            job_name='split_trim',
            job_type='transform',
            cpus_per_task=1),
        input=deduped,
        add_inputs=ruffus.add_inputs(ref_fa),
        filter=ruffus.formatter(
            "output/mark_duplicates_and_sort/(?P<LIB>.+).deduped.bam"),
        output=["{subdir[0][1]}/split_trim/{LIB[0]}.split.bam"])\
        .follows(fa_idx)

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
