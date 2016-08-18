#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

############################################
# five accessions variant calling pipeline #
############################################

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
        output=raw_files_subset)  # subset speficied here

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
            cpus_per_task=2),
        input=deduped,
        add_inputs=ruffus.add_inputs(ref_fa),
        filter=ruffus.formatter(
            "output/mark_duplicates_and_sort/(?P<LIB>.+).deduped.bam"),
        output=["{subdir[0][1]}/split_trim/{LIB[0]}.split.bam"])\
        .follows(fa_idx)


    # we're going to recycle call_variants and filter_variants, so let's get
    # the functions in advance
    call_variants = functions.generate_job_function(
        job_script='src/sh/call_variants',
        job_name='call_variants',
        job_type='merge',
        cpus_per_task=7)
    filter_variants = functions.generate_job_function(
        job_script='src/sh/filter_variants',
        job_name='filter_variants',
        job_type='transform',
        cpus_per_task=7)

    # call variants without recalibration tables
    uncalibrated_variants = main_pipeline.merge(
        name='uncalibrated_variants',
        task_func=call_variants,
        input=[split_and_trimmed, ref_fa],
        output='output/variants_uncalibrated/variants_uncalibrated.vcf')

    # filter variants on un-corrected bamfiles
    uncalibrated_variants_filtered = main_pipeline.transform(
        name='uncalibrated_variants_filtered',
        task_func=filter_variants,
        input=uncalibrated_variants,
        add_inputs=ruffus.add_inputs(ref_fa),
        filter=ruffus.suffix('_uncalibrated.vcf'),
        output='_uncalibrated_filtered.vcf')

    # recalibrate bases using filtered variants

    # call variants

    # filter variants

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
