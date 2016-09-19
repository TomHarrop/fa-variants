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

    # test function for checking input/output passed to job_script and parsing
    # by io_parser
    test_job_function = functions.generate_job_function(
        job_script='src/sh/io_parser',
        job_name='test')

    # initialise pipeline
    main_pipeline = ruffus.Pipeline.pipelines["main"]

    # bamfiles
    raw_files = [x.path for x in os.scandir('data/bam') if
                 x.name.endswith('.bam') and x.is_file]

    # subset the files while the pipeline is in development. Make this equal
    # to the raw_files to run the whole pipline.
    # active_raw_files = [x for x in raw_files if
    #                     'G1' in x or 'G4' in x or 'J1' in x or 'J4' in x]
    active_raw_files = raw_files

    # species short names for vcf splitting
    species_short_names = list(set(
        [os.path.basename(x)[0] for x in active_raw_files]))

    # check that the files exist
    mapped_raw = main_pipeline.originate(
        name='mapped_raw',
        task_func=os.path.isfile,
        output=active_raw_files)

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

    # convert annotation to .bed
    annot_bed = main_pipeline.transform(
        name='annot_bed',
        task_func=functions.generate_job_function(
            job_script='src/sh/annot_bed',
            job_name='annot_bed',
            job_type='transform',
            cpus_per_task=7),
        input=annot,
        filter=ruffus.suffix('.gtf'),
        output='.bed')

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
        output=(r"output/mark_duplicates_and_sort/\1.deduped.bam"))

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

    # we're going to recycle call_variants, merge_variants, filter_variants
    # and analyze_covar so we'll get the functions in advance
    call_variants = functions.generate_queue_job_function(
        job_script='src/sh/call_variants',
        job_name='call_variants')
    merge_variants = functions.generate_job_function(
        job_script='src/sh/merge_variants',
        job_name='merge_variants',
        job_type='transform',
        cpus_per_task=8)
    filter_variants = functions.generate_job_function(
        job_script='src/sh/filter_variants',
        job_name='filter_variants',
        job_type='transform',
        cpus_per_task=1)
    analyze_covar = functions.generate_queue_job_function(
        job_script='src/sh/analyze_covar',
        job_name='analyze_covar')

    # call variants without recalibration tables
    uncalibrated_variants = main_pipeline.transform(
        name='uncalibrated_variants',
        task_func=call_variants,
        input=split_and_trimmed,
        add_inputs=ruffus.add_inputs([ref_fa, annot_bed]),
        filter=ruffus.formatter('output/split_trim/(?P<LIB>.+).split.bam'),
        output='{subdir[0][1]}/variants_uncalibrated/{LIB[0]}.g.vcf.gz')

    # merge gVCF variants
    uncalibrated_variants_merged = main_pipeline.merge(
        name='uncalibrated_variants_merged',
        task_func=merge_variants,
        input=[uncalibrated_variants, ref_fa],
        output='output/variants_uncalibrated/variants_uncalibrated.vcf.gz')

    # filter variants on un-corrected bamfiles
    uncalibrated_variants_filtered = main_pipeline.transform(
        name='uncalibrated_variants_filtered',
        task_func=filter_variants,
        input=uncalibrated_variants_merged,
        add_inputs=ruffus.add_inputs(ref_fa),
        filter=ruffus.suffix('_uncalibrated.vcf.gz'),
        output='_uncalibrated_filtered.vcf.gz')

    # select variant (only recalibrate using passed SNPs)
    uncalibrated_variants_selected = main_pipeline.transform(
        name='uncalibrated_variants_selected',
        task_func=functions.generate_job_function(
            job_script='src/sh/select_variants',
            job_name='select_variants',
            job_type='transform'),
        input=uncalibrated_variants_filtered,
        add_inputs=ruffus.add_inputs(ref_fa),
        filter=ruffus.suffix('_uncalibrated_filtered.vcf.gz'),
        output='_uncalibrated_selected.vcf.gz')

    # create recalibration report with filtered variants
    covar_report = main_pipeline.merge(
        name='covar_report',
        task_func=analyze_covar,
        input=[split_and_trimmed, ref_fa, annot_bed,
               uncalibrated_variants_selected],
        output="output/covar_analysis/recal_data.table")

    # second pass to analyze covariation remaining after recalibration
    second_pass_covar_report = main_pipeline.merge(
        name='second_pass_covar_report',
        task_func=analyze_covar,
        input=[split_and_trimmed, ref_fa, annot_bed,
               uncalibrated_variants_filtered, covar_report],
        output="output/covar_analysis/post_recal_data.table")

    # plot effect of base recalibration
    recal_plot = main_pipeline.transform(
        name='recal_plot',
        task_func=functions.generate_job_function(
            job_script='src/R/recal_plot.R',
            job_name='recal_plot',
            job_type='transform',
            cpus_per_task=1),
        input=second_pass_covar_report,
        filter=ruffus.suffix('post_recal_data.table'),
        add_inputs=ruffus.add_inputs(covar_report),
        output='recalibration_plots.pdf')

    # recalibrate bases using recalibration report
    recalibrated = main_pipeline.transform(
        name='recalibrate',
        task_func=functions.generate_job_function(
            job_script='src/sh/recalibrate',
            job_name='recalibrate',
            job_type='transform',
            cpus_per_task=2),
        input=split_and_trimmed,
        add_inputs=ruffus.add_inputs([ref_fa, covar_report]),
        filter=ruffus.formatter('output/split_trim/(?P<LIB>.+).split.bam'),
        output='{subdir[0][1]}/recal/{LIB[0]}.recal.bam')

    # final variant calling
    variants = main_pipeline.transform(
        name='variants',
        task_func=call_variants,
        input=recalibrated,
        add_inputs=ruffus.add_inputs(ref_fa, annot_bed),
        filter=ruffus.formatter('output/recal/(?P<LIB>.+).recal.bam'),
        output='{subdir[0][1]}/variants/{LIB[0]}.g.vcf.gz')

    # merge gVCF variants
    variants_merged = main_pipeline.merge(
        name='variants_merged',
        task_func=merge_variants,
        input=[variants, ref_fa],
        output='output/variants/variants.vcf.gz')

    # variant filtering
    variants_filtered = main_pipeline.transform(
        name='variants_filtered',
        task_func=filter_variants,
        input=variants_merged,
        add_inputs=ruffus.add_inputs(ref_fa),
        filter=ruffus.suffix('.vcf.gz'),
        output='_filtered.vcf.gz')

    # variants by species
    split_variants = main_pipeline.subdivide(
        name='split_variants',
        task_func=functions.generate_job_function(
            job_script='src/sh/split_variants',
            job_name='split_variants',
            job_type='transform',
            cpus_per_task=1,
            ntasks=len(species_short_names)),
        input=variants_filtered,
        filter=ruffus.formatter(),
        add_inputs=ruffus.add_inputs(ref_fa),
        output=[('output/split_variants/' + x + '.variants_filtered.vcf.gz')
                for x in species_short_names])

    # count variants per gene per species
    cds_variants = main_pipeline.transform(
        name='cds_variants',
        task_func=functions.generate_job_function(
            job_script='src/R/cds_variants.R',
            job_name='cds_variants',
            job_type='transform'),
        input=split_variants,
        add_inputs=ruffus.add_inputs([ref_fa, annot]),
        filter=ruffus.formatter(
            'output/split_variants/(?P<LIB>.+).variants_filtered.vcf.gz'),
        output='{subdir[0][1]}/cds_variants/{LIB[0]}.cds_variants.Rds')

    # merge counted variants
    variants_per_gene = main_pipeline.merge(
        name='cds_merge',
        task_func=functions.generate_job_function(
            job_script='src/R/cds_merge.R',
            job_name='cds_merge',
            job_type='transform'),
        input=cds_variants,
        output='output/cds_variants/cds_variants.Rds')

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
