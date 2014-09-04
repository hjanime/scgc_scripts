#!/usr/bin/env python
# coding=utf-8
"""
Runs bcl2fastq creating fastqs and concatenates fastqs across lanes. Intended
to be used with NextSeq data and it does not do any cleanup! Original dumped
fastqs will remain along with all of the bcl files.
"""

import csv
import os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, REMAINDER
from sh import bcl2fastq, cat


def parse_samples(sample_sheet):
    """Parse SampleSheet.csv"""
    samples = []
    experiment = ""
    with open(sample_sheet, 'rU') as fh:
        reader = csv.reader(fh)
        start = False
        for toks in reader:
            try:
                if toks[0] == "Experiment Name":
                    experiment = toks[1]
                if toks[0] == "Sample_ID":
                    start = True
                    continue
            except IndexError:
                # blank lines are expected in this file
                pass
            if start:
                samples.append(toks[0])
    return samples, experiment


def main(runfolder, min_log_level, loading_threads, demultiplexing_threads,
            processing_threads, barcode_mismatches, args):
    samplesheet = "{runfolder}/SampleSheet.csv".format(runfolder=runfolder)
    assert os.path.exists(samplesheet), "%s was not found." % samplesheet
    samples, experiment = parse_samples(samplesheet)

    if not '-w' in args or '--writing-threads' in args:
        args.append('--writing-threads')
        args.append(len(samples))

    # call bcl2fastq
    bcl2fastq(runfolder_dir=runfolder, min_log_level=min_log_level,
        loading_threads=loading_threads, demultiplexing_threads=demultiplexing_threads,
        processing_threads=processing_threads, barcode_mismatches=barcode_mismatches,
        *args, _no_err=True)

    # location of fastqs from bcl2fastq
    if '-o' in args:
        output_dir = args[args.index('-o') + 1]
    elif '--output-dir' in args:
        output_dir = args[args.index('--output-dir') + 1]
    else:
        output_dir = "{runfolder}/Data/Intensities/BaseCalls".format(runfolder=runfolder)

    # cat files across lanes into new file per sample
    processes = []
    for i, sample in enumerate(samples, start=1):
        # bcl2fastq converts underscores (from samplesheet) to dashes
        sample = sample.replace("_", "-")
        for read in ['R1', 'R2']:
            paths = []
            result_file=("{output_dir}/{experiment}/{sample}_{read}.fastq.gz").format(
                output_dir=output_dir, experiment=experiment, sample=sample, read=read)
            for lane in xrange(1, 5):
                # build the files paths
                path = ("{output_dir}/{experiment}/{sample}_S{sample_num}"
                        "_L00{lane}_{read}_001.fastq.gz").format(output_dir=output_dir,
                            experiment=experiment, sample=sample, sample_num=i, lane=lane, read=read)
                paths.append(path)
            # execute file concatenation
            processes.append(cat(paths, _out=result_file, _no_err=True, _bg=True))
    for p in processes:
        p.wait()


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('runfolder', help="path to run folder")
    p.add_argument('-l', '--min-log-level', choices=['NONE', 'FATAL', 'ERROR',
                                                     'WARNING', 'INFO', 'DEBUG', 'TRACE'],
                    default='INFO', help="minimum log level")
    p.add_argument('-r', '--loading-threads', default=12, type=int,
                    help="threads used for loading BCL data")
    p.add_argument('-d', '--demultiplexing-threads', default=12, type=int,
                    help="threads used for demultiplexing")
    p.add_argument('-p', '--processing-threads', default=12, type=int,
                    help="threads used for processing demultiplexed data")
    p.add_argument('--barcode-mismatches', default=0, type=int,
                    help="number of allowed mismatches per index")
    p.add_argument('args', nargs=REMAINDER,
                    help="any additional bcl2fastq args and their values")
    args = vars(p.parse_args())
    main(**args)
