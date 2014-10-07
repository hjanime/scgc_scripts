#!/usr/bin/env python
# coding=utf-8
"""
Runs bcl2fastq creating fastqs and concatenates fastqs across lanes. Intended
to be used with NextSeq data and it does not do any cleanup! Original dumped
fastqs will remain along with all of the bcl files.
"""

from __future__ import print_function

import fileinput
import os
import subprocess as sp
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, REMAINDER
from itertools import izip_longest


def process_samplesheet(samplesheet):
    samples = []
    experiment = ""

    try:
        start = False
        # strip whitespace and rewrite file in place
        for toks in fileinput.input(samplesheet, mode='rU', backup='.bak', inplace=True):
            toks = toks.rstrip("\r\n").split(',')
            print(",".join([t.strip() for t in toks]))
            if toks[0] == "Sample_ID":
                start = True
                continue
            if start:
                # bcl2fastq converts underscores to dashes
                samples.append(toks[0].replace("_", "-").replace(".", "-"))
                # location of fastq output
                experiment = toks[8]
    finally:
        fileinput.close()

    return samples, experiment


def main(runfolder_dir, loading_threads, demultiplexing_threads,
            processing_threads, barcode_mismatches, args):

    # verify working directory
    samplesheet = os.path.join(runfolder_dir, "SampleSheet.csv")
    if not os.path.exists(samplesheet):
        sys.exit("%s was not found. Exiting." % samplesheet)

    # get sample names and experiment name from SampleSheet
    samples, experiment = process_samplesheet(samplesheet)
    if len(samples) == 0:
        sys.exit("No samples were found in the SampleSheet. Check formatting.")

    # set args for call to bcl2fastq
    if not '-w' in args or '--writing-threads' in args:
        args.append('--writing-threads')
        args.append(len(samples))
    args.extend(['-R', runfolder_dir, '-r', loading_threads, '-d',
        demultiplexing_threads, '-p', processing_threads,
        '--barcode-mismatches', barcode_mismatches])
    args.insert(0, 'bcl2fastq')

    # run bcl2fastq
    cmd = " ".join(map(str, args))
    print("Converting .bcl to .fastq using:")
    print("$> ", cmd)
    sp.check_call(cmd, shell=True)
    print(".bcl conversion was successful")

    # location of fastqs from bcl2fastq
    if '-o' in args:
        output_dir = args[args.index('-o') + 1]
    elif '--output-dir' in args:
        output_dir = args[args.index('--output-dir') + 1]
    else:
        output_dir = os.path.join(runfolder_dir, "Data", "Intensities", "BaseCalls")
    output_dir = os.path.join(output_dir, experiment)

    # build `cat` commands to join files across lanes
    commands = []
    for i, sample in enumerate(samples, start=1):
        for read in ['R1', 'R2']:
            cmd = ['cat']
            result_file = "%s/%s_%s.fastq.gz" % (output_dir, sample, read)
            for lane in [1, 2, 3, 4]:
                # build the files paths
                path = "%s/%s_S%d_L00%d_%s_001.fastq.gz" % \
                    (output_dir, sample, i, lane, read)
                if not os.path.exists(path):
                    sys.exit("Can't find %s. Concatenation failed." % path)
                cmd.append(path)

            # using shell to ease redirection
            commands.append(" ".join(cmd) + " > " + result_file)

    # execute concatenation 4 samples at a time
    print("Joining reads across lanes")
    groups = [(sp.Popen(cmd, shell=True) for cmd in commands)] * 4
    for processes in izip_longest(*groups):
        for p in filter(None, processes):
            p.wait()


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('-R', '--runfolder-dir', default=".", help="path to run folder")
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
