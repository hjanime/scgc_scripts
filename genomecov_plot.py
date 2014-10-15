#!/usr/bin/env python
# coding=utf-8
"""
Coverage plot over entire genome. Intended for single chromosomal organisms.
"""

from __future__ import print_function

import os
import subprocess
import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from rpy2 import robjects
from rpy2.robjects.vectors import FloatVector
from rpy2.robjects.packages import importr
from toolshed import reader


grdevices = importr('grDevices')
smooth = robjects.r['loess.smooth']
plot = robjects.r['plot']


def check_input(input_file):
    bedgraph = ""
    filename, ext = os.path.splitext(input_file)
    if filename.endswith(".bedgraph") or filename.endswith(".bg") \
            or ext == ".bedgraph" or ext == ".bg":
        bedgraph = input_file
    elif ext == ".bam":
        bedgraph = filename + ".bedgraph.gz"
        if not os.path.exists(bedgraph):
            print("Finding coverages")
            cmd = ("bedtools genomecov -bga -ibam {bam} "
                "| bedtools sort -i - "
                "| gzip > {bedgraph}").format(bam=input_file, bedgraph=bedgraph)
            subprocess.check_call(cmd, shell=True)
    else:
        sys.exit("Unable to determine input file type.")
    return bedgraph


def main(input_file, output_dir, sample, reference):
    bedgraph = check_input(input_file)
    plot_file = os.path.join(output_dir, "%s_%s.png" % (sample, reference))

    print("Reading count data", file=sys.stderr)
    pos = []
    cnt = []
    for toks in reader(bedgraph, header=['chrom', 'start', 'stop', 'count']):
        pos.append(toks['start'])
        cnt.append(toks['count'])

    print("Generating plot", file=sys.stderr)
    prediction = smooth(FloatVector(pos), FloatVector(cnt), span=0.01,
        degree=0, family="gaussian", evaluation=5000)
    grdevices.png(plot_file, width=1200, height=800)
    plot(prediction, type="l", xlab="Position", ylab="Count", bty="n",
        main="%s\n Genomic Coverage" % reference)
    grdevices.dev_off()


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=ADHF)
    p.add_argument('input', help='bam or bedgraph file')
    p.add_argument('output', help='location to write the plot file')
    p.add_argument('sample', help='name of the input sample')
    p.add_argument('reference', help='reference name to which sample is aligned')

    args = p.parse_args()

    if not os.path.exists(args.input):
        sys.exit("The input file does not exist. Exiting.")
    if not os.path.exists(args.output):
        os.makedirs(args.output)

    main(args.input, args.output, args.sample, args.reference)
