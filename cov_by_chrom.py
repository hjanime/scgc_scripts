#!/usr/bin/env python
# coding=utf-8
"""
avg. per chrom coverage for a bam file. defaults to using:

bedtools genomecov -d -split -ibam <bam>
"""

import os
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from itertools import groupby
from toolshed import reader


def main(bam, split=True):
    assert os.path.exists(bam)
    hdr = ['chrom', 'pos', 'count']
    # ensure sorted by chrom, position
    cmd = '| bedtools genomecov -d -split -ibam {bam} | sort -k1,1 -k2,2n'.format(bam=bam)
    if not split:
        cmd = cmd.replace(' -split', '')

    for chrom, grouped in groupby(reader(cmd, header=hdr), lambda x: x['chrom']):
        total_count = 0
        chrom_length = 0

        for toks in grouped:
            total_count += int(toks['count'])
            # length is final position
            chrom_length = toks['pos']

        avgcov = total_count / float(chrom_length)
        print "\t".join([chrom, "%0.3f" % avgcov])


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
    p.add_argument('bam', help="coordinate sorted bam file")
    p.add_argument('--no-split', action='store_false', help="do not split reads based on CIGAR")

    args = p.parse_args()
    main(args.bam, args.no_split)
