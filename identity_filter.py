#!/usr/bin/env python
# coding=utf-8

import click
import os.path as op
import sys
from collections import Counter, defaultdict
from numpy import mean
from pysam import Samfile
from subprocess import check_call


@click.command(help="""
Finding how many reads that mapped to a reference above an edit distance
threshold. Only the match section is analyzed. Hard clipped reads are ignored.
STDOUT contains some info on the process while <bam>_filtered.bam is written
to disk.
""")
@click.argument("bam", type=click.Path(exists=True))
@click.option("-i", "--identity-threshold", default=0.90, type=float,
    show_default=True, help="identity threshold between read and reference")
@click.option("-m", "--mapping-threshold", default=0.33, type=float,
    show_default=True, help="at least this much of the read must be matching (M)")
def main(bam, identity_threshold, mapping_threshold):
    try:
        check_call("samtools index %s" % bam, shell=True)
    except Exception:
        print >>sys.stderr, "Unable to create indexes for bam and fasta."
        sys.exit(1)

    remainder = 1 - identity_threshold
    filtered_bam = bam.rsplit(".bam", 1)[0] + "_filtered.bam"

    too_short = Counter()
    passing_aln = Counter()
    total_aln = Counter()

    fail_length = defaultdict(list)
    pass_length = defaultdict(list)

    with Samfile(bam, 'rb') as fh, Samfile(filtered_bam, 'wb', template=fh) as fo:
        total_mapped = fh.mapped
        processed = 1

        for chrom in fh.references:
            print >>sys.stderr, "processing reference", chrom

            for aln in fh.fetch(chrom):
                if processed % 100000 == 0:
                    sys.stderr.write("\r%d of %d" % (processed, total_mapped))
                    sys.stderr.flush()

                processed += 1
                total_aln[chrom] += 1

                # filter just in case
                if aln.is_unmapped: continue

                # 1:insertion; 2:deletion; 5:hardclip
                if not all(op != 5 for (op, l) in aln.cigar): continue

                # too much soft or hard clipping
                M = aln.opt('AS')
                if not M / float(aln.rlen) > mapping_threshold:
                    too_short[chrom] += 1
                    fail_length[chrom].append(M)
                    continue
                pass_length[chrom].append(M)

                # check edit distance
                if aln.opt('NM') / float(aln.alen) > remainder: continue

                passing_aln[chrom] += 1
                fo.write(aln)

        print "#Processed:", op.basename(bam)
        print "#Total mapped:", total_mapped
        print "#Identity_threshold:", identity_threshold
        print "#Mapping threshold:", mapping_threshold
        fields = ['Reference', 'Passing_reads', 'Percent_pass_by_reference',
            'Percent_pass_by_total', 'Short_mapped', 'Avg_passing_length',
            'Avg_failed_length']
        print "\t".join(fields)
        for chrom in fh.references:
            passing = passing_aln[chrom]
            perc_pass_ref = 0 if total_aln[chrom] == 0 else passing / float(total_aln[chrom]) * 100
            perc_pass_total = passing / float(total_mapped) * 100
            fields = [chrom, passing, perc_pass_ref, perc_pass_total,
                too_short[chrom], mean(pass_length[chrom]),
                mean(fail_length[chrom])]
            print "\t".join(map(str, fields))


if __name__ == '__main__':
    main()
