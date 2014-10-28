#!/usr/bin/env python
# coding=utf-8
"""
Build and execute metaquast command for a run folder over contig files that
match a given pattern.
"""

from __future__ import print_function

import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from fnmatch import fnmatch
from os import path, walk
from subprocess import check_call


def find(pattern, folder):
    result = []
    for root, dirs, files in walk(folder):
        for name in files:
            if fnmatch(name, pattern):
                result.append(path.join(root, name))
    return result


def main(folder, ref, output_dir, contig_thresholds, min_contig, pattern, threads):
    # find the files
    contigs = find(pattern, folder)

    # in case this gets run in the wrong directory
    if len(contigs) < 1:
        sys.exit("No contig files were found using %s." % pattern)
    elif len(contigs) > 50:
        sys.exit("Found too many contig files to run quast.")

    # build the labels
    # assumes file names are separated by '_' and sample name comes first
    sep = '_'
    labels = [path.basename(c).split(sep, 1)[0] for c in contigs]

    # run metaquast
    cmd = ("quast.py -f -o {out} -R {ref} -M {min_contig} -T {cpus} "
        "-l {labels} -t {thresholds} {fastas}").format(out=output_dir,
        ref=ref, min_contig=min_contig, cpus=threads, labels=','.join(labels),
        thresholds=contig_thresholds, fastas=' '.join(contigs))
    check_call(cmd, shell=True)

    # combine reports into single tsv
    combined_report = path.join(output_dir, 'report_combined.tsv')
    report = path.join(output_dir, 'report.tsv')
    misassembly_report = path.join(output_dir, 'contigs_reports', 'misassemblies_report.tsv')
    unaligned_report = path.join(output_dir, 'contigs_reports', 'unaligned_report.tsv')

    with open(combined_report, 'w') as com_rep, open(report, 'r') as rep:

        [print(line.strip("\r\n"), file=com_rep) for line in rep]

        if path.exists(misassembly_report):
            with open(misassembly_report, 'r') as mis_rep:
                [print(line.strip("\r\n"), file=com_rep) for i, line in enumerate(mis_rep) if i > 0]

        if path.exists(unaligned_report):
            with open(unaligned_report, 'r') as una_rep:
                [print(line.strip("\r\n"), file=com_rep) for i, line in enumerate(una_rep) if i > 0]


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=ADHF)
    p.add_argument('folder',
        help='a parent directory containing contig fastas within its tree')
    p.add_argument('ref', help='reference fasta')
    p.add_argument('-o', '--output-dir', default='quast',
        help='quast output directory')
    p.add_argument('--contig-thresholds', default='2000,3000,5000',
        help='contig length thresholds')
    p.add_argument('--min-contig', default=2000,
        help='lower threshold for contig length')
    # might be more useful if this accepted multiple patterns
    p.add_argument('-p', '--pattern', default='*all_contigs.fasta',
        help='file search pattern to use from input-dir')
    p.add_argument('--threads', default=10, help='processing threads')
    args = vars(p.parse_args())
    main(**args)
