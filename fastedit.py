#!/usr/bin/env python
# coding=utf-8

from __future__ import print_function

import click
import sys
from csv import reader
from toolshed import nopen


@click.group()
@click.version_option('0.1.0')
@click.pass_context
def cli(obj):
    """2-step fasta header editing"""


def readfx(fastx):
    """Fasta and fastq parser."""
    with nopen(fastx) as fp:
        last = None
        while True:
            if not last:
                for l in fp:
                    if l[0] in '>@':
                        last = l[:-1]
                        break
            if not last: break
            name, seqs, last = last[1:].partition(" ")[0], [], None
            for l in fp:
                if l[0] in '@+>':
                    last = l[:-1]
                    break
                seqs.append(l[:-1])
            if not last or last[0] != '+':
                yield name, ''.join(seqs), None
                if not last: break
            else:
                seq, leng, seqs = ''.join(seqs), 0, []
                for l in fp:
                    seqs.append(l[:-1])
                    leng += len(l) - 1
                    if leng >= len(seq):
                        last = None
                        yield name, seq, ''.join(seqs);
                        break
                if last:
                    yield name, seq, None
                    break


@cli.command('get', short_help='retrieve fasta headers')
@click.argument('fasta', type=click.Path(exists=True))
@click.option('-o', '--out', help='output file: default writes to STDOUT')
def get_headers(fasta, out):
    """Pull out the headers from <fasta> and write them to <out>."""
    total = 0
    with open(out, 'w') if out else sys.stdout as fh:
        for name, seq, qual in readfx(fasta):
            total += 1
            print (name, file=fh)
    print (total, "headers exported.", file=sys.stderr)


@cli.command('put', short_help='replace fasta headers')
@click.argument('fasta', type=click.Path(exists=True))
@click.argument('csv', type=click.Path(exists=True))
@click.option('-o', '--out', help='default writes to STDOUT')
def put_headers(fasta, csv, out):
    """
    Replace <fasta> headers with column 2 of <csv> where column 1
    matches the original header in <fasta>. Headers will only be altered for
    <csv> entries with 2 columns. Rows removed from <csv> will be removed
    from the output.
    """
    headers = {}
    with open(csv, 'rU') as fh:
        it = reader(fh)
        for r in it:
            try:
                headers[r[0]] = r[1][1:] if r[1].startswith(">") else r[1]
            except IndexError:
                headers[r[0]] = ''

    total = 0
    renamed = 0
    ignored = 0
    with open(out, 'w') if out else sys.stdout as fh:
        for name, seq, qual in readfx(fasta):
            total += 1
            try:
                newname = headers[name]
                if newname:
                    renamed += 1
                else:
                    # in the csv, but no new name
                    newname = name

                # print the record
                print (">%s" % newname, file=fh)
                for i in xrange(0, len(seq), 60):
                    print (seq[i:i + 60], file=fh)

            except KeyError:
                # removed from the csv
                ignored += 1

    print ("Input fasta contained", total, "records", file=sys.stderr)
    print ("You renamed", renamed, "and kept a total of", total - ignored,
        file=sys.stderr)


if __name__ == '__main__':
    cli()
