#!/usr/bin/env python
# coding=utf-8

from __future__ import print_function

import click
import os
import parmap
import sys
from itertools import count, groupby, izip
from multiprocessing import Pool
from numpy import mean
from signal import signal, SIGPIPE, SIG_DFL
from toolshed import nopen

# no IOError on pipe to head
signal(SIGPIPE,SIG_DFL)


@click.group()
@click.version_option('0.4.0')
@click.pass_context
def cli(obj):
    """Fasta and fastq tools."""


def readfx(fastx):
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


def stats_from_fasta(fasta):
    """returns number of contigs, their sizes as a list, and GC total"""
    total_contigs = 0
    contig_sizes = []
    gc_total = 0
    for i, (name, seq, qual) in enumerate(readfx(fasta), start=1):
        contig_sizes.append(len(seq))
        gc_total += seq.count('G') + seq.count('C')
        total_contigs = i
    return total_contigs, contig_sizes, gc_total


def print_fasta_record(name, seq, fh, wrap=60):
    print(">%s" % name, file=fh)
    if wrap:
        for i in xrange(0, len(seq), wrap):
            print(seq[i:i + wrap], file=fh)
    else:
        print(seq, file=fh)


@cli.command('assembly-stats', short_help='basic fasta stats')
@click.argument('fasta', type=click.Path(exists=True))
@click.option('-o', '--out', help='when None write to STDOUT, else append to this file')
def assembly_stats(fasta, out):
    """
    Report total bases, max contig size, number of contigs,
    GC percent, N50, and mean contig length to <out>.
    """
    total_contigs, contig_sizes, gc_total = stats_from_fasta(fasta)

    if total_contigs > 0:
        contig_sizes.sort(reverse=True)
        maxcontigsize = contig_sizes[0]
        mincontigsize = contig_sizes[-1]
        total_bases = sum(contig_sizes)
        gc_percent = float(gc_total) / total_bases * 100

        if total_contigs == 1:
            nfifty = contig_sizes[0]
            avg = contig_sizes[0]
        else:
            avg = mean(contig_sizes)
            testsum = 0
            half = total_bases / 2.
            for size in contig_sizes:
                testsum += size
                if testsum >= half:
                    nfifty = size
                    break

    else:
        total_bases = 0
        mincontigsize = 0
        maxcontigsize = 0
        gc_percent = 0
        nfifty = 0
        avg = 0

    ofh = open(out, 'a') if out else sys.stdout
    print(total_bases, maxcontigsize, total_contigs, gc_percent, nfifty, avg, sep=",", file=ofh)


@cli.command('length-filter', short_help="filter fasta by seq length")
@click.argument('fasta')
@click.option('--length', default=2000, type=int, show_default=True,
    help='passing sequence length')
@click.option('--wrap', default=60, type=int, show_default=True,
    help='length of fasta sequence per line in output')
@click.option('--operation', default='min', type=click.Choice(['min', 'max']),
    show_default=True, help='operation to perform using threshold <length>')
def length_filter(fasta, length, wrap, operation):
    """Filter <fasta> sequences by <length>."""
    import operator
    op_func = operator.gt if operation == 'max' else operator.lt

    path = 'length_filter' if fasta == "-" or fasta == "stdin" else fasta.rsplit(".", 1)[0]
    failed = "%s.failed.%dbp_%s.fasta" % (path, length, operation)

    with open(failed, 'w') as failfh:
        for name, seq, qual in readfx(fasta):
            if op_func(len(seq), length):
                print_fasta_record(name, seq, failfh, wrap)
            else:
                print_fasta_record(name, seq, sys.stdout, wrap)


def complexity_filter_sequence(seq, threshold=0.05, alphabet='ACGT'):
    """
    >>> complexity_filter_sequence("AAACCCGGGT", threshold=0.10)
    False
    >>> complexity_filter_sequence("AAACCCGGGT", threshold=0.09)
    True
    """
    seqlen = float(len(seq))
    return False if any([seq.count(x) / seqlen <= threshold for x in
        alphabet]) else True


def complexity_filter_se_record((name, seq, qual), threshold, alphabet):
    if complexity_filter_sequence(seq, threshold, alphabet):
        if qual:
            records = "@%s\n%s\n+\n%s" % (name, seq, qual)
        else:
            records = ">%s\n%s" % (name, "\n".join([seq[i:i+60] for i in
                xrange(0, len(seq), 60)]))
        return records
    return ""


def complexity_filter_pe_record(((n1, s1, q1), (n2, s2, q2)), threshold, alphabet):
    assert n1 == n2, "Pairing failed due to sync: %s %s" % (n1, n2)
    pass1 = complexity_filter_sequence(s1, threshold, alphabet)
    pass2 = complexity_filter_sequence(s2, threshold, alphabet)

    if pass1 and pass2:
        if not n1.endswith('/1'): n1 += '/1'
        if not n2.endswith('/2'): n2 += '/2'
        records = "@%s\n%s\n+\n%s\n@%s\n%s\n+\n%s" % (n1, s1, q1, n2, s2, q2)
        return records
    return ""


def multiprocess(f, iterable, *args, **kwargs):
    """
    Map an iterable to a function. Default key function chunks iterable by
    1000s.

    :param f: function
    :param iterable: any iterable where each item is sent to f
    :param *args: arguments passed to mapped function
    :param **kwargs: additional arguments for parmap.map
    """
    chunksize = kwargs.pop('chunksize', 1000)
    key = kwargs.pop('key', lambda k, l=count(): next(l)//chunksize)
    for k, g in groupby(iterable, key=key):
        yield parmap.map(f, g, *args, **kwargs)


@cli.command('complexity-filter', short_help='fastx complexity filter')
@click.argument('fastx', nargs=-1)
@click.option('-t', '--threshold', default=0.05, type=float, show_default=True,
    help='expected fraction threshold')
@click.option('-a', '--alphabet', default='ACTG', show_default=True,
    help='nucleotides or other acceptable letters')
@click.option('-p', '--pool', default=1, type=int, show_default=True,
    help='cpu pool size to use')
@click.option('--sample', default=None,
    help='optional name used rather than <fastx>; useful when using stdin')
def complexity_filter(fastx, threshold, alphabet, pool, sample):
    """
    Remove low complexity reads (each letter in <alphabet> above <threshold>).
    When running paired-ends read, both can be specified as fastx and are
    interwoven to stdout.
    """

    if len(fastx) > 2:
        print("Unexpected number of files to process.", file=sys.stderr)
        sys.exit(1)

    fname = sample if sample else (",".join(os.path.basename(x) for x in fastx))
    total_records = 0
    filtered_records = 0
    p = Pool(pool)

    if len(fastx) == 1:
        for results in multiprocess(complexity_filter_se_record, readfx(fastx[0]),
                threshold, alphabet, pool=p):
            for record in results:
                total_records += 1
                if not record:
                    filtered_records += 1
                    continue
                print(record)

    else:
        for results in multiprocess(complexity_filter_pe_record,
                izip(readfx(fastx[0]), readfx(fastx[1])), threshold, alphabet, pool=p):
            for record in results:
                total_records += 2
                if not record:
                    filtered_records += 2
                    continue
                print(record)

    print(("#Complexity filtered input (threshold=%0.2f,alphabet=%s)"
        "\ttotal\tpassing\n%s\t%d\t%d") % (threshold, alphabet, fname,
        total_records, total_records - filtered_records), file=sys.stderr)


@cli.command('header', short_help='prepend or append to header')
@click.argument('fastx')
@click.argument('text', type=str)
@click.option('--intent', default='prepend',
    type=click.Choice(['prepend', 'append']), show_default=True,
    help='action to perform')
@click.option('--sep', default='_', show_default=True, help='field separator')
def munge_header(fastx, text, intent, sep):
    """
    Prepend or append <text> to name field using <sep>.
    """
    for name, seq, qual in readfx(fastx):

        if intent == 'prepend':
            name = "%s%s%s" % (text, sep, name)
        else:
            name = "%s%s%s" % (name, sep, text)

        # fastq
        if qual:
            print("@" + name, seq, "+", qual, sep="\n")
        # fasta
        else:
            print_fasta_record(name, seq, sys.stdout)


def seq_sliding_gc(seq, window_size):
    """sliding window implementation of gc_skew and gc_content.

    >>> seq = 'TTATCCTATTCAGCCACCCGATTTGGAAACCCGGATTGCAATTCTCCAAA'
    >>> window_size = 40
    >>> vals = [(mid, skew, content) for mid, skew, content in gc_skew_and_content(seq, window_size)]
    >>> vals[0][0]
    20
    >>> vals[0][1]
    -0.263...
    >>> vals[0][2]
    0.47...
    """
    seq_len = len(seq)
    if seq_len < window_size:
        yield None, None, None

    half_window = window_size / 2
    seq = seq.upper()

    for start in xrange(0, seq_len - window_size + 1):
        stop = start + window_size
        mid = stop - half_window
        s = seq[start:stop]

        g = s.count('G')
        c = s.count('C')
        gc = g + c

        content = float(gc) / window_size
        skew = 0 if gc == 0 else (g - c) / float(g + c)
        yield mid, skew, content


@cli.command('sliding-gc', short_help='sliding GC and skew calculations')
@click.argument('fasta')
@click.option('--window_size', default=500, type=int, show_default=True,
    help='size of sliding window')
def sliding_gc(fasta, window_size):
    """
    Calculate sliding GC content and skew for sequences over <window_size>.
    """
    header = ["name", "position", "skew", "content"]
    print(*header, sep="\t")
    for name, seq, qual in readfx(fasta):
        for point, skew, content in seq_sliding_gc(seq, 500):
            if point:
                print("%s\t%i\t%0.3f\t%0.3f" % (name, point, skew, content))
            else:
                print("skipped (seq length < window_size): %s" % name, file=sys.stderr)


@cli.command('count', short_help='count the reads')
@click.argument('fastx', type=click.Path(exists=True))
def read_count(fastx):
    """
    count the number of reads present in fastq or fasta.
    """

    total = 0
    fq = True
    for name, seq, qual in readfx(fastx):
        if not qual:
            fq = False
        break

    if fastx.endswith("gz"):
        count_file = fastx.rsplit(".gz", 1)[0] + ".count"
        cat = "gunzip -c"
    else:
        count_file = fastx + ".count"
        cat = "cat"

    if os.path.exists(count_file):
        cmd = "cat %s" % count_file
    elif not fq:
        cmd = '%s %s | grep -c "^>"' % (cat, fastx)
    else:
        cmd = "%s %s | wc -l" % (cat, fastx)

    for count in nopen("|" + cmd):
        total = int(count)
        if fq and not cmd.endswith('.count'):
            assert total % 4 == 0
            total = total / 4

    if not os.path.exists(count_file):
        with open(count_file, 'w') as fh:
            print (total, file=fh)

    print(total)


@cli.command('split-merged', short_help='unmerge interweaved fastq file')
@click.argument('fastq')
@click.argument('r1', type=str)
@click.argument('r2', type=str)
@click.option('threads', type=int, default=8)
def split_merged(fastq, r1, r2, threads):
    """
    Unmerges <fastq> into gzipped (uses pigz!) <r1> and <r2>.
    """
    import subprocess

    assert r1 != r2
    r1 = r1[:-3] if r1.endswith('.gz') else r1
    r2 = r2[:-3] if r2.endswith('.gz') else r2

    print("Splitting reads (1/2).")
    cat = "gunzip -c %s" % fastq if fastq.endswith('gz') else "cat %s" % fastq
    cmd = """%s | awk 'BEGIN{o=1}{r++;if(o%%2==0){r2++; print $0 > "%s"}\
        else{r1++;print $0 > "%s"}if(r%%4==0){o++}}END{print "R1 records: "r1/4\
        "\\nR2 records: "r2/4}'""" % (cat, r2, r1)
    subprocess.check_call(cmd, shell=True)
    print("Gzipping reads (2/2).")
    cmd = "pigz -fp%d %s %s" % (threads, r1, r2)
    subprocess.check_call(cmd, shell=True)


@cli.command('merge-pe', short_help='merge R1 and R2 into interleaved fastq')
@click.argument('r1', type=click.Path(exists=True))
@click.argument('r2', type=click.Path(exists=True))
@click.option('-n','--name-trim', default=2, type=int, show_default=True,
    help='number of characters to remove from name before comparing R1 and R2, eg. 2 for /1 /2 convention')
def merge_pe(r1, r2, name_trim):
    """
    Merge <r1> and <r2> into interleaved fastq. This function (1) assumes you're
    dealing with small fastqs as they are read into memory (2) you don't care
    about the singletons.
    """

    def fq2dict(fq, n):
        d = {}
        s = set()
        for name, seq, qual in readfx(fq):
            d[name[:-n]] = (name, seq, qual)
            s.add(name[:-n])
        return d, s

    r1d, r1names = fq2dict(r1, name_trim)
    r2d, r2names = fq2dict(r2, name_trim)
    i = r1names & r2names
    remaining = 0

    for name in i:
        remaining += 2
        n1, s1, q1 = r1d[name]
        n2, s2, q2 = r2d[name]
        if not n1.endswith('/1'): n1 += '/1'
        if not n2.endswith('/2'): n2 += '/2'
        print('@' + n1, s1, '+', q1, '@' + n2, s2, '+', q2, sep="\n")

    print("Total: %d" % (len(r1names) + len(r2names)),
        "Remaining: %d" % remaining, sep=",", file=sys.stderr)


if __name__ == '__main__':
    cli()
