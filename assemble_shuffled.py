#!/usr/bin/env python
# coding=utf-8
"""
assemble interweaved, paired-end reads; the result of shuffleSequences_fastq.pl.
"""

import logging
import os
import os.path as op
import shutil
import sys
import tempfile as tf
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from itertools import groupby, islice
from subprocess import call, CalledProcessError, Popen
from toolshed import nopen


# R function wanted a log file name
LOG = None


def verbose_logging(**kwargs):
    logging.info("Record of arguments for %s" % os.path.basename(__file__))
    for k, v in kwargs.items():
        logging.info(">>> %s = %s", k, v)


def get_sample(filename):
    name, ext = op.splitext(filename)
    if not ext == ".fastq" and not ext == ".fq":
        # gzipped input
        name, ext = op.splitext(name)
        if not ext == ".fastq" and not ext == ".fq":
            logging.error("Looking for .fastq file; found %s.%s", name, ext)
            raise IOError
    return name


def readfa(fa):
    with nopen(fa) as fh:
        for header, group in groupby(fh, lambda line: line[0] == '>'):
            if header:
                line = group.next()
                name = line[1:].strip()
            else:
                seq = ''.join(line.strip() for line in group)
                yield name, seq


def readfq(fq):
    with nopen(fq) as fh:
        fqclean = (x.strip("\r\n") for x in fh if x.strip())
        while True:
            r = [x for x in islice(fqclean, 4)]
            if not r: raise StopIteration
            assert all(r) and len(r) == 4
            yield r[0][1:], r[1], r[3]


def fqreads(fq):
    # `wc -l / 4` would be faster
    total = 0
    for name, seq, qual in readfq(fq):
        total += 1
    return total


def copyfiles(src, dst, **kwargs):
    cmd = "rsync" + kwargs_to_flag_string(kwargs) + " %s/* %s" % (src, dst)
    return runcmd(cmd)


def copytree(src, dst, symlinks=False, ignore=None):
    """copying to a cifs mount: copy2 will not work."""
    if not op.exists(dst):
        os.makedirs(dst)
        shutil.copystat(src, dst)
    lst = os.listdir(src)
    if ignore:
        excl = ignore(src, lst)
        lst = [x for x in lst if x not in excl]
    for item in lst:
        s = op.join(src, item)
        d = op.join(dst, item)
        if symlinks and op.islink(s):
            if op.lexists(d):
                os.remove(d)
            os.symlink(os.readlink(s), d)
            try:
                st = os.lstat(s)
                mode = stat.S_IMODE(st.st_mode)
                os.lchmod(d, mode)
            except:
                pass # lchmod not available
        elif op.isdir(s):
            copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)


def runcmd(cmd, log=True):
    """
    >>> import sys
    >>> import logging
    >>> from subprocess import call, CalledProcessError, Popen
    >>> logging.basicConfig(stream=sys.stdout)
    >>> runcmd("date > /dev/null", False)
    0
    >>> runcmd(["date > /dev/null", "date > /dev/null"], False)
    0
    """
    if isinstance(cmd, basestring):
        if log: logging.debug("Running command: %s", cmd)
        returncode = call(cmd, shell=True)
        if not returncode == 0:
            logging.error("Execution failed.")
            raise CalledProcessError

    else:
        assert isinstance(cmd, list)
        if log: logging.debug("Running commands: %s", cmd)

        # this will be bad if the list of commands is really long
        processes = [Popen(c, shell=True) for c in cmd]
        for p in processes:
            p.wait()
            returncode = p.returncode
            if returncode != 0:
                logging.error("Execution failed.")
                raise CalledProcessError

    return returncode


def kwargs_to_flag_string(kwargs, ignore=[]):
    """
    >>> kwargs = {"d":"d_test", "p":"p_test", "sep":","}
    >>> kwargs_to_flag_string(kwargs, ignore=['sep'])
    ' -p p_test -d d_test'
    """
    s = ""
    for k, v in kwargs.items():
        if k in ignore: continue
        s += (" -" if len(k) == 1 else " --") + k + ("" if v is None else (" " + str(v)))
    return s


def run_kmernorm(fastq, **kwargs):
    outfile = fastq.rsplit(".", 1)[0] + ".kmernorm.fastq"
    kwargs['k'] = 19 if not 'k' in kwargs else kwargs['k']
    kwargs['t'] = 80 if not 't' in kwargs else kwargs['t']
    kwargs['c'] = 2 if not 'c' in kwargs else kwargs['c']
    cmd = "kmernorm" + kwargs_to_flag_string(kwargs) + " " + fastq + " > " + outfile
    runcmd(cmd)
    return outfile


def lowcomp_filter(**kwargs):
    out = kwargs['input'].rsplit(".", 1)[0] + ".compfilter.fastq"
    kwargs['output'] = out
    # need to check the dependencies of this script
    cmd = "Rscript --vanilla /opt/scgc/rscript/fastq-low-complexity-filter.Rscript" \
                + kwargs_to_flag_string(kwargs)
    runcmd(cmd)
    return out


def spades(**kwargs):
    out = kwargs['o']
    # written for spades 3.0.0; installed in /usr/local/bin
    cmd = "spades.py" + kwargs_to_flag_string(kwargs)
    contigs = op.join(out, "contigs.fasta")
    runcmd(cmd)
    return contigs


def postprocess_spades(fasta, sample, outdir):
    out = op.join(outdir, sample + ".fasta")
    with open(out, 'wb') as fo, nopen(fasta) as fh:
        for line in fh:
            line = line.strip("\r\n")
            if line.startswith(">"):
                line = ">" + sample + "_" + line[1:]
            print >>fo, line
    return out


def assembly_stats(fasta):
    """calculate min_contig_size, max_contig_size, N50, total_bases,
    num_contigs, and gc_percent. writes to log.
    """
    contig_sizes = []
    gc_total = 0

    for name, seq in readfa(fasta):
        contig_sizes.append(len(seq))
        gc_total += seq.count('G') + seq.count('C')

    contigs = len(contig_sizes)

    if contigs > 0:
        contig_sizes.sort(reverse=True)
        maxcontigsize = contig_sizes[0]
        mincontigsize = contig_sizes[-1]
        total_bases = sum(contig_sizes)
        gc_percent = float(gc_total) / total_bases * 100

        if contigs == 1:
            nfifty = contig_sizes[0]
        else:
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

    logging.info("Stats for FASTA: %s", op.basename(fasta))
    logging.info("Total Bases: %s", total_bases)
    logging.info("Number of Contigs: %s", contigs)
    logging.info("Min Contig Size: %s", mincontigsize)
    logging.info("Max Contig Size: %s", maxcontigsize)
    logging.info("GC Percent: %s", gc_percent)
    logging.info("N50: %s", nfifty)


def read_counts(inputfq, kmer, compfiltered):
    count = fqreads(inputfq)
    logging.info("Input reads: %s", count)
    count = fqreads(kmer)
    logging.info("Normalized read count: %s", count)
    count = fqreads(compfiltered)
    logging.info("Low Complexity Filtered reads: %s", count)


def filter_fasta_by_size(fasta, size=2000):
    path = fasta.rsplit(".", 1)[0]
    passed = "%s.passed.%dbp_min.fasta" % (path, size)
    failed = "%s.failed.%dbp_min.fasta" % (path, size)
    with open(passed, 'w') as passfh, open(failed, 'w') as failfh:
        for name, seq in readfa(fasta):
            if len(seq) > size:
                # may need to word wrap for some applications
                print >>passfh, ">%s\n%s" % (name, seq)
            else:
                print >>failfh, ">%s\n%s" % (name, seq)
    return passed


def send_email(to="scgc@bigelow.org", subject="", message="", attachment=None):
    """ Send a small message via email
        @param to the email address to send the message to
        @param subject the subject line
        @param message a single line message for the email body
        @param attachment a fully qualifed filename of a file to attach (None by default)
        @return the return code for mutt application
        """
    if isinstance(message, list):
        message = " ".join(message)

    muttstub = "echo \"" + message + "\" | mutt"
    if attachment:
        muttstub = muttstub + " -a \"" + attachment + "\""
    cmd = muttstub + " -s \"" + subject + "\" -- " + to
    return runcmd(cmd)


def tar(src):
    dst = src + ".tgz"
    cmd = "tar -cvzf %s %s" % (dst, src)
    runcmd(cmd)
    shutil.rmtree(src)
    return dst


def gzip_all(src):
    cmds = []
    for f in os.listdir(src):
        if f.endswith("gz"): continue
        cmds.append("gzip -f %s" % op.join(src, f))
    runcmd(cmds)


def main(fastq, output, kmernorm, complexity_filter, email, threads=16):
    verbose_logging(**locals())

    fastq_name = op.basename(fastq)
    sample = get_sample(fastq_name)

    try:
        tmpdir = tf.mkdtemp("_tmp", "%s_" % sample, tf.tempdir)
        tmpfastq = op.join(tmpdir, fastq_name)
        shutil.copyfile(fastq, tmpfastq)

        kmerfq = run_kmernorm(tmpfastq, k=19, t=80, c=2)
        spades_input = kmerfq

        if complexity_filter:
            # output is defined in function
            compfilteredfq = lowcomp_filter(**{'paired':None, 'input':kmerfq,
                                                'logfile':LOG, 'min_complexity':0.01})
            spades_input = compfilteredfq
            if op.getsize(spades_input) <= 0:
                logging.info("low complexity filtered has removed all of the reads")
                sys.exit(1)

        spades_dir = tmpdir + "/spades"
        spades_fasta = spades(**{'pe1-12':spades_input, 'threads':threads,
                                 'o':spades_dir, 'sc':None, 'careful':None})

        # moves spades contigs.fasta and adds sample to header name
        renamed_hdrs = postprocess_spades(spades_fasta, sample, tmpdir)

        read_counts(tmpfastq, kmerfq, compfilteredfq)
        assembly_stats(renamed_hdrs)

        sizefilteredfq = filter_fasta_by_size(renamed_hdrs, 2000)
        assembly_stats(sizefilteredfq)

    finally:
        # archive/gzip all of the spades output
        spades_dir = tar(spades_dir)
        # gzip all of the files in the temp dir
        gzip_all(tmpdir)
        # copy over the files
        copyfiles(tmpdir, output, r=None, v=None, h=None, u=None, progress=None)
        # delete the temp working directory
        shutil.rmtree(tmpdir)

        if email:
            send_email(to=email, subject=op.basename(__file__),
                message="finished processing %s; results were copied to %s" % (fastq, output))

        logging.info("Complete.")


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('fastq', help="interweaved PacBio reads")
    p.add_argument('output', help="location to store output files")
    p.add_argument('--kmernorm', action='store_true', help="run kmer normalization prior to spades")
    p.add_argument('--complexity-filter', action='store_true', help="filter out low complexity reads before running spades")
    p.add_argument('--email', default="", help="send completion alert")
    p.add_argument('--threads', default=16, type=int, help="threads for spades to utilize")
    args = p.parse_args()

    if not op.exists(args.output):
        os.makedirs(args.output)

    LOG = "%s/%s.log" % (args.output, op.basename(__file__).rstrip(".py"))

    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',
                        filename=LOG,
                        level=logging.DEBUG,
                        filemode='wb')

    tf.tempdir = tf.gettempdir()
    main(args.fastq, args.output, args.kmernorm, args.complexity_filter, args.email, args.threads)
