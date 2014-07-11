#!/usr/bin/env python
# coding=utf-8
"""
1.  prodigal
2.  tRNAscan-SE
3.  calculate gc content and skew
4.  tetramerPCA (/opt/scgc/rscript/pipeline-tetramerPCA.Rscript)
5.  blastp
6.  blastx viral
7.  blastx bacterial
8.  coverage
9.  mpileup

It's recommended that you set environmental variable TMPDIR to something local,
eg. export TMPDIR=/tmp
"""

import fnmatch
import logging
import os
import os.path as op
import parmap
import shutil
import stat
import subprocess as sp
import sys
import tempfile as tf
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from collections import Counter
from Bio.Blast import NCBIXML
from itertools import groupby
from toolshed import nopen, reader


def verbose_logging(**kwargs):
    logging.info("Record of arguments for %s" % os.path.basename(__file__))
    for k, v in kwargs.items():
        logging.info(">>> %s = %s", k, v)


def find_files(path, pattern):
    """returns all files (recursively) under path that match pattern or
    patterns.
    """
    found_files = []
    if isinstance(pattern, basestring):
        for root, dirs, files in os.walk(path):
            for filename in fnmatch.filter(files, pattern):
                found_files.append(op.join(root, filename))
    else:
        assert isinstance(pattern, list)
        for root, dirs, files in os.walk(path):
            for p in pattern:
                for filename in fnmatch.filter(files, p):
                    found_files.append(op.join(root, filename))
    return found_files


def copytree(src, dst, symlinks=False, ignore=None):
    if not os.path.exists(dst):
        os.makedirs(dst)
        shutil.copystat(src, dst)
    lst = os.listdir(src)
    if ignore:
        excl = ignore(src, lst)
        lst = [x for x in lst if x not in excl]
    for item in lst:
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if symlinks and os.path.islink(s):
            if os.path.lexists(d):
                os.remove(d)
            os.symlink(os.readlink(s), d)
            try:
                st = os.lstat(s)
                mode = stat.S_IMODE(st.st_mode)
                os.lchmod(d, mode)
            except:
                pass # lchmod not available
        elif os.path.isdir(s):
            copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)


def get_sample(filename):
    name, ext = op.splitext(filename)
    if not ext == ".fasta" and not ext == ".fa":
        name, ext = op.splitext(name)
        if not ext == ".fasta" and not ext == ".fa":
            logging.error("Looking for .fasta file; found %s.%s", name, ext)
            raise IOError
    return name


def preprocess_fasta(fasta, tmpdir):
    logging.debug("preprocessing: %s", fasta)
    tmpfasta = op.join(tmpdir, op.basename(fasta))
    with open(tmpfasta, 'wb') as fh:
        for name, seq in readfa(fasta):
            print >>fh, ">%s\n%s" % (name, seq)
    return tmpfasta


def runcmd(cmd, log=True):
    """
    >>> import sys
    >>> import logging
    >>> import subprocess as sp
    >>> logging.basicConfig(stream=sys.stdout)
    >>> runcmd("date > /dev/null", False)
    >>> runcmd(["date > /dev/null", "date > /dev/null"], False)
    """
    if isinstance(cmd, basestring):
        if log: logging.debug("Running command: %s", cmd)
        # will raise sp.CalledProcessError on retcode != 0
        sp.check_call(cmd, shell=True)

    else:
        assert isinstance(cmd, list)
        if log: logging.debug("Running commands: %s", cmd)

        # this will be bad if the list of commands is really long
        processes = [sp.Popen(c, shell=True) for c in cmd]
        for i, p in enumerate(processes):
            p.wait()
            if p.returncode != 0:
                raise sp.CalledProcessError(p.returncode, cmd[i])


def prodigal(fasta, sample, outdir):
    """Runs prodigal to call genes and returns 4 file names; proteins, genes,
    genbank, and score. Simpler to leave this one as is and not change it like
    tRNAscan due to having to track output files.

    @param fasta - input fasta file
    @param sample - identifier for this sample or run
    @param outdir - where to place output files
    """
    proteins = op.join(outdir, sample + ".proteins.fasta")
    genes = op.join(outdir, sample + ".genes.fasta")
    genbank = op.join(outdir, sample + ".gbk")
    score = op.join(outdir, sample + ".scores")
    cmd = "prodigal -a {proteins} -d {genes} -i {fasta} -o {genbank} -p meta -s {score}".format(**locals())
    runcmd(cmd)
    return proteins, genes, genbank, score


def readfa(fa):
    """Fasta iterator."""
    with nopen(fa) as fh:
        for header, group in groupby(fh, lambda line: line[0] == '>'):
            if header:
                line = group.next()
                name = line[1:].strip().split(" ", 1)[0]
            else:
                seq = ''.join(line.strip() for line in group)
                yield name, seq


def fa_seqs(fasta):
    """count num of sequences"""
    total = 0
    for name, seq in readfa(fasta):
        total += 1
    return total


def kwargs_to_flag_string(kwargs, ignore=[], dash="-"):
    """
    >>> kwargs = {"d":"d_test", "p":"p_test", "sep":","}
    >>> kwargs_to_flag_string(kwargs, ignore=['sep'], dash="-")
    ' -p p_test -d d_test'
    >>> kwargs_to_flag_string(kwargs, dash="--")
    ' --p p_test --d d_test --sep ,'
    """
    s = ""
    for k, v in kwargs.items():
        if k in ignore: continue
        s += " " + dash + k + ("" if v is None else (" " + str(v)))
    return s


def trnascan(fasta, **kwargs):
    """
    Run tRNAscan-SE and return the file name. Q kwarg is always used.

    @param fasta - contigs to run
    """

    def _rename_fa_headers(fa):
        """renames the headers and returns tmp fasta and dictionary."""
        # copy incoming fasta over to tmp
        tmp_copied_fa = "%s/%s" % (tf.tempdir, os.path.basename(fa))
        shutil.copyfile(fa, tmp_copied_fa)

        try:
            d = {}
            tmp_fa = open(tf.mkstemp(suffix=".fasta")[1], 'wb')
            for i, (name, seq) in enumerate(readfa(tmp_copied_fa)):
                d[str(i)] = name
                print >> tmp_fa, ">%d\n%s" % (i, seq)
            tmp_fa.close()

        finally:
            # remove the unnamed fasta copied over to tmp
            os.remove(tmp_copied_fa)

        return tmp_fa.name, d

    if kwargs.has_key('o'):
        out = kwargs['o']
    else:
        # user did not specify output; give temp in curdir
        out = tf.mkstemp(suffix=".trnascan", dir=".")[1]

    # need to redirect uncorrected output
    kwargs['o'] = tf.mkstemp(suffix=".trnascan")[1]
    # never attempt to prompt if user wants overwrite existing file
    kwargs['Q'] = None

    # read the fasta; fasta gets renamed to a temp file name
    fasta, fasta_dict = _rename_fa_headers(fasta)

    try:
        cmd = "tRNAscan-SE" + kwargs_to_flag_string(kwargs) + " " + fasta
        runcmd(cmd)

        # rename first column of tRNAscan output back to full header name
        with open(out, 'wb') as o:
            for i, line in enumerate(open(kwargs['o'], 'rb')):
                # first 3 lines of tRNAscan-SE comprise the header
                line = line.rstrip("\r\n")
                if i < 3:
                    print >>o, line
                else:
                    line = line.rstrip('\r\n').split()
                    # KeyError possibility; *should* never happen
                    line[0] = fasta_dict[line[0]]
                    print >>o, "\t".join(line)

    finally:
        os.remove(fasta)
        os.remove(kwargs['o'])

    return out


def gc_skew_and_content(seq, window_size=500):
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
    half_window = window_size / 2
    seq = seq.upper()
    seq_len = len(seq)
    assert seq_len >= window_size

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


def gc_content(fasta):
    """Calculate GC content and skew, write to file, and return output file
    name.
    """
    output = fasta.rsplit(".", 1)[0] + ".gc_content.tsv"
    header = ["SEQUENCE_NAME", "POSITION", "SKEW", "CONTENT"]

    with open(output, 'wb') as fout:
        print >>fout, "\t".join(header)
        for name, seq in readfa(fasta):
            # preserve contig order in the calculation
            for point, skew, content in gc_skew_and_content(seq, 500):
                print >>fout, "%s\t%i\t%0.3f\t%0.3f" % (name, point, skew, content)

    return output


def tetramerPCA(fasta, **kwargs):
    output_dir = op.dirname(fasta) + "/tetramerPCA"
    kwargs['input'] = fasta
    kwargs['output_dir'] = output_dir
    cmd = "Rscript --vanilla /opt/scgc/rscript/pipeline-tetramerPCA.Rscript" + \
                kwargs_to_flag_string(kwargs, dash="--")
    runcmd(cmd)


def blastp(**kwargs):
    assert kwargs.has_key('out')
    cmd = "parallel-blastp" + kwargs_to_flag_string(kwargs)
    runcmd(cmd)


def makeblastdb(**kwargs):
    # TODO: make a generic command and remove all these other ones
    # runcmd seems like the ideal place, but i also like it uncoupled
    cmd = "makeblastdb" + kwargs_to_flag_string(kwargs)
    runcmd(cmd)


def batch_blastx(blast_db, viral_db, tmpdir, evalue=0.001, threads=20, outfmt=5):
    """returns list of result files and list of used viral databases"""
    # could/should be a parameter if this function leaves this script
    valid_db_exts = ["*.fna.gz", "*.fasta.gz", "*.fa.gz", "*.fa", "*.fna"]
    res = []
    queries = find_files(viral_db, valid_db_exts)

    # run blastx for each query
    for query in queries:
        xml_output = op.join(tmpdir, "%s.viral_fraction.xml" % op.basename(query).split(".", 1)[0])
        res.append(xml_output)
        if not op.exists(xml_output):
            blastx(query=query, db=blast_db, evalue=evalue, max_target_seqs=1,
                    num_threads=threads, outfmt=5, query_gencode=11, out=xml_output)
    return res, queries


def blastx(**kwargs):
    # TODO rewrite parallel-blast
    cmd = "parallel-blastx" + kwargs_to_flag_string(kwargs)
    runcmd(cmd)


def construct_cigar(query_seq, subj_seq, query_len, query_from, query_to):
    """given two blast alignments, construct a valid CIGAR string for SAM
    http://samtools.github.io/hts-specs/SAMv1.pdf
    """
    cigar_string = []
    last_op = None
    last_op_length = 0
    assert query_from < query_to

    # hard clipping prior to alignment
    if query_from > 1:
        cigar_string.append(str(query_from - 1) + "H")

    for qseq, sseq in zip(query_seq, subj_seq):

        # query sequence has a gap; deletion compared to reference
        if qseq == '-':
            current_op = 'D'

        # subject sequence has a gap; insertion compared to reference
        elif sseq == '-':
            current_op = 'I'

        # sequence match or mismatch
        else:
            current_op = 'M'

        # change in operation
        if current_op != last_op:
            if last_op:
                cigar_string.append(str(last_op_length) + last_op)
            last_op = current_op
            last_op_length = 1

        else:
            last_op_length += 1

    cigar_string.append(str(last_op_length) + last_op)

    # hard clipping after alignment
    if query_to < query_len:
        cigar_string.append(str(query_len - query_to) + "H")

    return cigar_string


def xml_to_bam(blast_xml, query_fasta, reference_fasta, threads=12):
    """
    convert protein sequence XML to protein sequence bam.
    """
    filepath = blast_xml.rsplit(".", 1)[0]
    sam_file = filepath + ".sam"
    bam_file = filepath + ".bam"
    faidx = reference_fasta + ".fai"

    try:
        with open(sam_file, 'wb') as sam, nopen(blast_xml) as xml:

            for record in NCBIXML.parse(nopen(xml)):
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        # still a list for future reversing if anything ever
                        # maps to negative strand; assertion in function
                        cigar = construct_cigar(hsp.query, hsp.sbjct,
                                        record.query_letters, hsp.query_start,
                                        hsp.query_end)
                        # 11 fields
                        sam_fields = [record.query.split(" ", 1)[0],
                                        0,
                                        alignment.title.split(" ", 1)[0],
                                        hsp.sbjct_start,
                                        255,
                                        "".join(cigar),
                                        '*',
                                        0,
                                        0,
                                        hsp.query.replace("-", ""),
                                        '*']
                        print >>sam, "\t".join(map(str, sam_fields))

        # index the reference fasta
        # there should be a step here to gunzip this file
        runcmd("samtools faidx %s" % reference_fasta)

        # convert sam to bam
        cmd = ("cat {sam} | samtools view -S - -b -t {fasta} | "
                "samtools sort -@ {threads} -m 1G - {out}").format(sam=sam_file,
                    fasta=faidx, threads=threads, out=bam_file.split(".bam")[0])
        runcmd(cmd)

        # index the bam
        cmd = "samtools index " + bam_file
        runcmd(cmd)

    finally:
        if op.exists(sam_file): os.remove(sam_file)
        if op.exists(faidx): os.remove(faidx)

    return bam_file


def contig_coverage_ratio(bam, min_cov=1):
    """
    return coverage ratio dictionary of:

        (bases covered by at least min_cov / contig length)

    key is fasta name prior to first space.
    """
    coverage_ratios = {}
    cmd = ("| samtools view -u {bam} "
           "| bedtools genomecov -d -split -ibam - "
           "| sort -k1,1 -k2,2n").format(bam=bam)

    for name, igroup in groupby(reader(cmd, header=['name', 'pos', 'count']), lambda x: x['name']):
        # not sure if necessary
        if name == 'genome': continue
        covered = 0
        total = 0.0
        for toks in igroup:
            total += 1
            count = int(toks['count'])
            if count >= min_cov:
                covered += 1
        coverage_ratios[toks['name']] = covered / total

    return coverage_ratios


def per_contig_coverage(bam, fasta, min_cov=1):
    """
    given the reference fasta and the bam, calculate average coverage for each
    contig and return new file name that is placed in same dir as bam.
    """
    out_hdr = ['CONTIG', 'LENGTH', 'AVERAGE_COVERAGE', 'Coverage_ratio_at_%sx' % min_cov]
    out = bam.rsplit(".bam", 1)[0] + ".coverage.tsv"
    coverage_ratios = contig_coverage_ratio(bam)

    total_reads = int(nopen("| samtools view -c -F4 " + bam).next().strip("\r\n"))
    if total_reads == 0: return out

    # consensus of contigs that are present
    contig_sizes = {}
    fasta_contigs = set()
    for name, seq in readfa(fasta):
        fasta_contigs.add(name)
        contig_sizes[name] = len(seq)

    # contigs with coverage and their coverage
    coverages = Counter()
    genomecov_hdr = ['name', 'coverage', 'bases_at_coverage', 'total_bases', 'fraction']
    for toks in reader("| bedtools genomecov -ibam %s" % bam, header=genomecov_hdr):
        coverages[toks['name']] += int(toks['coverage']) * int(toks['bases_at_coverage'])

    coverages_set = set(coverages.keys())

    with open(out, 'wb') as fh:
        print >>fh, "\t".join(out_hdr)
        # should not need to iterate over this again;
        # add groupby in previous loop
        for name, cov in coverages.iteritems():
            if name == "genome": continue
            avg_cov = cov / float(contig_sizes[name])
            cov_ratio = coverage_ratios[name]
            print >>fh, "%s\t%d\t%0.3f\t%0.3f" % (name, contig_sizes[name], avg_cov, cov_ratio)

        # nonzero coverage; not in bedtools genomecov output
        for name in fasta_contigs - coverages_set:
            print >>fh, "%s\t%d\t0\t0" % (name, contig_sizes[name])

    return out


def samtools_mpileup(bam):
    out = bam.rsplit(".bam", 1)[0] + ".mpileup.tsv"
    runcmd("samtools mpileup %s | cut -f 1,2,4 > %s" % (bam, out))
    return out


def gzip_all(src, ignore=['.gz']):
    if not '.gz' in ignore: ignore.append('.gz')
    cmds = []

    for f in os.listdir(src):
        f = op.join(src, f)
        if op.isdir(f):
            gzip_all(f, ignore=ignore)
            continue

        append = True
        for extension in ignore:
            if f.endswith(extension):
                append = False

        if append:
            cmds.append("gzip -f %s" % f)

    if len(cmds) > 0:
        runcmd(cmds)


def send_email(to="scgc@bigelow.org", subject="", message="", attachment=None):
    if isinstance(message, list):
        message = " ".join(message)

    muttstub = "echo \"" + message + "\" | mutt"
    if attachment:
        muttstub = muttstub + " -a \"" + attachment + "\""
    cmd = muttstub + " -s \"" + subject + "\" -- " + to
    return runcmd(cmd)


def main(fasta, output_dir, bacterial_query, email, viral_db,
            keep_tmp=False, threads=20, evalue=0.001, outfmt=5):

    verbose_logging(**locals())

    fasta_name = op.basename(fasta)
    sample = get_sample(fasta_name)

    try:
        tmpdir = tf.mkdtemp("_tmp", "%s_" % sample, tf.tempdir)
        # rename the fa headers while writing to temp working dir
        tmpfasta = preprocess_fasta(fasta, tmpdir)

        blastp_xml = op.join(tmpdir, "%s.blastp.xml" % sample)
        blast_db = op.join(tmpdir, "%s.blastdb" % sample)

        tmp_viral_db = op.join(tmpdir, op.basename(viral_db))
        copytree(viral_db, tmp_viral_db)

        tmpquery = op.join(tmpdir, op.basename(bacterial_query))
        shutil.copyfile(bacterial_query, tmpquery)

        # prodigal
        p_proteins, p_genes, p_genbank, p_score = prodigal(tmpfasta, sample, tmpdir)

        # tRNAscan-SE
        trna_output = os.path.join(tmpdir, "%s.trnascan" % sample)
        trna_output = trnascan(tmpfasta, o=trna_output, B=None)

        # gc content and skew
        gc_output = gc_content(tmpfasta)

        # tetramerPCA
        tetramerPCA_output = tetramerPCA(tmpfasta, window=1600, step=200)

        # blastp
        blastp(query=p_proteins, db='nr', out=blastp_xml, evalue=evalue,
                num_alignments=10, num_threads=threads, outfmt=outfmt)

        # blast_db
        makeblastdb(**{'in':p_proteins, 'parse_seqids':None, 'dbtype':'prot', 'out':blast_db})

        # blastx; viral
        viral_xmls, viral_queries = batch_blastx(blast_db, tmp_viral_db, tmpdir, evalue, threads, outfmt)
        viral_bams = parmap.starmap(xml_to_bam, zip(viral_xmls, viral_queries), p_proteins, processes=12)
        viral_coverages = parmap.map(per_contig_coverage, viral_bams, p_proteins, processes=12)
        viral_pileups = parmap.map(samtools_mpileup, viral_bams, processes=12)

        # blastx; bacterial
        bacterial_xml = os.path.join(tmpdir, "%s.bacterial_fraction.xml" % sample)
        blastx(query=tmpquery, db=blast_db, out=bacterial_xml, evalue=evalue,
                max_target_seqs=1, num_threads=threads, outfmt=5, query_gencode=11)
        bacterial_bam = xml_to_bam(bacterial_xml, bacterial_query, p_proteins)
        bacterial_coverage = per_contig_coverage(bacterial_bam, p_proteins)
        bacterial_pileup = samtools_mpileup(bacterial_bam)

    except Exception as e:
        print e.__doc__
        print e.message
        raise

    finally:
        # always remove viral fastas; don't copy back from ram
        if op.exists(tmp_viral_db):
            shutil.rmtree(tmp_viral_db)
        # gzip all of the files in the temp dir
        gzip_all(tmpdir, ignore=['pdf', 'bam', 'bai'])
        # copy over the files
        runcmd("cp -R -v {src}/* {dst}".format(src=tmpdir, dst=output_dir))

        if not keep_tmp:
            # delete the temp working directory
            shutil.rmtree(tmpdir)

        if email:
            send_email(to=email, subject=op.basename(__file__),
                message="finished processing %s; results were copied to %s" % (fasta, output_dir))

        logging.info("Complete.")


if __name__=='__main__':
    import doctest
    p = ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
    p.add_argument('fasta', metavar="FASTA", help="sequences to analyze")
    p.add_argument('output', metavar="OUTPUT", help="location for output files")
    p.add_argument('query', metavar="QUERY", help="bacterial query fasta")

    p.add_argument("--email", default='',
            help="send an email to this address when complete (default: <skip>)")
    p.add_argument("--keep-tmp", action="store_true",
            help="preserve temp files on exist (default: %(default)s)")
    p.add_argument("--threads", default=90, type=int,
            help="max number of threads to use (default: %(default)s)")
    p.add_argument("--viral-db", default="/mnt/scgc/viral_dbs",
            help="path to viral databases (default: %(default)s)")
    args = p.parse_args()

    if args.output == ".":
        args.output = op.abspath(os.curdir)
    elif not op.exists(args.output):
        os.makedirs(args.output)

    # TODO parameterize log level
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',
                        filename="%s/%s.log" % (args.output, op.basename(__file__).rstrip(".py")),
                        level=logging.DEBUG,
                        filemode='wb')

    if not op.exists(args.fasta):
        logging.error("Reads file does not exist: %s", args.fasta)
        sys.exit(1)

    if not op.exists(args.query):
        logging.error("Bacterial query does not exist: %s", args.query)
        sys.exit(1)

    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        tf.tempdir = tf.gettempdir()
        main(args.fasta, args.output, args.query, args.email, args.viral_db, \
                args.keep_tmp, args.threads)
