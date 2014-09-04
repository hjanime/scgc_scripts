#SCGC scripts
A place to stick random scripts. Most will require `toolshed` which can be
installed via `pip`.

##bcl2fastq.py
```
usage: bcl2fastq.py [-h] [-l {NONE,FATAL,ERROR,WARNING,INFO,DEBUG,TRACE}]
                    [-r LOADING_THREADS] [-d DEMULTIPLEXING_THREADS]
                    [-p PROCESSING_THREADS]
                    [--barcode-mismatches BARCODE_MISMATCHES]
                    runfolder ...

Runs bcl2fastq creating fastqs and concatenates fastqs across lanes. Intended
to be used with NextSeq data and it does not do any cleanup! Original dumped
fastqs will remain along with all of the bcl files.

positional arguments:
  runfolder             path to run folder
  args                  any additional bcl2fastq args and their values

optional arguments:
  -h, --help            show this help message and exit
  -l {NONE,FATAL,ERROR,WARNING,INFO,DEBUG,TRACE}, --min-log-level {NONE,FATAL,ERROR,WARNING,INFO,DEBUG,TRACE}
                        minimum log level (default: INFO)
  -r LOADING_THREADS, --loading-threads LOADING_THREADS
                        threads used for loading BCL data (default: 12)
  -d DEMULTIPLEXING_THREADS, --demultiplexing-threads DEMULTIPLEXING_THREADS
                        threads used for demultiplexing (default: 12)
  -p PROCESSING_THREADS, --processing-threads PROCESSING_THREADS
                        threads used for processing demultiplexed data
                        (default: 12)
  --barcode-mismatches BARCODE_MISMATCHES
                        number of allowed mismatches per index (default: 0)
```

##cov_by_chrom.py

Requires:
+ bedtools

```
usage: cov_by_chrom.py [-h] [--no-split] bam

avg. per chrom coverage for a bam file. defaults to using:

bedtools genomecov -d -split -ibam <bam>

positional arguments:
  bam         coordinate sorted bam file

optional arguments:
  -h, --help  show this help message and exit
  --no-split  do not split reads based on CIGAR
```

##fastx.py

Requires: `pip install click parmap`

```
$ fastx.py --help
Usage: fastx.py [OPTIONS] COMMAND [ARGS]...

  Fasta and fastq tools.

Options:
  --version  Show the version and exit.
  --help     Show this message and exit.

Commands:
  assembly-stats     basic fasta stats
  complexity-filter  fastx complexity filter
  count              count the reads
  header             prepend or append to header
  length-filter      filter fasta by seq length
  merge-pe           merge R1 and R2 into interleaved fastq
  sliding-gc         sliding GC and skew calculations
  split-merged       unmerge interweaved fastq file
```
