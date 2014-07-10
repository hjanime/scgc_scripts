#SCGC scripts
A place to stick random scripts.

##assemble_shuffled.py

Requires:
+ kmernorm
+ mutt
+ spades
+ pip install toolshed

```
usage: assemble_shuffled.py [-h] [--kmernorm] [--complexity-filter]
                            [--email EMAIL] [--threads THREADS]
                            fastq output

assemble interweaved, paired-end reads; the result of
shuffleSequences_fastq.pl.

positional arguments:
  fastq                interweaved, paired-end reads in fastq format
  output               location to store output files

optional arguments:
  -h, --help           show this help message and exit
  --kmernorm           run kmer normalization prior to spades (default: False)
  --complexity-filter  filter out low complexity reads before running spades
                       (default: False)
  --email EMAIL        send completion alert (default: )
  --threads THREADS    threads for spades to utilize (default: 16)
```

##cov_by_chrom.py

Requires:
+ bedtools
+ pip install toolshed

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
