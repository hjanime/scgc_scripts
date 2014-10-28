#!/usr/bin/env python
# coding=utf-8
"""
Lorenz curve.
"""
from __future__ import print_function

import sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from os import path
from toolshed import reader

from rpy2 import robjects
from rpy2.robjects.vectors import IntVector
from rpy2.robjects.packages import importr


grdevices = importr('grDevices')
ineq = importr('ineq')
plot = robjects.r['plot']
bquote = robjects.r['bquote']
legend = robjects.r['legend']


def main(bam, output):
    sample = path.basename(bam).rsplit(".bam", 1)[0]
    plot_file = output if output else bam.rsplit(".bam", 1)[0] + "_lorenz_curve.png"

    coverages = []
    print("Calculating coverages", file=sys.stderr)
    for toks in reader("|bedtools genomecov -5 -d -ibam %s | sort -k3,3n" % bam,
        header=['name', 'start', 'coverage']):
        coverages.append(int(toks['coverage']))

    coverages_r = IntVector(coverages)

    print("Generating Lorenz curve", file=sys.stderr)

    # Gini coefficient
    G = ineq.Gini(coverages_r)
    l = "G = %.3f" % G[0]

    grdevices.png(plot_file, width=1200, height=800)
    # draw the plot
    plot(ineq.Lc(coverages_r), xlab="Genome Fraction",
        ylab="Coverage Fraction", bty="n", lwd=1, main="Lorenz Curve of %s" % sample,
        col="black", xaxs="r", yaxs="r")

    # add the Gini coefficient to the plot
    legend('topleft', legend=l, bty='n', cex=1.3)
    grdevices.dev_off()

    print("Gini Coefficient = %f" % G[0])


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('bam', help="aligned reads")
    p.add_argument('-o', '--output', help="optional output file name (.png)")

    args = p.parse_args()
    main(args.bam, args.output)
