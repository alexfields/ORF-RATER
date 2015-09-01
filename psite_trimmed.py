#! /usr/bin/env python

import argparse
import sys
import pysam
from collections import defaultdict
import numpy as np
from yeti.genomics.genome_array import read_length_nmis
import itertools
import multiprocessing as mp

parser = argparse.ArgumentParser()

parser.add_argument('offsetfile', help='Output file. Two columns of tab-delimited text; first column indicates read length, second column indicates '
                                       'offset to apply. Read lengths are calculated after trimming 5\' mismatches.')
parser.add_argument('bamfiles', nargs='+', help='Path to transcriptome-aligned BAM file(s) for read data')
parser.add_argument('--cdsbed', type=argparse.FileType('rU'), default=sys.stdin,
                    help='BED-file containing annotated CDSs whose start codons are to be used to identify P-site offsets (Default: stdin)')
parser.add_argument('--minrdlen', type=int, default=27, help='Minimum permitted read length (Default: 27)')
parser.add_argument('--maxrdlen', type=int, default=34, help='Maximum permitted read length, inclusive (Default: 34)')
parser.add_argument('--max5mis', type=int, default=1, help='Maximum 5\' mismatches to trim. Reads with more than this number will be excluded. '
                                                           '(Default: 1)')
parser.add_argument('--tallyfile', help='Output file for tallied offsets as a function of read length. First column indicates read length for that '
                                        'row; columns are different offset values, starting at 0.')
parser.add_argument('-p', '--numproc', type=int, default=1, help='Number of processes to run. Defaults to 1 but recommended to use more (e.g. 12-16)')
opts = parser.parse_args()

inbams = [pysam.Samfile(infile) for infile in opts.bamfiles]  # defaults to read mode, and will figure out if it's BAM or SAM - though we require BAM

gcoorddict = defaultdict(set)  # use sets to avoid repeating start codons
for bedline in opts.cdsbed:
    ls = bedline.strip().split()
    chrom = ls[0]
    strand = ls[5]
    thickstart = int(ls[6])
    thickend = int(ls[7])
    if thickstart != thickend:  # coding region is not empty
        if strand == '+':
            gcoorddict[(chrom, strand)].add(thickstart)
        else:
            gcoorddict[(chrom, strand)].add(thickend-1)


def offset_to_gcoord(read, gcoord):
    if read.is_reverse:
        return next((len(read.positions)-idx-1 for (idx, pos) in enumerate(read.positions) if pos == gcoord), None)
    else:
        return next((idx for (idx, pos) in enumerate(read.positions) if pos == gcoord), None)


def get_reads(chrom, strand, gcoord):
    if strand == '+':
        return itertools.ifilter(lambda rd: not rd.is_reverse,
                                 itertools.chain.from_iterable((bamfile.fetch(reference=chrom,
                                                                              start=gcoord,
                                                                              end=gcoord+1) for bamfile in inbams)))
    else:  # strand == '-'
        return itertools.ifilter(lambda rd: rd.is_reverse,
                                 itertools.chain.from_iterable((bamfile.fetch(reference=chrom,
                                                                              start=gcoord,
                                                                              end=gcoord+1) for bamfile in inbams)))


def map_start_sites((chrom, strand)):
    offset_tallies = np.zeros((opts.maxrdlen+1-opts.minrdlen, opts.maxrdlen), dtype=np.uint32)
    for gcoord in gcoorddict[(chrom, strand)]:
        for read in get_reads(chrom, strand, gcoord):
            (rdlen, nmis) = read_length_nmis(read)
            if opts.minrdlen <= rdlen <= opts.maxrdlen:
                curr_offset = offset_to_gcoord(read, gcoord)
                if curr_offset is not None:
                    offset_tallies[rdlen-opts.minrdlen, curr_offset-nmis] += 1
    return offset_tallies

workers = mp.Pool(opts.numproc)
offset_tallies = sum(workers.map(map_start_sites, gcoorddict.keys()))
workers.close()

if opts.tallyfile:
    with open(opts.tallyfile, 'w') as outfile:
        for rdlen in xrange(opts.minrdlen, opts.maxrdlen+1):
            outfile.write('\t'.join([str(rdlen)]+[str(tally) for tally in offset_tallies[rdlen-opts.minrdlen, :]])+'\n')
with open(opts.offsetfile, 'w') as outfile:
    for rdlen in xrange(opts.minrdlen, opts.maxrdlen+1):
        outfile.write('%d\t%d\n' % (rdlen, np.argmax(offset_tallies[rdlen-opts.minrdlen, :])))
for inbam in inbams:
    inbam.close()
