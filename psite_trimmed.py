#! /usr/bin/env python

import argparse
import sys
import os
import pysam
from collections import defaultdict
import numpy as np
from hashed_read_genome_array import read_length_nmis
import itertools
import multiprocessing as mp
from time import strftime

parser = argparse.ArgumentParser(description='Find most common P-site offset for each read length in a ribosome profiling experiment. If multiple '
                                             'ribosome profiling datasets are to be analyzed separately (e.g. if they were collected under different '
                                             'drug treatments), then this program should be run separately for each, ideally in separate subfolders '
                                             'indicated by SUBDIR.')

parser.add_argument('bamfiles', nargs='+', help='Path to transcriptome-aligned BAM file(s) for read data')
parser.add_argument('--subdir', default=os.path.curdir,
                    help='Convenience argument when dealing with multiple datasets. In such a case, set SUBDIR to an appropriate name (e.g. HARR, '
                         'CHX) to avoid file conflicts. (Default: current directory)')
parser.add_argument('--offsetfile', default='offsets.txt',
                    help='Output file. Two columns of tab-delimited text; first column indicates read length, second column indicates offset to '
                         'apply. Read lengths are calculated after trimming 5\' mismatches. If SUBDIR is set, this file will be placed in that '
                         'directory. (Default: offsets.txt)')
parser.add_argument('--cdsbed', type=argparse.FileType('rU'), default=sys.stdin,
                    help='BED-file containing annotated CDSs whose start codons are to be used to identify P-site offsets (Default: stdin)')
parser.add_argument('--minrdlen', type=int, default=27, help='Minimum permitted read length, inclusive (Default: 27)')
parser.add_argument('--maxrdlen', type=int, default=34, help='Maximum permitted read length, inclusive (Default: 34)')
parser.add_argument('--max5mis', type=int, default=1, help='Maximum 5\' mismatches to trim. Reads with more than this number will be excluded. '
                                                           '(Default: 1)')
parser.add_argument('--tallyfile', help='Optional output file for tallied offsets as a function of read length. First column indicates read length '
                                        'for that row; columns are different offset values, starting at 0. Will be placed in SUBDIR automatically.')
parser.add_argument('-v', '--verbose', action='store_true', help='Output a log of progress and timing (to stdout)')
parser.add_argument('-p', '--numproc', type=int, default=1, help='Number of processes to run. Defaults to 1 but more recommended if available.')
parser.add_argument('-f', '--force', action='store_true', help='Force file overwrite')
opts = parser.parse_args()

outfilename = os.path.join(opts.subdir, opts.offsetfile)
if not opts.force and os.path.exists(outfilename):
    raise IOError('%s exists; use --force to overwrite' % outfilename)

if opts.tallyfile:
    tallyfilename = os.path.join(opts.subdir, opts.tallyfile)
    if not opts.force and os.path.exists(tallyfilename):
        raise IOError('%s exists; use --force to overwrite' % tallyfilename)

if not os.path.isdir(opts.subdir):
    os.mkdir(opts.subdir)

if opts.verbose:
    sys.stdout.write(' '.join(sys.argv) + '\n')

    def logprint(nextstr):
        sys.stdout.write('[%s] %s\n' % (strftime('%Y-%m-%d %H:%M:%S'), nextstr))
        sys.stdout.flush()

    logprint('Identifying reads near annotated translation start sites')

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


def _offset_to_gcoord(read, gcoord):
    """Identify the offset from the 5' end of a pysam.AlignedRead corresponding to a specified coordinate"""
    if read.is_reverse:
        return next((len(read.positions)-idx-1 for (idx, pos) in enumerate(read.positions) if pos == gcoord), None)
    else:
        return next((idx for (idx, pos) in enumerate(read.positions) if pos == gcoord), None)


def _get_reads(chrom, strand, gcoord, inbams):
    """Get all of the reads overlapping a specific genomic position"""
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


def _map_start_sites((chrom, strand)):
    """Tally reads by read length and offset from translation start sites, for a particular chromsome and strand. Each read contributes the same
    amount, so a start position with more reads will contribute more than one with fewer"""

    inbams = [pysam.Samfile(infile, 'rb') for infile in opts.bamfiles]
    offset_tallies = np.zeros((opts.maxrdlen+1-opts.minrdlen, opts.maxrdlen), np.uint32)
    for gcoord in gcoorddict[(chrom, strand)]:
        for read in _get_reads(chrom, strand, gcoord, inbams):
            (rdlen, nmis) = read_length_nmis(read)
            if opts.minrdlen <= rdlen <= opts.maxrdlen and nmis <= opts.max5mis:
                curr_offset = _offset_to_gcoord(read, gcoord)
                if curr_offset is not None and curr_offset >= nmis:
                    offset_tallies[rdlen-opts.minrdlen, curr_offset-nmis] += 1
    for inbam in inbams:
        inbam.close()
    return offset_tallies

workers = mp.Pool(opts.numproc)
offset_tallies = sum(workers.map(_map_start_sites, gcoorddict.keys())).astype(np.float64)  # convert to float, or at least int64, for convolution
workers.close()

if opts.verbose:
    logprint('Saving results')

if opts.tallyfile:
    with open(tallyfilename, 'w') as outfile:
        for rdlen in xrange(opts.minrdlen, opts.maxrdlen+1):
            outfile.write('\t'.join([str(rdlen)]+[str(tally) for tally in offset_tallies[rdlen-opts.minrdlen, :]])+'\n')

# Strategy for determining P-site offsets: find the P-site for the most common read length the simple way (argmax), then see how much each other
# read length's profile is shifted relative to that. The naive argmax approach can result in some read lengths shifted by e.g. 3 extra nucleotides;
# this is more robust.
master_rdlen_idx = np.argmax(offset_tallies.sum(1))  # this is the row with the most reads, to be used as the master
master_offset = np.argmax(offset_tallies[master_rdlen_idx, :])
offsets = [master_offset+np.argmax(np.correlate(row, offset_tallies[master_rdlen_idx, :], 'full'))+1-opts.maxrdlen for row in offset_tallies]

with open(outfilename, 'w') as outfile:
    for (i, offset) in enumerate(offsets):
        outfile.write('%d\t%d\n' % (opts.minrdlen+i, offset))

if opts.verbose:
    logprint('Tasks complete')
