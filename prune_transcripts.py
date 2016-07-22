#! /usr/bin/env python

import argparse
import os
import sys
from plastid.genomics.genome_array import BAMGenomeArray, FivePrimeMapFactory, SizeFilterFactory
from plastid.genomics.roitools import SegmentChain
import pysam
from collections import defaultdict
from Bio import SeqIO
from string import maketrans
import pandas as pd
import numpy as np
import multiprocessing as mp
from time import strftime

parser = argparse.ArgumentParser(description='Use ribosome profiling data to remove unwanted transcripts from a transcriptome. Transcripts will be '
                                             'removed if they do not have enough reads of the desired length, if they have greater than some '
                                             'fraction of their reads from a single position, if they are annotated as pseudogenes (and multimap '
                                             'with another transcript), or if they have more multimapping reads beyond what would be expected based '
                                             'on the number of multimapping positions. It is recommended that this file be run in an empty directory '
                                             'and that OUTBED remain at the default value ("transcripts.bed") for consistency with later scripts.')
parser.add_argument('--inbed', type=argparse.FileType('rU'), default=sys.stdin, help='Transcriptome BED-file (Default: stdin)')
parser.add_argument('genomefasta', help='Path to genome FASTA-file')
parser.add_argument('bamfiles', nargs='+', help='Path to transcriptome-aligned BAM file(s) to use for transcript filtering purposes. Alignment '
                                                'should report all possible multimapping positions for each read. Ideally, should be ribosome '
                                                'profiling data sets collected in the absence of initiation inhibitors (e.g. CHX or no drug).')
parser.add_argument('--summarytable',
                    help='Filename to use for (optional) tab-delimited text output (including column titles). First column is transcript IDs, '
                         'followed by summary information such as number of mappable positions and reads, maximum reads from any one position, '
                         'and why the transcript was dropped (if it was dropped).')
parser.add_argument('--outbed', default='transcripts.bed',
                    help='File to which to output BED-formatted transcripts that passed all filters (Default: transcripts.bed)')
parser.add_argument('--minlen', type=int, default=29,
                    help='Minimum length of read to be considered when evaluating transcripts. '
                         'Also serves as the size of the segment when identifying multimapping positions. Must be <= 31 (Default: 29)')
parser.add_argument('--maxlen', type=int, default=30,
                    help='Maximum length (inclusive) of read to be considered when evaluating transcripts. Must be >= MINLEN (Default: 30)')
parser.add_argument('--minreads', type=int, default=64, help='Minimum number of reads demanded for each transcript (Default: 64)')
parser.add_argument('--peakfrac', type=float, default=1./5,
                    help='Maximum fraction of a transcript\'s reads coming from any one position (Default: 0.2)')
parser.add_argument('--pseudogenes', help='List of transcript IDs annotated as pseudogenes, one per line. To be removed if >PSEUDOFRAC fraction of '
                                          'reads mapping to each are multimappers.')
parser.add_argument('--pseudofrac', type=float, default=1./3, help='Maximum allowable multimapping fraction of an annotated pseudogene\'s reads. '
                                                                   'Ignored if list of pseudogenes is not provided. (Default: 0.333)')
parser.add_argument('--multiexcess', type=float, default=1./3,
                    help='Maximum disparity in multimapping reads versus multimapping positions for any transcript (Default: 0.333)')
parser.add_argument('--keeptempfiles', action='store_true', help='Keep the generated intermediate files (useful for debugging)')
parser.add_argument('-v', '--verbose', action='count', help='Output a log of progress and timing (to stdout). Repeat for higher verbosity level.')
parser.add_argument('-p', '--numproc', type=int, default=1, help='Number of processes to run. Defaults to 1 but more recommended if available.')
parser.add_argument('-f', '--force', action='store_true', help='Force file overwrite')
opts = parser.parse_args()

if not opts.force and os.path.exists(opts.outbed):
    raise IOError('%s exists; use --force to overwrite' % opts.outbed)

if opts.minlen > 31:
    raise ValueError('MINLEN must be <= 31 (currently %d)' % opts.minlen)

if opts.minlen > opts.maxlen:
    raise ValueError('MINLEN must be <= MAXLEN (currently %d and %d, respectively)' % (opts.minlen, opts.maxlen))

if opts.verbose:
    sys.stdout.write(' '.join(sys.argv) + '\n')

    def logprint(nextstr):
        sys.stdout.write('[%s] %s\n' % (strftime('%Y-%m-%d %H:%M:%S'), nextstr))
        sys.stdout.flush()

    logprint('Reading transcriptome and genome')
    log_lock = mp.Lock()

fpsize = opts.minlen
psite = (fpsize + 1) / 2

if opts.pseudogenes:
    with open(opts.pseudogenes, 'rU') as inpseudo:
        pseudotids = {line.strip() for line in inpseudo}
else:
    pseudotids = {}

# Parse through all transcripts once to hash them by (chrom,strand) for easy reference later
bedlinedict = defaultdict(dict)
ntids = 0
for line in opts.inbed:
    ls = line.split()
    bedlinedict[(ls[0], ls[5])][ls[3]] = line
    ntids += 1

if ntids == 0:
    raise EOFError('Insufficient input or empty file provided')

genome = SeqIO.to_dict(SeqIO.parse(opts.genomefasta, 'fasta'))
str_dict = maketrans('ACGTacgt', '01230123')

temp_folder = 'tid_seq_info_temp'
if os.path.exists(temp_folder):
    num = 1
    while os.path.exists('tid_seq_info_temp_%d' % num):
        num += 1
    temp_folder = 'tid_seq_info_temp_%d' % num
os.mkdir(temp_folder)

# seq_info_hdf_orig = os.path.join(temp_folder, 'tid_seq_%s%s_orig.h5')
seq_info_hdf = os.path.join(temp_folder, 'tid_seq_%s%s.h5')


def _get_tid_info(tup):
    """For each transcript on this chromosome/strand, identifies every sub-sequence of the appropriate length (fpsize), converts it to an integer,
    identifies the number of reads mapping to that position, and outputs all of that information to a pandas HDF store."""
    (chrom, strand) = tup
    inbams = [pysam.Samfile(infile, 'rb') for infile in opts.bamfiles]
    gnd = BAMGenomeArray(inbams, mapping=FivePrimeMapFactory(psite))
    # map to roughly the center of each read so that identical sequences that cross different splice sites
    # (on different transcripts) still end up mapping to the same place
    gnd.add_filter('size', SizeFilterFactory(opts.minlen, opts.maxlen))

    tid_seq_info = []
    tid_summary = pd.DataFrame(
        {'chrom': chrom, 'strand': strand, 'n_psite': -1, 'n_reads': -1, 'peak_reads': -1, 'dropped': ''},
        index=pd.Index(bedlinedict[(chrom, strand)].keys(), name='tid'))
    for (tid, line) in bedlinedict[(chrom, strand)].iteritems():
        currtrans = SegmentChain.from_bed(line)
        curr_pos_list = currtrans.get_position_list()  # not in stranded order!
        if strand == '-':
            curr_pos_list = curr_pos_list[::-1]
        n_psite = len(curr_pos_list) + 1 - fpsize
        tid_summary.at[tid, 'n_psite'] = n_psite
        if n_psite > 0:
            curr_counts = np.array(currtrans.get_counts(gnd))[psite:n_psite + psite]
            #                if((curr_counts>0).any()):
            sumcounts = curr_counts.sum()
            maxcounts = curr_counts.max()
            tid_summary.at[tid, 'n_reads'] = sumcounts
            tid_summary.at[tid, 'peak_reads'] = maxcounts
            if sumcounts >= opts.minreads:
                if maxcounts < sumcounts * opts.peakfrac:
                    numseq = np.array(list(currtrans.get_sequence(genome).upper().translate(str_dict)))
                    curr_seq = ''.join(numseq)
                    tid_seq_info.append(pd.DataFrame({'tid': tid,
                                                      'genpos': curr_pos_list[psite:n_psite + psite],
                                                      'seq': np.array([(int(curr_seq[i:i + fpsize], 4) if 'N' not in curr_seq[i:i + fpsize] else -1)
                                                                       for i in xrange(n_psite)], dtype=np.int64),
                                                      'reads': curr_counts}))
                else:
                    tid_summary.at[tid, 'dropped'] = 'peakfrac'
            else:
                tid_summary.at[tid, 'dropped'] = 'lowreads'
    if tid_seq_info:  # don't bother saving anything if there's nothing to save
        pd.concat(tid_seq_info, ignore_index=True).to_hdf(seq_info_hdf % (chrom, strand), 'tid_seq_info', format='t',
                                                          data_columns=True, complevel=1, complib='blosc')
    #    sp.call(['ptrepack', orig_store_name, seq_info_hdf%(chrom,strand)])  # repack for efficiency
    #    os.remove(orig_store_name)
    if opts.verbose > 1:
        with log_lock:
            logprint('%s (%s strand) complete' % (chrom, strand))

    for inbam in inbams:
        inbam.close()

    return tid_summary


if opts.verbose:
    logprint('Parsing sequence and count information')

workers = mp.Pool(opts.numproc)
tid_summary = pd.concat(workers.map(_get_tid_info, bedlinedict.keys()))
workers.close()

if not (tid_summary['dropped'] == '').any():  # all transcripts dropped
    if not opts.keeptempfiles:
        try:
            os.rmdir(temp_folder)
        except OSError:
            pass  # try to remove the folder, but don't die if it's not empty
    lowreads_dropped = (tid_summary['dropped'] == 'lowreads').sum()
    raise ValueError('All %d transcripts dropped due to too few reads (%d) or too many reads coming from one position (%d). Consider decreasing '
                     'MINREADS (currently %d) or increasing PEAKFRAC (currently %f), or check validity of input BAM file.'
                     % (len(tid_summary), lowreads_dropped, len(tid_summary)-lowreads_dropped, opts.minreads, opts.peakfrac))

min_numseq = 0
max_numseq = 4 ** fpsize

npart = 64  # Divide sequences into this many separate dataframes based on sequence
# E.g. if npart == 4, then one dataframe will handle 'A'-initiated sequences, one 'C', etc
# This divides up the job naturally - because a sequence beginning with 'A' cannot multimap
# with one beginning with 'T' etc

partitions = np.ceil(np.linspace(min_numseq, max_numseq, npart + 1)).astype(np.int64)

# outname_orig = os.path.join(temp_folder,'tid_seq_mm_part_%d_orig.h5')
outname = os.path.join(temp_folder, 'tid_seq_mm_part_%d.h5')


def _find_mm_in_range(partnum):
    """Using the pandas HDF stores saved by _get_tid_info(), this function partitions reads based on their starting sequence, to divide up the
    problem of identifying multimapping positions. Each set of partitioned reads is saved to its own pandas HDF store."""
    seq_df = []
    for (chrom, strand) in bedlinedict.keys():
        fname = seq_info_hdf % (chrom, strand)
        if os.path.isfile(fname):
            seq_df.append(pd.read_hdf(fname, 'tid_seq_info', mode='r',
                                      where="seq >= %d & seq < %d" % (partitions[partnum], partitions[partnum + 1])))
            seq_df[-1]['chrom'] = chrom
            seq_df[-1]['strand'] = strand
    seq_df = pd.concat(seq_df, ignore_index=True)
    # Only care if the multimap is to a different genomic position
    uniqpos = seq_df.drop_duplicates(['chrom', 'strand', 'genpos', 'seq'])
    # Now, if a sequence appears multiple times, it means it's a multimapper (differing in chrom, strand, or position)
    seq_df = seq_df[seq_df['seq'].isin(uniqpos.loc[uniqpos['seq'].duplicated(), 'seq'].unique())]
    seq_df.to_hdf(outname % partnum, 'tid_seq_mm', format='t',
                  data_columns=True)  # save the result to avoid having to recalculate (after dropping pseudogenes)
    # seq_df.to_hdf(outname_orig%partnum,'tid_seq_mm',format='t',data_columns=True)
    # #  save the result to avoid having to recalculate (after dropping pseudogenes)
    # sp.call(['ptrepack', outname_orig%partnum, outname%partnum])  # repack for efficiency
    # os.remove(outname_orig%partnum)
    mm_df = seq_df.groupby('tid').agg({'genpos': len, 'reads': np.sum})
    if opts.verbose > 1:
        with log_lock:
            logprint('Partition %d of %d complete' % (partnum + 1, npart))
    return mm_df


if opts.verbose:
    logprint('Partitioning sequences to identify multimappers')

workers = mp.Pool(opts.numproc)
mm_df = workers.map(_find_mm_in_range, range(npart))
workers.close()
# Each dataframe in mm_df now contains the number of multimapping positions and reads, for sequences within that partition only
# Need to add them together to get the results for all sequences
mm_df_res = mm_df[0]
for df in mm_df[1:]:
    mm_df_res = mm_df_res.add(df, fill_value=0)  # this can handle empty dataframes (treated as all 0s), even if mm_df[0] is empty

if not opts.keeptempfiles:
    for (chrom, strand) in bedlinedict.keys():
        try:
            os.remove(seq_info_hdf % (chrom, strand))  # no longer needed
        except OSError:
            pass  # some files may not exist, in which case...no problem

tid_info = tid_summary[['chrom', 'strand', 'n_psite', 'n_reads']] \
    .join(mm_df_res.rename(columns={'genpos': 'mm_psite', 'reads': 'mm_reads'})) \
    .fillna({'mm_psite': 0,
             'mm_reads': 0})  # if a transcript isn't listed in mm_df_res, it doesn't have any multimapping positions
tid_info['reads_mm_frac'] = tid_info['mm_reads'] / tid_info['n_reads']
tid_info['psite_mm_frac'] = tid_info['mm_psite'] / tid_info['n_psite']
pseudos = tid_info.index[tid_info['reads_mm_frac'] > opts.pseudofrac].intersection(pd.Index(pseudotids))
if pseudos.size > 0:
    tid_summary.loc[pseudos, 'dropped'] = 'pseudo'
    kept_tids = tid_summary.index[tid_summary['dropped'] == '']


    def _find_kept_mm_in_range(partnum):
        """Load the multimapping positions on kept_tids from each partitioned dataframe"""
        try:
            seq_df = pd.read_hdf(outname % partnum, 'tid_seq_mm', mode='r')
        except IOError:
            return pd.DataFrame()
        seq_df = seq_df[seq_df['tid'].isin(kept_tids)]
        # Only care if the multimap is to a different genomic position
        uniqpos = seq_df.drop_duplicates(['chrom', 'strand', 'genpos', 'seq'])
        # Now, if a sequence appears multiple times, it means it's a multimapper (differing in chromosome, strand, and/or position)
        seq_df = seq_df[seq_df['seq'].isin(uniqpos.loc[uniqpos['seq'].duplicated(), 'seq'].unique())]
        return seq_df.groupby('tid').agg({'genpos': len, 'reads': np.sum})


    if opts.verbose:
        logprint('Recalculating multimappers after eliminating pseudogenes')
    workers = mp.Pool(opts.numproc)
    mm_df = workers.map(_find_kept_mm_in_range, range(npart))
    workers.close()
    mm_df_res = mm_df[0]
    for df in mm_df[1:]:
        mm_df_res = mm_df_res.add(df, fill_value=0)
    tid_info = tid_summary[['chrom', 'strand', 'n_psite', 'n_reads']] \
        .join(mm_df_res.rename(columns={'genpos': 'mm_psite', 'reads': 'mm_reads'})) \
        .fillna({'mm_psite': 0,
                 'mm_reads': 0})  # if a transcript isn't listed in mm_df_res, it doesn't have any multimapping positions, or it has been dropped
    tid_info['reads_mm_frac'] = tid_info['mm_reads'] / tid_info['n_reads']
    tid_info['psite_mm_frac'] = tid_info['mm_psite'] / tid_info['n_psite']

if not opts.keeptempfiles:
    for partnum in xrange(npart):
        try:
            os.remove(outname % partnum)  # no longer needed
        except OSError:
            pass
    try:
        os.rmdir(temp_folder)
    except OSError:
        pass  # try to remove the folder, but don't die if it's not empty

# Fraction of multimapped reads should be comparable to the fraction of multimapped positions
# Significantly greater multimapped read fraction than position fraction suggests that the reads "belong" in the other place, not here
# Really, should be doing a regression for the whole shebang - but that would be a computational nightmare
excess_mm = tid_info.index[(tid_info['reads_mm_frac'] > (tid_info['psite_mm_frac'] + opts.multiexcess))]
tid_summary.loc[excess_mm, 'dropped'] = 'multi'
if opts.summarytable:
    tid_summary.to_csv(opts.summarytable, sep='\t')
with open(opts.outbed, 'w') as outbed:
    for (tid, chrom, strand) in tid_summary.loc[tid_summary['dropped'] == '', ['chrom', 'strand']].itertuples(True):
        outbed.write(bedlinedict[(chrom, strand)][tid])

if opts.verbose:
    logprint('Tasks complete')
