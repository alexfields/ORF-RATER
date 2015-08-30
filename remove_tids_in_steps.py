#! /usr/bin/env python

import argparse
import os
import sys
from yeti.genomics.genome_array import BAMGenomeArray, FivePrimeMapFactory, SizeFilterFactory
from yeti.genomics.roitools import SegmentChain  # , positionlist_to_segments
import pysam
from collections import defaultdict
from Bio import SeqIO
from string import maketrans
import pandas as pd
import numpy as np
import multiprocessing as mp
from time import strftime

parser = argparse.ArgumentParser()
parser.add_argument('--inbed', type=argparse.FileType('rU'), default=sys.stdin,
                    help='Transcriptome BED-file (Default: stdin)')
parser.add_argument('genomefasta', help='Path to genome FASTA-file')
parser.add_argument('bamfiles', nargs='+',
                    help='Path to transcriptome-aligned BAM file(s) to use for transcript filtering purposes')
parser.add_argument('--outfile',
                    help='Filename to use for tab-delimited text output (including column titles). First column is transcript IDs, '
                         'followed by summary information such as number of mappable positions and reads, maximum reads from any one position, '
                         'and why the transcript was dropped (if it was dropped)')
parser.add_argument('--outbed', type=argparse.FileType('w'), default=sys.stdout,
                    help='File to which to output BED-formatted transcripts that passed all filters (Default: stdout)')
parser.add_argument('--minlen', type=int, default=29,
                    help='Minimum length of read to be considered when evaluating transcripts. '
                         'Also serves as the size of the segment when identifying multimapping positions. (Default: 29)')
parser.add_argument('--maxlen', type=int, default=30,
                    help='Maximum length (inclusive) of read to be considered when evaluating transcripts. (Default: 30)')
parser.add_argument('-r', '--minreads', type=int, default=64,
                    help='Minimum number of reads demanded for each transcript (Default: 64)')
parser.add_argument('-f', '--peakfrac', type=float, default=1. / 5,
                    help='Maximum fraction of a transcript\'s reads coming from any one position (Default: 0.2)')
parser.add_argument('--pseudogenes',
                    help='List of transcript IDs annotated as pseudogenes, one per line. To be removed if >PSEUDOFRAC fraction of reads mapping to '
                         'each are multimappers.')
parser.add_argument('--pseudofrac', type=float, default=1. / 3,
                    help='Maximum allowable multimapping fraction of an annotated pseudogene\'s reads. '
                         'Ignored if list of pseudogenes is not provided. (Default: 0.333)')
parser.add_argument('--multiexcess', type=float, default=1. / 3,
                    help='Maximum disparity in multimapping reads versus multimapping positions for any transcript (Default: 0.333)')
parser.add_argument('-p', '--numproc', type=int, default=1,
                    help='Number of processes to run. Defaults to 1 but recommended to use more (e.g. 12-16)')
parser.add_argument('--keeptempfiles', action='store_true',
                    help='Keep the generated intermediate files (useful for debugging)')
parser.add_argument('--logfile', help='Generate a log file to record progress and timing')
opts = parser.parse_args()

if opts.logfile:
    with open(opts.logfile, 'w') as logfile:
        logfile.write(' '.join(sys.argv) + '\n')


    def logprint(nextstr):
        with open(opts.logfile, 'a') as logfile:
            logfile.write('[%s] %s\n' % (strftime('%Y-%m-%d %H:%M:%S'), nextstr))


    logprint('Reading transcriptome and genome')
    log_lock = mp.Lock()

fpsize = opts.minlen
psite = (fpsize + 1) / 2

if opts.pseudogenes:
    with open(opts.pseudogenes, 'rU') as inpseudo:
        pseudotids = {line.strip() for line in inpseudo}
else:
    pseudotids = {}

inbams = [pysam.Samfile(infile) for infile in opts.bamfiles]  # defaults to read mode, and will figure out if it's BAM or SAM
gnd = BAMGenomeArray(inbams, FivePrimeMapFactory(psite))
# map to roughly the center of each read so that identical sequences that cross different splice sites
# (on different transcripts) still end up mapping to the same place
gnd.add_filter('size', SizeFilterFactory(opts.minlen, opts.maxlen))

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

seq_info_hdf_orig = os.path.join(temp_folder, 'tid_seq_%s%s_orig.h5')
seq_info_hdf = os.path.join(temp_folder, 'tid_seq_%s%s.h5')


def _get_tid_info((chrom, strand)):
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
                                                      'seq': [(int(curr_seq[i:i + fpsize], 4) if 'N' not in curr_seq[i:i + fpsize] else -1)
                                                              for i in xrange(n_psite)],
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
    if opts.logfile:
        with log_lock:
            logprint('%s (%s strand) complete' % (chrom, strand))
    return tid_summary


if opts.logfile:
    logprint('Parsing sequence and count information')

workers = mp.Pool(opts.numproc)
tid_summary = pd.concat(workers.map(_get_tid_info, bedlinedict.keys()))
workers.close()

for inbam in inbams:
    inbam.close()

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
    seq_df = []
    for (chrom, strand) in bedlinedict.keys():
        fname = seq_info_hdf % (chrom, strand)
        if os.path.isfile(fname):
            seq_df.append(pd.read_hdf(fname, 'tid_seq_info',
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
    if opts.logfile:
        with log_lock:
            logprint('Partition %d of %d complete' % (partnum + 1, npart))
    return mm_df


if opts.logfile:
    logprint('Partitioning sequences to identify multimappers')

workers = mp.Pool(opts.numproc)
mm_df = workers.map(_find_mm_in_range, range(npart))
workers.close()
# Each dataframe in mm_df now contains the number of multimapping positions and reads, for sequences within that partition only
# Need to add them together to get the results for all sequences
mm_df_res = mm_df[0]
for df in mm_df[1:]:
    mm_df_res = mm_df_res.add(df, fill_value=0)

if not opts.keeptempfiles:
    for (chrom, strand) in bedlinedict.keys():
        try:
            os.remove(seq_info_hdf % (chrom, strand))  # no longer needed
        except:
            pass  # some files may not exist, in which case...no problem

tid_info = pd.concat((tid_summary[['chrom', 'strand', 'n_psite', 'n_reads']],
                      mm_df_res.rename(columns={'genpos': 'mm_psite', 'reads': 'mm_reads'})), axis=1) \
    .fillna({'mm_psite': 0,
             'mm_reads': 0})  # if a transcript isn't listed in mm_df_res, it doesn't have any multimapping positions
tid_info['reads_mm_frac'] = tid_info['mm_reads'] / tid_info['n_reads']
tid_info['psite_mm_frac'] = tid_info['mm_psite'] / tid_info['n_psite']
pseudos = tid_info.index[tid_info['reads_mm_frac'] > opts.pseudofrac].intersection(pd.Index(pseudotids))
if pseudos.size > 0:
    tid_summary.loc[pseudos, 'dropped'] = 'pseudo'
    kept_tids = tid_summary.index[tid_summary['dropped'] == '']


    def _find_kept_mm_in_range(partnum):
        seq_df = pd.read_hdf(outname % partnum, 'tid_seq_mm')
        seq_df = seq_df[seq_df['tid'].isin(kept_tids)]
        # Only care if the multimap is to a different genomic position
        uniqpos = seq_df.drop_duplicates(['chrom', 'strand', 'genpos', 'seq'])
        # Now, if a sequence appears multiple times, it means it's a multimapper (differing in chromosome, strand, and/or position)
        seq_df = seq_df[seq_df['seq'].isin(uniqpos.loc[uniqpos['seq'].duplicated(), 'seq'].unique())]
        return seq_df.groupby('tid').agg({'genpos': len, 'reads': np.sum})


    if opts.logfile:
        logprint('Recalculating multimappers after eliminating pseudogenes')
    workers = mp.Pool(opts.numproc)
    mm_df = workers.map(_find_kept_mm_in_range, range(npart))
    workers.close()
    mm_df_res = mm_df[0]
    for df in mm_df[1:]:
        mm_df_res = mm_df_res.add(df, fill_value=0)
    tid_info = pd.concat((tid_summary[['chrom', 'strand', 'n_psite', 'n_reads']],
                          mm_df_res.rename(columns={'genpos': 'mm_psite', 'reads': 'mm_reads'})), axis=1) \
        .fillna({'mm_psite': 0,
                 'mm_reads': 0})  # if a transcript isn't listed in mm_df_res, it doesn't have any multimapping positions, or it has been dropped
    tid_info['reads_mm_frac'] = tid_info['mm_reads'] / tid_info['n_reads']
    tid_info['psite_mm_frac'] = tid_info['mm_psite'] / tid_info['n_psite']

if not opts.keeptempfiles:
    for partnum in xrange(npart):
        try:
            os.remove(outname % partnum)  # no longer needed
        except:
            pass
    try:
        os.rmdir(temp_folder)
    except:
        pass  # try to remove the folder, but don't die if it's not empty

# Fraction of multimapped reads should be comparable to the fraction of multimapped positions
# Significantly greater multimapped read fraction than position fraction suggests that the reads "belong" in the other place, not here
# Really, should be doing a regression for the whole shebang - but that would be computationally challenging
excess_mm = tid_info.index[(tid_info['reads_mm_frac'] > (tid_info['psite_mm_frac'] + opts.multiexcess))]
tid_summary.loc[excess_mm, 'dropped'] = 'multi'
if opts.outfile:
    tid_summary.to_csv(opts.outfile, sep='\t')
# with open(opts.outbed,'rU') as outbed:
for (tid, chrom, strand) in tid_summary.loc[tid_summary['dropped'] == '', ['chrom', 'strand']].itertuples(True):
    opts.outbed.write(bedlinedict[(chrom, strand)][tid])

if opts.logfile:
    logprint('Tasks complete')
