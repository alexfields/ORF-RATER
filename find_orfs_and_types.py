#! /usr/bin/env python

import argparse
from Bio import SeqIO
from yeti.genomics.seqtools import seq_to_regex, IUPAC_TABLE_DNA
from yeti.genomics.roitools import Transcript, SegmentChain
import re
from collections import defaultdict
import pandas as pd
import numpy as np
import multiprocessing as mp
import subprocess as sp
import os
import sys
from time import strftime

parser = argparse.ArgumentParser(description='Identify all possible ORFs in a transcriptome. ORF-RATER will evaluate translation of only these ORFs.')
parser.add_argument('genomefasta', help='Path to genome FASTA-file')
parser.add_argument('--tfamstem', default='tfams', help='Transcript family information generated by make_tfams.py. Both TFAMSTEM.txt and '
                                                        'TFAMSTEM.bed should exist. (Default: tfams)')
parser.add_argument('--orfstore', default='orf.h5',
                    help='File to which to output the final table of identified ORFs. Will be formatted as a pandas HDF store (table name is '
                         '"all_orfs"). Different columns of the table indicate various of each ORF, such as start codon, length, etc. '
                         '(Default: orf.h5)')
parser.add_argument('--inbed', default='transcripts.bed', help='Transcriptome BED-file. Annotated CDSs are assumed to be bona fide CDSs, unless '
                                                               '--ignoreannotations is set. (Default: transcripts.bed)')
parser.add_argument('--codons', nargs='+', default=['ATG'],
                    help='Codons to consider as possible translation initiation sites. All must be 3 nucleotides long. Standard IUPAC nucleotide '
                         'codes are recognized; for example, to query all NTG codons, one could input "NTG" or "ATG CTG GTG TTG" (Default: ATG)')
parser.add_argument('--ignoreannotations', action='store_true', help='If flag is set, CDS annotations in INBED will be ignored. Typically used in '
                                                                     'conjunction with --extracdsbeds')
parser.add_argument('--extracdsbeds', nargs='+', help='Extra bed file(s) containing additional annotated CDSs beyond (or instead of) those in inbed. '
                                                      'Requires pybedtools.')
parser.add_argument('-v', '--verbose', action='store_true', help='Output a log of progress and timing (to stdout)')
parser.add_argument('-p', '--numproc', type=int, default=1, help='Number of processes to run. Defaults to 1 but more recommended if available.')
parser.add_argument('-f', '--force', action='store_true', help='Force file overwrite')
opts = parser.parse_args()

if not opts.force and os.path.exists(opts.orfstore):
    raise IOError('%s exists; use --force to overwrite' % opts.orfstore)

for codon in opts.codons:
    if len(codon) != 3 or any(x not in IUPAC_TABLE_DNA for x in codon.upper()):
        raise ValueError('%s is an invalid codon sequence' % codon)

if opts.verbose:
    sys.stdout.write(' '.join(sys.argv) + '\n')

    def logprint(nextstr):
        sys.stdout.write('[%s] %s\n' % (strftime('%Y-%m-%d %H:%M:%S'), nextstr))
        sys.stdout.flush()

    logprint('Reading transcriptome and genome')

START_RE = seq_to_regex('|'.join(opts.codons), nucleotide_table=IUPAC_TABLE_DNA)
STOP_RE = re.compile(r'(?:...)*?(?:TAG|TAA|TGA)')

# hash transcripts by ID for easy reference later
with open(opts.inbed, 'rU') as inbed:
    bedlinedict = {line.split()[3]: line for line in inbed}

tfamtids = defaultdict(list)
with open('%s.txt' % opts.tfamstem, 'rU') as tfamtable:
    for line in tfamtable:
        ls = line.strip().split()
        tfamtids[ls[1]].append(ls[0])

with open('%s.bed' % opts.tfamstem, 'rU') as tfambed:
    tfambedlines = {line.split()[3]: line for line in tfambed}

genome = SeqIO.to_dict(SeqIO.parse(opts.genomefasta, 'fasta'))

if not opts.ignoreannotations:
    annot_tfam_lookups = [tfamtids]
    annot_tid_lookups = [bedlinedict]
else:
    annot_tfam_lookups = []
    annot_tid_lookups = []

if opts.extracdsbeds:
    if opts.verbose:
        logprint('Identifying tfams for extra CDS annotations')
    import pybedtools  # to handle identifying which tfams get the extra CDSs - otherwise would need to replicate a lot of intersection functionality

    tfambedtool = pybedtools.BedTool('%s.bed' % opts.tfamstem)
    for cdsbedfile in opts.extracdsbeds:
        with open(cdsbedfile, 'rU') as cdsbed:
            annot_tid_lookups.append({line.split()[3]: line for line in cdsbed})  # as usual, hash bed lines by transcript ID
        annot_tfam_lookups.append(defaultdict(list))
        for line in tfambedtool.intersect(pybedtools.BedTool(cdsbedfile), split=True, s=True, wa=True, wb=True):
            annot_tfam_lookups[-1][line[3]].append(line[15])
# after this has finished, each element of annot_tfam_lookup will be a dictionary mapping tfams to lists of transcript IDs in the annotation bed files
# similarly, each element of annot_tid_lookup will map transcript IDs to BED lines

tfams_with_annots = set(sum([x.keys() for x in annot_tfam_lookups], []))


def _find_all_orfs(myseq):
    """Identify ORFs, or at least starts.
    Returns list of (start, stop, codon), where stop == 0 if no valid stop codon is present and codon is e.g. 'ATG'.
    Starts and stops are defined by START_RE and STOP_RE, respectively
    """
    result = []
    for i in range(len(myseq)-2):
        if START_RE.match(myseq[i:i+3]):
            m = STOP_RE.match(myseq[i:])
            if m:
                result.append((i, m.end()+i, myseq[i:i+3]))
            else:
                result.append((i, 0, myseq[i:i+3]))
    return result


def _name_orf(tfam, gcoord, AAlen):
    """Assign a usually unique identifier for each ORF. If not unique, a number should be appended to the end."""
    return '%s_%d_%daa' % (tfam, gcoord, AAlen)


def _identify_tfam_orfs((tfam, tids)):
    """Identify all of the possible ORFs within a family of transcripts. Relevant information such as genomic start and stop positions, amino acid
    length, and initiation codon will be collected for each ORF. Additionally, each ORF will be assigned a unique 'orfname', such that if it occurs
    on multiple transcripts, it can be recognized as the same ORF."""
    currtfam = SegmentChain.from_bed(tfambedlines[tfam])
    chrom = currtfam.chrom
    strand = currtfam.strand
    tfam_genpos = np.array(currtfam.get_position_list(stranded=True))
    tmask = np.empty((len(tids), len(tfam_genpos)), dtype=np.bool)  # True if transcript covers that position, False if not
    tfam_orfs = []
    tidx_lookup = {}
    for tidx, tid in enumerate(tids):
        tidx_lookup[tid] = tidx
        curr_trans = Transcript.from_bed(bedlinedict[tid])
        tmask[tidx, :] = np.in1d(tfam_genpos, curr_trans.get_position_list(stranded=True), assume_unique=True)
        trans_orfs = _find_all_orfs(curr_trans.get_sequence(genome).upper())
        if trans_orfs:
            (startpos, stoppos, codons) = zip(*trans_orfs)
            startpos = np.array(startpos)
            stoppos = np.array(stoppos)

            gcoords = np.array(curr_trans.get_genomic_coordinate(startpos)[1], dtype='u4')

            stop_present = (stoppos > 0)
            gstops = np.zeros(len(trans_orfs), dtype='u4')
            gstops[stop_present] = curr_trans.get_genomic_coordinate(stoppos[stop_present] - 1)[1] + (strand == '+')*2 - 1
            # the decrementing/incrementing stuff preserves half-openness regardless of strand

            AAlens = np.zeros(len(trans_orfs), dtype='u4')
            AAlens[stop_present] = (stoppos[stop_present] - startpos[stop_present])/3 - 1
            tfam_orfs.append(pd.DataFrame.from_items([('tfam', tfam),
                                                      ('tid', tid),
                                                      ('tcoord', startpos),
                                                      ('tstop', stoppos),
                                                      ('chrom', chrom),
                                                      ('gcoord', gcoords),
                                                      ('gstop', gstops),
                                                      ('strand', strand),
                                                      ('codon', codons),
                                                      ('AAlen', AAlens),
                                                      ('orfname', '')]))
    if any(x is not None for x in tfam_orfs):
        orf_pos_dict = {}
        tfam_orfs = pd.concat(tfam_orfs, ignore_index=True)
        for ((gcoord, AAlen), gcoord_grp) in tfam_orfs.groupby(['gcoord', 'AAlen']):  # group by genomic start position and length
            if len(gcoord_grp) == 1:
                tfam_orfs.loc[gcoord_grp.index, 'orfname'] = _name_orf(tfam, gcoord, AAlen)
            else:
                orf_gcoords = np.vstack(np.flatnonzero(tmask[tidx_lookup[tid], :])[tcoord:tstop]
                                        for (tid, tcoord, tstop) in gcoord_grp[['tid', 'tcoord', 'tstop']].itertuples(False))
                if (orf_gcoords == orf_gcoords[0, :]).all():  # all of the grouped ORFs are identical, so should receive the same name
                    orfname = _name_orf(tfam, gcoord, AAlen)
                    tfam_orfs.loc[gcoord_grp.index, 'orfname'] = orfname
                    orf_pos_dict[orfname] = tfam_genpos[orf_gcoords[0, :]]
                else:
                    named_so_far = 0
                    unnamed = np.ones(len(gcoord_grp), dtype=np.bool)
                    basename = _name_orf(tfam, gcoord, AAlen)
                    while unnamed.any():
                        next_gcoords = orf_gcoords[unnamed, :][0, :]
                        identicals = (orf_gcoords == next_gcoords).all(1)
                        orfname = '%s_%d' % (basename, named_so_far)
                        tfam_orfs.loc[gcoord_grp.index[identicals], 'orfname'] = orfname
                        orf_pos_dict[orfname] = tfam_genpos[next_gcoords]
                        unnamed[identicals] = False
                        named_so_far += 1

        # Now that the ORFs have been found and named, figure out their orftype
        tfam_orfs['annot_start'] = False
        tfam_orfs['annot_stop'] = False  # start out assuming all are False; replace with True as needed
        tfam_orfs['orftype'] = 'new'
        tfam_orfs['untyped'] = tfam_orfs['tstop'] > 0
        tfam_orfs.loc[~tfam_orfs['untyped'], 'orftype'] = 'nonstop'  # no stop codon
        if tfam in tfams_with_annots:
            cds_info = []
            all_annot_pos = set()
            for (annot_fidx, (annot_tfam_lookup, annot_tid_lookup)) in enumerate(zip(annot_tfam_lookups, annot_tid_lookups)):
                if tfam in annot_tfam_lookup:
                    for (annot_tidx, annot_tid) in enumerate(annot_tfam_lookup[tfam]):
                        curr_trans = Transcript.from_bed(annot_tid_lookup[annot_tid])
                        if curr_trans.cds_start is not None and curr_trans.cds_end is not None:
                            curr_cds_pos_set = curr_trans.get_cds().get_position_set()
                            curr_len = len(curr_cds_pos_set)
                            if curr_len % 3 == 0:
                                curr_gcoord = curr_trans.get_genomic_coordinate(curr_trans.cds_start)[1]
                                curr_gstop = curr_trans.get_genomic_coordinate(curr_trans.cds_end - 1)[1] + (strand == '+') * 2 - 1
                                in_tfam = curr_cds_pos_set.issubset(tfam_genpos)
                                cds_info.append((curr_gcoord, curr_gstop, (curr_len-3)/3, in_tfam, annot_fidx, annot_tid, curr_cds_pos_set))
                                all_annot_pos.update(curr_cds_pos_set)
            if cds_info:  # False means no annotated CDSs or none are multiples of 3 in length
                cds_info = pd.DataFrame(cds_info, columns=['gcoord', 'gstop', 'AAlen', 'in_tfam', 'annot_fidx', 'annot_tid', 'pos']) \
                    .groupby(['gcoord', 'gstop', 'AAlen', 'in_tfam'], as_index=False) \
                    .apply(lambda x: x if len(x) == 1 else x[[not any(pos == x['pos'].iat[j] for j in xrange(i))
                                                              for (i, pos) in enumerate(x['pos'])]]) \
                    .set_index(['annot_fidx', 'annot_tid'])
                # this operation organizes cds_info into a dataframe and effectively drops duplicates
                # pandas drop_duplicates() is incompatible with sets so have to do it this manual way
                # the combination of annot_fidx (the number of the file if more than one annotation file provided) and annot_tid should be a unique ID
                tfam_orfs['annot_start'] = tfam_orfs['gcoord'].isin(cds_info['gcoord'])
                tfam_orfs['annot_stop'] = tfam_orfs['gstop'].isin(cds_info['gstop'])

                def _get_orf_pos(orfname, tid=None, tcoord=None, tstop=None):
                    """Helper function that identifies the genomic coordinates of an ORF (in stranded order) and caches them by orfname"""
                    if orfname in orf_pos_dict:
                        return orf_pos_dict[orfname]
                    else:
                        if tid is None or tcoord is None or tstop is None:
                            (tid, tcoord, tstop) = tfam_orfs.loc[tfam_orfs['orfname'] == orfname, ['tid', 'tcoord', 'tstop']].iloc[0]
                        res = tfam_genpos[np.flatnonzero(tmask[tidx_lookup[tid], :])[tcoord:tstop]]
                        orf_pos_dict[orfname] = res
                        return res

                # ANNOTATED and XISO
                cds_info['found'] = False
                possible_annot = tfam_orfs.drop_duplicates('orfname').merge(cds_info[cds_info['in_tfam']].reset_index())
                # merges on gcoord, gstop, and len - need to reset_index to preserve annot_fidx and annot_tid
                for ((orfname, tid, tcoord, tstop), cds_grp) in possible_annot.groupby(['orfname', 'tid', 'tcoord', 'tstop']):
                    orf_pos = _get_orf_pos(orfname, tid, tcoord, tstop)
                    for (annot_fidx, annot_tid, cds_pos_set) in cds_grp[['annot_fidx', 'annot_tid', 'pos']].itertuples(False):
                        if cds_pos_set.issubset(orf_pos):
                            tfam_orfs.loc[tfam_orfs['orfname'] == orfname, ['orftype', 'untyped']] = ['annotated', False]
                            cds_info.loc[(annot_fidx, annot_tid), 'found'] = True
                            break
                    else:
                        tfam_orfs.loc[tfam_orfs['orfname'] == orfname, ['orftype', 'untyped']] = ['Xiso', False]
                        # matching start and stop but differing in between
                if tfam_orfs['untyped'].any():
                    tfam_orfs.loc[tfam_orfs['orfname'].isin(tfam_orfs[tfam_orfs['untyped']].merge(cds_info[['gcoord', 'gstop']])['orfname']),
                                  ['orftype', 'untyped']] = ['Xiso', False]
                    # matching start and stop, but must differ somewhere, otherwise would have been identified as annotated (Xiso => "exact isoform")

                # SISO
                tfam_orfs.loc[tfam_orfs['annot_start'] & tfam_orfs['annot_stop'] & tfam_orfs['untyped'], ['orftype', 'untyped']] = ['Siso', False]
                # start and stop each match at least one CDS, but not the same one (Siso => "spliced isoform")

                # CISO
                tfam_orfs.loc[tfam_orfs['annot_start'] & tfam_orfs['untyped'], ['orftype', 'untyped']] = ['Ciso', False]
                # start is annotated, but stop is not - so must be on a new transcript (Ciso => "C-terminal isoform")

                # TRUNCATION
                if tfam_orfs['untyped'].any():
                    found_matched_stop = tfam_orfs[tfam_orfs['untyped']].merge(tfam_orfs[tfam_orfs['orftype'] == 'annotated'],
                                                                               on=['tid', 'tstop'], suffixes=('', '_annot'))
                    tfam_orfs.loc[tfam_orfs['orfname'].isin(found_matched_stop.loc[found_matched_stop['tcoord'] > found_matched_stop['tcoord_annot'],
                                                                                   'orfname']), ['orftype', 'untyped']] = ['truncation', False]
                # on the same transcript with an annotated CDS, with matching stop codon, initiating downstream - must be a truncation
                # still some missing truncations, if the original CDS was not on a transcript in the present transcriptome
                if tfam_orfs['untyped'].any() and not cds_info['found'].all():
                    possible_truncs = tfam_orfs[tfam_orfs['untyped']].drop_duplicates('orfname') \
                        .merge(cds_info.loc[~cds_info['found'], ['gstop', 'pos', 'AAlen']], on='gstop', suffixes=('', '_annot'))
                    possible_truncs = possible_truncs[possible_truncs['AAlen'] < possible_truncs['AAlen_annot']]
                    for ((orfname, tid, tcoord, tstop, gcoord), cds_pos_sets) in \
                            possible_truncs.groupby(['orfname', 'tid', 'tcoord', 'tstop', 'gcoord'])['pos']:
                        orf_pos = _get_orf_pos(orfname, tid, tcoord, tstop)
                        if strand == '-':
                            if any(cds_pos_set.issuperset(orf_pos) and
                                   all(pos in orf_pos for pos in cds_pos_set if pos <= gcoord) for cds_pos_set in cds_pos_sets):
                                tfam_orfs.loc[tfam_orfs['orfname'] == orfname, ['orftype', 'untyped']] = ['truncation', False]
                        else:
                            if any(cds_pos_set.issuperset(orf_pos) and
                                   all(pos in orf_pos for pos in cds_pos_set if pos >= gcoord) for cds_pos_set in cds_pos_sets):
                                tfam_orfs.loc[tfam_orfs['orfname'] == orfname, ['orftype', 'untyped']] = ['truncation', False]
                        # matching stop codon, contained within, and all positions in the annotation past the orf start codon are included in the orf

                # EXTENSION
                if tfam_orfs['untyped'].any():
                    found_matched_stop = tfam_orfs[tfam_orfs['untyped']].merge(tfam_orfs[tfam_orfs['orftype'] == 'annotated'],
                                                                               on=['tid', 'tstop'], suffixes=('', '_annot'))
                    assert (found_matched_stop['tcoord'] < found_matched_stop['tcoord_annot']).all()  # other possibilities should be done by now
                    tfam_orfs.loc[tfam_orfs['orfname'].isin(found_matched_stop['orfname']), ['orftype', 'untyped']] = ['extension', False]
                # on the same transcript with an annotated CDS, with matching stop codon, initiating upstream - must be an extension
                # no possibility for an "unfound" extension - if the extension is in the transcriptome, the CDS it comes from must be as well
                # (except for a few edge cases e.g. annotated CDS is a CUG initiator, but not considering CUG ORFs)

                # NISO
                tfam_orfs.loc[tfam_orfs['annot_stop'] & (tfam_orfs['untyped']), ['orftype', 'untyped']] = ['Niso', False]
                # stop is annotated, but start is not, and it's not a truncation or extension - so must be an isoform (Niso => "N-terminal isoform")

                # NCISO
                if tfam_orfs['untyped'].any():
                    orf_codons = []
                    for (orfname, tid, tcoord, tstop) in \
                            tfam_orfs.loc[tfam_orfs['untyped'], ['orfname', 'tid', 'tcoord', 'tstop']].drop_duplicates('orfname').itertuples(False):
                        orf_codons.append(pd.DataFrame(_get_orf_pos(orfname, tid, tcoord, tstop).reshape((-1, 3))))
                        orf_codons[-1]['orfname'] = orfname
                    orf_codons = pd.concat(orf_codons, ignore_index=True)
                    if strand == '-':
                        annot_codons = pd.DataFrame(np.vstack([np.reshape(sorted(cds_pos_set, reverse=True), (-1, 3))
                                                               for cds_pos_set in cds_info['pos'] if len(cds_pos_set) % 3 == 0])).drop_duplicates()
                    else:
                        annot_codons = pd.DataFrame(np.vstack([np.reshape(sorted(cds_pos_set, reverse=False), (-1, 3))
                                                               for cds_pos_set in cds_info['pos'] if len(cds_pos_set) % 3 == 0])).drop_duplicates()
                    tfam_orfs.loc[tfam_orfs['orfname'].isin(orf_codons.merge(annot_codons)['orfname']), ['orftype', 'untyped']] = ['NCiso', False]
                    # ORFs that have at least one full codon overlapping (in-frame) with a CDS are isoforms (NCiso => "N- and C-terminal isoform")
                    # Note that these must already differ at N- and C- termini, otherwise they would already have been classified

                # INTERNAL
                if tfam_orfs['untyped'].any():
                    sametrans = tfam_orfs[tfam_orfs['untyped']].merge(tfam_orfs[tfam_orfs['orftype'] == 'annotated'],
                                                                      on='tid', suffixes=('', '_annot'))
                    sametrans_internal = (sametrans['tcoord'] > sametrans['tcoord_annot']) & (sametrans['tstop'] < sametrans['tstop_annot'])
                    tfam_orfs.loc[tfam_orfs['orfname'].isin(sametrans.loc[sametrans_internal, 'orfname']),
                                  ['orftype', 'untyped']] = ['internal', False]
                # ORFs completely contained within a CDS on the same transcript, and not containing any full codon overlaps, must be internal
                # Still could be other ORFs internal to a CDS on a transcript not in the current transcriptome - need to check manually

                if tfam_orfs['untyped'].any() and not cds_info['found'].all():
                    for (orfname, gcoord, gstop) in \
                            tfam_orfs.loc[tfam_orfs['untyped'], ['orfname', 'gcoord', 'gstop']].drop_duplicates('orfname').itertuples(False):
                        orf_pos = _get_orf_pos(orfname)  # should be cached by now
                        if strand == '-':
                            if any(cds_pos_set.issuperset(orf_pos) and all(pos in orf_pos for pos in cds_pos_set if gcoord >= pos > gstop)
                                   for cds_pos_set in cds_info.loc[(~cds_info['found'])
                                                                   & (cds_info['gcoord'] > gcoord)
                                                                   & (cds_info['gstop'] < gstop), 'pos']):
                                tfam_orfs.loc[tfam_orfs['orfname'] == orfname, ['orftype', 'untyped']] = ['internal', False]
                        else:
                            if any(cds_pos_set.issuperset(orf_pos) and all(pos in orf_pos for pos in cds_pos_set if gcoord <= pos < gstop)
                                   for cds_pos_set in cds_info.loc[(~cds_info['found'])
                                                                   & (cds_info['gcoord'] < gcoord)
                                                                   & (cds_info['gstop'] > gstop), 'pos']):
                                tfam_orfs.loc[tfam_orfs['orfname'] == orfname, ['orftype', 'untyped']] = ['internal', False]

                # STOP_OVERLAP
                if tfam_orfs['untyped'].any():
                    sametrans = tfam_orfs[tfam_orfs['untyped']].merge(tfam_orfs[tfam_orfs['orftype'] == 'annotated'],
                                                                      on='tid', suffixes=('', '_annot'))
                    sametrans_stopover = (sametrans['tcoord'] > sametrans['tcoord_annot']) & (sametrans['tcoord'] < sametrans['tstop_annot'])
                    tfam_orfs.loc[tfam_orfs['orfname'].isin(sametrans.loc[sametrans_stopover, 'orfname']),
                                  ['orftype', 'untyped']] = ['stop_overlap', False]
                    # starts within a CDS and not an internal - must be a stop_overlap
                    # do not need to check for unfounds - requiring that stop_overlap must be on same transcript as cds

                # START_OVERLAP
                if tfam_orfs['untyped'].any():
                    sametrans = tfam_orfs[tfam_orfs['untyped']].merge(tfam_orfs[tfam_orfs['orftype'] == 'annotated'],
                                                                      on='tid', suffixes=('', '_annot'))
                    sametrans_startover = (sametrans['tstop'] > sametrans['tcoord_annot']) & (sametrans['tstop'] < sametrans['tstop_annot'])
                    tfam_orfs.loc[tfam_orfs['orfname'].isin(sametrans.loc[sametrans_startover, 'orfname']),
                                  ['orftype', 'untyped']] = ['start_overlap', False]
                    # ends within a CDS and not an internal - must be a start_overlap
                    # do not need to check for unfounds - requiring that start_overlap must be on same transcript as cds

                # LOOF
                if tfam_orfs['untyped'].any():
                    sametrans = tfam_orfs[tfam_orfs['untyped']].merge(tfam_orfs[tfam_orfs['orftype'] == 'annotated'],
                                                                      on='tid', suffixes=('', '_annot'))
                    sametrans_loof = (sametrans['tcoord'] < sametrans['tcoord_annot']) & (sametrans['tstop'] > sametrans['tstop_annot'])
                    tfam_orfs.loc[tfam_orfs['orfname'].isin(sametrans.loc[sametrans_loof, 'orfname']), ['orftype', 'untyped']] = ['LOOF', False]
                    # starts upstream of a CDS and ends downstream of it - must be a LOOF (long out-of-frame)
                    # don't need to check for unfounds because the CDS must be on the same transcript as the ORF if the ORF completely contains it

                # UPSTREAM
                if tfam_orfs['untyped'].any():
                    sametrans = tfam_orfs[tfam_orfs['untyped']].merge(tfam_orfs[tfam_orfs['orftype'] == 'annotated'],
                                                                      on='tid', suffixes=('', '_annot'))
                    sametrans_upstream = (sametrans['tstop'] <= sametrans['tcoord_annot'])
                    tfam_orfs.loc[tfam_orfs['orfname'].isin(sametrans.loc[sametrans_upstream, 'orfname']),
                                  ['orftype', 'untyped']] = ['upstream', False]
                    # ends upstream of a CDS - must be an upstream (uORF)
                    # cannot check manually for unfounds because those are not on well-defined transcripts

                # DOWNSTREAM
                if tfam_orfs['untyped'].any():
                    sametrans = tfam_orfs[tfam_orfs['untyped']].merge(tfam_orfs[tfam_orfs['orftype'] == 'annotated'],
                                                                      on='tid', suffixes=('', '_annot'))
                    sametrans_downstream = (sametrans['tstop_annot'] <= sametrans['tcoord'])
                    tfam_orfs.loc[tfam_orfs['orfname'].isin(sametrans.loc[sametrans_downstream, 'orfname']),
                                  ['orftype', 'untyped']] = ['downstream', False]
                    # starts downstream of a CDS - must be an upstream (uORF)
                    # cannot check manually for unfounds because those are not on well-defined transcripts

                # NEW_ISO and GISO
                for orfname in tfam_orfs.loc[tfam_orfs['untyped'], 'orfname'].drop_duplicates():
                    if all_annot_pos.isdisjoint(_get_orf_pos(orfname)):
                        tfam_orfs.loc[tfam_orfs['orfname'] == orfname, ['orftype', 'untyped']] = ['new_iso', False]
                        # no overlaps whatsoever with any annotated CDS, but in a tfam that has annotations: new_iso
                    else:
                        tfam_orfs.loc[tfam_orfs['orfname'] == orfname, ['orftype', 'untyped']] = ['Giso', False]
                        # overlaps out-of-frame with a CDS, and not on the same transcript with a CDS: Giso => "genomic isoform"

                assert not tfam_orfs['untyped'].any()
        return tfam_orfs.drop('untyped', axis=1)
    else:
        return None

if opts.verbose:
    logprint('Identifying ORFs within each transcript family')

workers = mp.Pool(opts.numproc)
all_orfs = pd.concat(workers.map(_identify_tfam_orfs, tfamtids.iteritems()), ignore_index=True)
workers.close()

for catfield in ['chrom', 'strand', 'codon', 'orftype']:
    all_orfs[catfield] = all_orfs[catfield].astype('category')  # saves disk space and read/write time

if opts.verbose:
    logprint('Saving results')

origname = opts.orfstore+'.tmp'
all_orfs.to_hdf(origname, 'all_orfs', format='t', data_columns=True)
sp.call(['ptrepack', origname, opts.orfstore])  # repack for efficiency
os.remove(origname)

if opts.verbose:
    logprint('Tasks complete')
