import sys, os

import h5py

import numpy as np

from itertools import izip
from collections import defaultdict, namedtuple

from Bio import SeqIO

readData = namedtuple('readData', (
    'start', 'end', 'segs', 'read_start_rel_to_raw',
    'strand', 'means', 'stdevs', 'fn', 'corr_group'))
channelInfo = namedtuple(
    'channelInfo',
    ('offset', 'range', 'digitisation', 'number', 'sampling_rate'))
scaleValues = namedtuple(
    'scaleValues',
    ('shift', 'scale', 'lower_lim', 'upper_lim'))

NORM_TYPES = ('none', 'ont', 'median', 'robust_median')

# got quantiles from analysis of stability after shift-only normalization
robust_quantiles = (46.5, 53.5)

COMP_BASES = {'A':'T', 'C':'G', 'G':'C', 'T':'A', '-':'-', 'N':'N'}
def comp_base(base):
    # replace non-ACGT bases with dash
    try:
        return COMP_BASES[base]
    except KeyError:
        return 'N'
def rev_comp(seq):
    return ''.join(comp_base(b) for b in seq[::-1])

def parse_obs_filter(obs_filter):
    if obs_filter is None:
        return None

    # parse obs_filter
    try:
        obs_filter = [
            (float(pctl_nobs.split(':')[0]),
             int(pctl_nobs.split(':')[1]))
            for pctl_nobs in obs_filter]
    except:
        raise RuntimeError, 'Invalid format for observation filter'

    return obs_filter

def filter_reads(raw_read_coverage, obs_filter):
    if obs_filter is None:
        return raw_read_coverage

    num_reads = len([None for chrm_reads in raw_read_coverage.values()
                     for _ in chrm_reads])
    filt_raw_read_cov = {}
    for chrm, chrm_reads in raw_read_coverage.items():
        chrm_filt_reads = [
            r_data for r_data in chrm_reads if not any(
                np.percentile(np.diff(r_data.segs), pctl) > thresh
                for pctl, thresh in obs_filter)]
        if len(chrm_filt_reads) > 0:
            filt_raw_read_cov[chrm] = chrm_filt_reads
    num_filt_reads = len([
        None for chrm_reads in filt_raw_read_cov.values()
        for _ in chrm_reads])
    if num_filt_reads < num_reads:
        sys.stderr.write(
            'Filtered ' + str(num_reads - num_filt_reads) +
            ' reads due to observations per base filter from a ' +
            'total of ' + str(num_reads) + ' reads.\n')

    return filt_raw_read_cov


def parse_fast5s(files, corrected_group, basecall_subgroups,
                 get_means=False, get_stdevs=False):
    raw_read_coverage = defaultdict(list)
    for read_fn, basecall_subgroup in [
            (fn, bc_grp)for fn in files
            for bc_grp in basecall_subgroups]:
        try:
            read_data = h5py.File(read_fn, 'r')
        except IOError:
            # probably truncated file
            continue

        corr_slot = '/'.join((
            '/Analyses', corrected_group, basecall_subgroup))
        if corr_slot not in read_data:
            continue
        corr_data = read_data[corr_slot]

        try:
            align_data = dict(corr_data['Alignment'].attrs.items())
            read_start_rel_to_raw = corr_data['Events'].attrs[
                'read_start_rel_to_raw']
            event_data = corr_data['Events'].value
            events_end = event_data[-1]['start'] + event_data[-1][
                'length']
            segs = np.concatenate([event_data['start'], [events_end,]])
            base_means = event_data['norm_mean'] if get_means else None
            base_stdevs = event_data['norm_stdev'] if get_stdevs \
                          else None
        except:
            sys.stderr.write(
                '********** WARNING: Corrupted corrected events in ' +
                read_fn + '. ********\n')
            continue
        raw_read_coverage[align_data['mapped_chrom']].append(
            readData(
                align_data['mapped_start'],
                align_data['mapped_start'] + len(segs) - 1,
                segs, read_start_rel_to_raw,
                align_data['mapped_strand'], base_means, base_stdevs,
                read_fn, corrected_group + '/' + basecall_subgroup))

        read_data.close()

    return raw_read_coverage

def get_channel_info(fast5_data):
    try:
        fast5_info = fast5_data['UniqueGlobalKey/channel_id'].attrs
    except:
        raise RuntimeError, ("No channel_id group in HDF5 file. " +
                             "Probably mux scan HDF5 file.")

    channel_info = channelInfo(
        fast5_info['offset'], fast5_info['range'],
        fast5_info['digitisation'], fast5_info['channel_number'],
        fast5_info['sampling_rate'].astype('int_'))

    return channel_info

def normalize_raw_signal(
        all_raw_signal, read_start_rel_to_raw, read_obs_len,
        norm_type=None, channel_info=None, outlier_thresh=None,
        shift=None, scale=None, lower_lim=None, upper_lim=None):
    if norm_type not in NORM_TYPES and (shift is None or scale is None):
        raise NotImplementedError, (
            'Normalization type ' + norm_type + ' is not a valid ' +
            'option and shift or scale parameters were not provided.')

    raw_signal = np.array(all_raw_signal[
        read_start_rel_to_raw:
        read_start_rel_to_raw + read_obs_len])
    if shift is None or scale is None:
        if norm_type == 'none':
            shift, scale = 0, 1
        elif norm_type == 'ont':
            # correct raw signal as described here:
            #     https://community.nanoporetech.com
            #           /posts/squiggle-plot-for-raw-data
            shift, scale = (
                -1 * channel_info.offset,
                channel_info.digitisation / channel_info.range)
        elif norm_type == 'median':
            shift = np.median(raw_signal)
            scale = np.median(np.abs(raw_signal - shift))
        elif norm_type == 'robust_median':
            shift = np.mean(np.percentile(
                raw_signal, robust_quantiles))
            scale = np.median(np.abs(raw_signal - read_robust_med))

    raw_signal = (raw_signal - shift) / scale

    if outlier_thresh is not None or (
            lower_lim is not None and upper_lim is not None):
        if outlier_thresh is not None:
            read_med = np.median(raw_signal)
            read_mad = np.median(np.abs(raw_signal - read_med))
            lower_lim = read_med - (read_mad * outlier_thresh)
            upper_lim = read_med + (read_mad * outlier_thresh)
        raw_signal = np.array([
            upper_lim if outlier_high else
            (lower_lim if outlier_low else x)
            for x, outlier_high, outlier_low in izip(
                    raw_signal, raw_signal > upper_lim,
                    raw_signal < lower_lim)])

    return raw_signal, scaleValues(shift, scale, lower_lim, upper_lim)

def parse_fasta(fasta_fp):
    # could consider a conditional dependence on pyfaix if on-memory
    # indexing is required for larger genomes
    # Biopython's record level indexing will do for now...
    return SeqIO.index(fasta_fp,'fasta')

# Old fasta parser method
""" 
    curr_id = None
    curr_seq = None
    for line in fasta_fp:
        if line.startwith('>'):
            if (curr_id is not None and
                curr_seq is not None):
                fasta_records[curr_id] = curr_seq
                curr_id = line.replace(">","").split()[0]
                curr_seq = ''
        else:
            curr_seq += line.strip()

    if (curr_id is not None and
        curr_seq is not None):
        fasta_records[curr_id] = curr_seq

    return fasta_records
"""
