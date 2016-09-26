import h5py

import numpy as np

from itertools import izip
from collections import defaultdict, namedtuple

readData = namedtuple('readData', (
    'start', 'end', 'segs', 'read_start_rel_to_raw',
    'strand', 'means', 'fn', 'corr_group'))

NORM_TYPES = ('none', 'ont', 'median', 'robust_median')

# got quantiles from analysis of stability after shift-only normalization
robust_quantiles = (46.5, 53.5)


def parse_fast5s(files, corrected_group, get_means=False):
    raw_read_coverage = defaultdict(list)
    for read_fn in files:
        try:
            read_data = h5py.File(read_fn, 'r')
        except IOError:
            # probably truncated file
            continue

        if 'Analyses/' + corrected_group not in read_data:
            continue
        align_data = read_data['Analyses/' + corrected_group +
                               '/Alignment/'].attrs
        seg_grp = read_data[
            'Analyses/' + corrected_group + '/template/Segments']
        read_start_rel_to_raw = seg_grp.attrs[
            'read_start_rel_to_raw']
        segs = seg_grp.value
        base_means = read_data[
            'Analyses/' + corrected_group +
            '/template/Events']['norm_mean'] if get_means else None
        raw_read_coverage[align_data['mapped_chrom']].append(
            readData(
                align_data['mapped_start'],
                align_data['mapped_start'] + len(segs) - 1,
                segs, read_start_rel_to_raw,
                align_data['mapped_strand'], base_means, read_fn,
                corrected_group))

        read_data.close()

    return raw_read_coverage

def normalize_raw_signal(
        all_raw_signal, read_start_rel_to_raw, read_obs_len,
        norm_type, channel_info, outlier_threshold):
    if norm_type not in NORM_TYPES:
        raise NotImplementedError, (
            'Normalization type ' + norm_type +
            ' is not a valid option. Please select a ' +
            'valid normalization method.')

    raw_signal = np.array(all_raw_signal[
        read_start_rel_to_raw:
        read_start_rel_to_raw + read_obs_len])
    shift, scale = 0, 1
    if norm_type == 'ont':
        # correct raw signal as described here:
        #     https://community.nanoporetech.com
        #           /posts/squiggle-plot-for-raw-data
        raw_signal = (((raw_signal + channel_info.offset) *
                       channel_info.range) / channel_info.digitisation)
        shift, scale = (
            channel_info.offset * channel_info.range /
            channel_info.digitisation,
            channel_info.range / channel_info.digitisation)
    elif norm_type == 'median':
        read_med = np.median(raw_signal)
        read_mad = np.median(np.abs(raw_signal - read_med))
        raw_signal = (raw_signal - read_med) / read_mad
        shift, scale = read_med, read_mad
    elif norm_type == 'robust_median':
        read_robust_med = np.mean(np.percentile(
            raw_signal, robust_quantiles))
        read_mad = np.median(np.abs(raw_signal - read_robust_med))
        raw_signal = (raw_signal - read_robust_med) / read_mad
        shift, scale = read_robust_med, read_mad

    if outlier_threshold is not None:
        read_med = np.median(raw_signal)
        read_mad = np.median(np.abs(raw_signal - read_med))
        lower_lim = read_med - (read_mad * outlier_threshold)
        upper_lim = read_med + (read_mad * outlier_threshold)
        raw_signal = np.array([
            upper_lim if outlier_high else
            (lower_lim if outlier_low else x)
            for x, outlier_high, outlier_low in izip(
                    raw_signal, raw_signal > upper_lim,
                    raw_signal < lower_lim)])

    return raw_signal, shift, scale

def parse_fasta(fasta_fp):
    fasta_records = {}
    try:
        from Bio import SeqIO
        seqio_fasta_recs = SeqIO.parse(fasta_fp,'fasta')
        for fasta in seqio_fasta_recs:
            fasta_records[fasta.id] = fasta.seq.tostring()

        return fasta_records
    except ImportError:
        sys.stderr.write(
            'WARNING: Could not load "Bio" python module so ' +
            'using less robust parser.\n')
        pass

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
