import sys, os

import re
import h5py

import numpy as np

from itertools import izip
from collections import defaultdict, namedtuple

NANORAW_VERSION = '0.4.2'

SMALLEST_PVAL = 1e-20

readData = namedtuple('readData', (
    'start', 'end', 'segs', 'read_start_rel_to_raw',
    'strand', 'fn', 'corr_group'))
channelInfo = namedtuple(
    'channelInfo',
    ('offset', 'range', 'digitisation', 'number', 'sampling_rate'))
scaleValues = namedtuple(
    'scaleValues',
    ('shift', 'scale', 'lower_lim', 'upper_lim'))

NORM_TYPES = ('none', 'pA', 'pA_raw', 'median', 'robust_median')

# single base conversion for motifs
SINGLE_LETTER_CODE = {
    'A':'A', 'C':'C', 'G':'G', 'T':'T', 'B':'[CGT]',
    'D':'[AGT]', 'H':'[ACT]', 'K':'[GT]', 'M':'[AC]',
    'N':'[ACGT]', 'R':'[AG]', 'S':'[CG]', 'V':'[ACG]',
    'W':'[AT]', 'Y':'[CT]'}

VERBOSE = False

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

def get_files_list(basedirs):
    return [os.path.join(base_dir, fn)
            for base_dir in basedirs
            for fn in os.listdir(base_dir)]
def get_files_lists(basedirs1, basedirs2):
    files1 = get_files_list(basedirs1)
    files2 = get_files_list(basedirs2) if basedirs2 else None

    return files1, files2

def parse_motif(motif):
    invalid_chars = re.findall('[^ACGTBDHKMNRSVWY]', motif)
    if len(invalid_chars) > 0:
       sys.stderr.write(
           '********* ERROR ********* Invalid characters in motif: ' +
           ', '.join(invalid_chars) + ' *********\n')
       sys.exit()

    return re.compile(''.join(
        SINGLE_LETTER_CODE[letter] for letter in motif))

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

def parse_fast5s(files, corrected_group, basecall_subgroups):
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
        except:
            sys.stderr.write(
                '********** WARNING: Corrupted corrected events in ' +
                read_fn + '. ********\n')
            continue
        raw_read_coverage[(
            align_data['mapped_chrom'],
            align_data['mapped_strand'])].append(
                readData(
                    align_data['mapped_start'],
                    align_data['mapped_start'] + len(segs) - 1,
                    segs, read_start_rel_to_raw,
                    align_data['mapped_strand'], read_fn,
                    corrected_group + '/' + basecall_subgroup))

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

def parse_pore_model(pore_model_fn):
    pore_model = {'mean':{}, 'inv_var':{}}
    with open(pore_model_fn) as fp:
        for line in fp:
            if line.startswith('#'): continue
            try:
                kmer, lev_mean, lev_stdev = line.split()[:3]
                lev_mean, lev_stdev = map(float, (lev_mean, lev_stdev))
            except ValueError:
                # header or other non-kmer field
                continue
            pore_model['mean'][kmer] = lev_mean
            pore_model['inv_var'][kmer] = 1 / (lev_stdev * lev_stdev)

    return pore_model

def calc_kmer_fitted_shift_scale(pore_model, events_means, events_kmers):
    r_model_means = np.array([pore_model['mean'][kmer]
                              for kmer in events_kmers])
    r_model_inv_vars = np.array([pore_model['inv_var'][kmer]
                                 for kmer in events_kmers])
    model_mean_var = r_model_means * r_model_inv_vars
    # prep kmer model coefficient matrix for the k-mers from this read
    model_mean_var_sum = model_mean_var.sum()
    coef_mat = np.array((
        (r_model_inv_vars.sum(), model_mean_var_sum),
        (model_mean_var_sum, (model_mean_var * r_model_means).sum())))

    # prep dependent values from this reads true events
    r_event_var = events_means * r_model_inv_vars
    r_event_var_mean = r_event_var * r_model_means
    dep_vect = np.array((r_event_var.sum(), r_event_var_mean.sum()))

    shift, scale = np.linalg.solve(coef_mat, dep_vect)

    return shift, scale

def normalize_raw_signal(
        all_raw_signal, read_start_rel_to_raw, read_obs_len,
        norm_type=None, channel_info=None, outlier_thresh=None,
        shift=None, scale=None, lower_lim=None, upper_lim=None,
        pore_model=None, event_means=None, event_kmers=None):
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
        elif norm_type in ('pA_raw', 'pA'):
            # correct raw signal as described here:
            #     https://community.nanoporetech.com
            #           /posts/squiggle-plot-for-raw-data
            shift, scale = (
                -1 * channel_info.offset,
                channel_info.digitisation / channel_info.range)
            if norm_type == 'pA':
                # perform k-mer model fitted correction as in
                # nanocorr/nanopolish/ONT
                fit_shift, fit_scale = calc_kmer_fitted_shift_scale(
                    pore_model, event_means, event_kmers)
                # apply shift and scale values fitted from kmer
                # conditional model after raw DAC scaling
                shift = shift + (fit_shift * scale)
                scale = scale * fit_scale
            # print fitted shift and scale for comparisons
            #print 'shift: ' + str(fit_shift) + \
            #  '\tscale: ' + str(fit_scale)
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

def parse_fasta(fasta_fn):
    # Tried Biopython index and that opened the fail again for each
    # record access request and was thus far too slow

    # could consider a conditional dependence on pyfaix if on-memory
    # indexing is required for larger genomes
    # testing shows that human genome only takes 3.2 GB with raw parser
    # so raw parsing is probably fine
    fasta_fp = open(fasta_fn)

    fasta_records = {}
    curr_id = None
    curr_seq = ''
    for line in fasta_fp:
        if line.startswith('>'):
            if (curr_id is not None and
                curr_seq is not ''):
                fasta_records[curr_id] = curr_seq
            curr_seq = ''
            curr_id = line.replace(">","").strip().split()[0]
        else:
            curr_seq += line.strip()

    # add last record
    if (curr_id is not None and
        curr_seq is not ''):
        fasta_records[curr_id] = curr_seq

    fasta_fp.close()

    return fasta_records

def get_coverage(raw_read_coverage):
    if VERBOSE: sys.stderr.write('Calculating read coverage.\n')
    read_coverage = {}
    for (chrm, strand), reads_data in raw_read_coverage.items():
        max_end = max(r_data.end for r_data in reads_data)
        chrm_coverage = np.zeros(max_end, dtype=np.int_)
        for r_data in reads_data:
            chrm_coverage[r_data.start:r_data.end] += 1
        read_coverage[(chrm, strand)] = chrm_coverage

    return read_coverage

def get_read_base_means(r_data):
    try:
        read_data = h5py.File(r_data.fn, 'r')
    except IOError:
        # probably truncated file
        return None
    events_slot = '/'.join((
        '/Analyses', r_data.corr_group, 'Events'))
    if events_slot not in read_data:
        return None
    read_means = read_data[events_slot]['norm_mean']

    if r_data.strand == '-':
        read_means = read_means[::-1]

    return read_means

def get_reads_base_means(chrm_strand_reads, chrm_len, rev_strand):
    base_sums = np.zeros(chrm_len)
    base_cov = np.zeros(chrm_len, dtype=np.int_)
    for r_data in chrm_strand_reads:
        # extract read means data so data across all chrms is not
        # in RAM at one time
        read_means = get_read_base_means(r_data)
        if read_means is None: continue

        base_sums[r_data.start:
                  r_data.start + len(read_means)] += read_means
        base_cov[r_data.start:r_data.start + len(read_means)] += 1

    return base_sums / base_cov

def get_base_means(raw_read_coverage, chrm_sizes):
    # ignore divide by zero errors that occur where there is no
    # coverage. Need to correct nan values after subtracting two sets of
    # coverage so leave as nan for now
    old_err_settings = np.seterr(all='ignore')
    # take the mean over all signal overlapping each base
    mean_base_signal = {}
    for chrm, strand in [(c, s) for c in chrm_sizes.keys()
                         for s in ('+', '-')]:
        if (chrm, strand) in raw_read_coverage:
            cs_base_means = get_reads_base_means(
                raw_read_coverage[(chrm, strand)], chrm_sizes[chrm],
                strand == '-')
        else:
            cs_base_means = np.empty(chrm_sizes[chrm])
            cs_base_means[:] = np.nan
        mean_base_signal[(chrm, strand)] = cs_base_means
    _ = np.seterr(**old_err_settings)

    return mean_base_signal

def get_reads_base_sds(chrm_strand_reads, chrm_len, rev_strand):
    base_sd_sums = np.zeros(chrm_len)
    base_cov = np.zeros(chrm_len, dtype=np.int_)
    for r_data in chrm_strand_reads:
        # extract read means data so data across all chrms is not
        # in RAM at one time
        try:
            read_data = h5py.File(r_data.fn, 'r')
        except IOError:
            # probably truncated file
            continue
        events_slot = '/'.join((
            '/Analyses', r_data.corr_group, 'Events'))
        if events_slot not in read_data:
            continue
        read_sds = read_data[events_slot]['norm_stdev']

        if rev_strand:
            read_sds = read_sds[::-1]
        base_sd_sums[r_data.start:
                     r_data.start + len(read_sds)] += read_sds
        base_cov[r_data.start:r_data.start + len(read_sds)] += 1

    return base_sd_sums / base_cov

def get_base_sds(raw_read_coverage, chrm_sizes):
    # ignore divide by zero errors that occur where there is no
    # coverage. Need to correct nan values after subtracting two sets of
    # coverage so leave as nan for now
    old_err_settings = np.seterr(all='ignore')
    # take the mean over all signal overlapping each base
    mean_base_sds = {}
    for chrm, strand in [(c, s) for c in chrm_sizes.keys()
                         for s in ('+', '-')]:
        mean_base_sds[(chrm, strand)] = get_reads_base_sds(
            raw_read_coverage[(chrm, strand)], chrm_sizes[chrm],
            strand == '-')
    _ = np.seterr(**old_err_settings)

    return mean_base_sds

def get_reads_base_lengths(chrm_strand_reads, chrm_len, rev_strand):
    base_length_sums = np.zeros(chrm_len)
    base_cov = np.zeros(chrm_len, dtype=np.int_)
    for r_data in chrm_strand_reads:
        # extract read means data so data across all chrms is not
        # in RAM at one time
        try:
            read_data = h5py.File(r_data.fn, 'r')
        except IOError:
            # probably truncated file
            continue
        events_slot = '/'.join((
            '/Analyses', r_data.corr_group, 'Events'))
        if events_slot not in read_data:
            continue
        read_lengths = read_data[events_slot]['length']

        if rev_strand:
            read_lengths = read_lengths[::-1]
        base_length_sums[
            r_data.start:
            r_data.start + len(read_lengths)] += read_lengths
        base_cov[r_data.start:r_data.start + len(read_lengths)] += 1

    return base_length_sums / base_cov

def get_base_lengths(raw_read_coverage, chrm_sizes):
    # ignore divide by zero errors that occur where there is no
    # coverage. Need to correct nan values after subtracting two sets of
    # coverage so leave as nan for now
    old_err_settings = np.seterr(all='ignore')
    # take the mean over all signal overlapping each base
    mean_base_lengths = {}
    for chrm, strand in [(c, s) for c in chrm_sizes.keys()
                         for s in ('+', '-')]:
        mean_base_lengths[(chrm, strand)] = get_reads_base_lengths(
            raw_read_coverage[(chrm, strand)], chrm_sizes[chrm],
            strand == '-')
    _ = np.seterr(**old_err_settings)

    return mean_base_lengths

def get_reads_events(chrm_strand_reads, rev_strand):
    # note that this function assumes that all reads come from the same
    # chromosome and strand
    chrm_strand_base_means = []
    for r_data in chrm_strand_reads:
        # extract read means data so data across all chrms is not
        # in RAM at one time
        try:
            read_data = h5py.File(r_data.fn, 'r')
        except IOError:
            # probably truncated file
            continue
        events_slot = '/'.join((
            '/Analyses', r_data.corr_group, 'Events'))
        if events_slot not in read_data:
            continue
        read_means = read_data[events_slot]['norm_mean']
        if rev_strand:
            read_means = read_means[::-1]
        chrm_strand_base_means.append((
            read_means, r_data.start, r_data.end))

    if len(chrm_strand_base_means) == 0: return None

    chrm_signal = np.concatenate(zip(*chrm_strand_base_means)[0])
    chrm_pos = np.concatenate(
        [np.arange(r_data[1], r_data[2])
         for r_data in chrm_strand_base_means])
    # get order of all bases from position array
    as_chrm_pos = np.argsort(chrm_pos)
    # then sort the signal array by genomic position and
    # split into event means by base
    chrm_strand_base_events = dict(zip(
        np.unique(chrm_pos[as_chrm_pos]),
        np.split(chrm_signal[as_chrm_pos], np.where(
            np.concatenate([[0,], np.diff(
                chrm_pos[as_chrm_pos])]) > 0)[0])))

    return chrm_strand_base_events

if __name__ == '__main__':
    raise NotImplementedError, (
        'This is a module. See commands with `nanoraw -h`')
