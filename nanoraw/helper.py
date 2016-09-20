import numpy as np

from itertools import izip

NORM_TYPES = ('none', 'ont', 'median', 'robust_median')

robust_quantiles = (46.5, 53.5)

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


"""
    parser.add_argument('--normalization-type', default='median',
                        choices=['raw', 'ont', 'median'],
                        help='Raw nanopore normalization ' +
                        'method. Should be one of (raw, ont and ' +
                        'median). Default: median')
"""
