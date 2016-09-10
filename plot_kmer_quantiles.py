import sys, os

import h5py

import numpy as np

from itertools import combinations
from collections import defaultdict

def plot_kmer_dist(files, corrected_group, kmer_val, kmer_thresh, fn_base):
    all_raw_data = []
    for fn in files:
        read_data = h5py.File(fn)
        if 'Analyses/RawGenomeCorrected' not in read_data:
            continue
        seq = ''.join(read_data['Analyses/' + corrected_group +
                                '/template/Events']['base'])
        means = np.array(read_data['Analyses/' + corrected_group +
                                   '/template/Events']['mean'])
        all_raw_data.append((seq, means))

    def get_mean_sd(norm_quantiles):
        all_trimers = defaultdict(list)
        for seq, means in all_raw_data:
            read_robust_med = np.percentile(means, norm_quantiles).mean()
            #read_mad = np.median(np.abs(means - read_robust_med))
            #norm_means = (means - read_robust_med) / read_mad
            norm_means = means - read_robust_med

            read_trimers = defaultdict(list)
            for trimer, event_mean in zip(
                    [''.join(bs) for bs in zip(
                        seq[:-2], seq[1:-1], seq[2:])],
                    norm_means[2:]):
                read_trimers[trimer].append(event_mean)
            if min(len(x) for x in read_trimers.values()) > kmer_thresh:
                for trimer, event_means in read_trimers.items():
                    all_trimers[trimer].append(np.mean(event_mean))

    return

def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(
        description='Plot raw signal corrected with correct_raw.' )
    parser.add_argument('fast5_basedir',
                        help='Directory containing fast5 files.')

    parser.add_argument('--corrected-group', default='RawGenomeCorrected_000',
                        help='FAST5 group to plot created by correct_raw ' +
                        'script. Default: %(default)s')

    parser.add_argument('--kmer-value', default=3, type=int, options=(2,3,4),
                        help='Value of K to analyze. Should be one of ' +
                        '{2,3,4}. Default: %(default)d')
    parser.add_argument('--num-trimer-threshold', default=4, type=int,
                        help='Number of each kmer required to include ' +
                        'a read in read level averages. Default: %(default)d')

    parser.add_argument('--pdf-filebase', default='Nanopore_read_coverage',
                        help='Base for PDF to store plots (suffix depends ' +
                        'on plot type to avoid overwriting plots). ' +
                        'Default: %(default)s')
    args = parser.parse_args()

    return (args.fast5_basedir, args.corrected_group,
            args.kmer_value, args.num_trimers_thresh, args.pdf_filebase)

def main():
    (filebase, corrected_group, kmer_val,
     kmer_thresh, fn_base) = parse_arguments()

    files = [os.path.join(filebase, fn) for fn in os.listdir(filebase)]
    plot_kmer_dist(files, corrected_group, kmer_val, kmer_thresh, fn_base)

    return

if __name__ == '__main__':
    main()
