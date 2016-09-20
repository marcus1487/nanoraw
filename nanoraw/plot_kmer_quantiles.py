import sys, os

import h5py

import numpy as np

from itertools import repeat
from collections import defaultdict

VERBOSE = False

try:
    import rpy2.robjects as r
    from rpy2.robjects.packages import importr
    ggplot = importr("ggplot2")
    r.r('''
plotKmerDist <- function(dat){
    print(ggplot(dat) + geom_boxplot(aes(x=Trimer, y=Signal, color=Base)) +
    theme_bw() + theme(axis.text.x=element_text(angle=60, hjust=1, size=8)) +
    scale_color_manual(values=c('#00CC00', '#0000CC', '#FFB300', '#CC0000')))
}
    ''')
    plotKmerDist = r.globalenv['plotKmerDist']
    r.r('''
    plotKmerDistWReadPath <- function(dat){
print(ggplot(dat) + geom_boxplot(aes(x=Trimer, y=Signal, color=Base)) +
    theme_bw() + theme(axis.text.x=element_text(angle=60, hjust=1, size=8)) +
    scale_color_manual(values=c('#00CC00', '#0000CC', '#FFB300', '#CC0000')))
print(ggplot(dat) + geom_path(aes(x=Trimer, y=Signal, group=Read), alpha=0.1) +
    theme_bw() + theme(axis.text.x=element_text(angle=60, hjust=1, size=8)) +
    scale_color_manual(values=c('#00CC00', '#0000CC', '#FFB300', '#CC0000')))
}
    ''')
    plotKmerDistWReadPath = r.globalenv['plotKmerDistWReadPath']
except:
    sys.stderr.write(
        '*' * 60 + '\nERROR: Must have rpy2, R and ' +
        'R package ggplot2 installed in order to plot.\n' +
        '*' * 60 + '\n\n')
    raise


def plot_kmer_dist(files, corrected_group, read_mean, kmer_len,
                   kmer_thresh, fn_base):
    all_raw_data = []
    for fn in files:
        read_data = h5py.File(fn)
        if 'Analyses/' + corrected_group not in read_data:
            continue
        seq = ''.join(read_data['Analyses/' + corrected_group +
                                '/template/Events']['base'])
        means = np.array(read_data['Analyses/' + corrected_group +
                                   '/template/Events']['norm_mean'])
        all_raw_data.append((seq, means))

    all_trimers = defaultdict(list)
    for read_i, (seq, means) in enumerate(all_raw_data):
        read_trimers = defaultdict(list)
        for trimer, event_mean in zip(
                [''.join(bs) for bs in zip(*[
                    seq[i:] for i in range(kmer_len)])],
                means[kmer_len - 1:]):
            read_trimers[trimer].append(event_mean)
        if min(len(x) for x in read_trimers.values()) > kmer_thresh:
            for trimer, trimer_means in read_trimers.items():
                if read_mean:
                    all_trimers[trimer].append((np.mean(trimer_means), read_i))
                else:
                    all_trimers[trimer].extend(zip(trimer_means, repeat(read_i)))

    kmer_levels = [kmer for means, kmer in sorted([
        (np.mean(zip(*means)[0]), kmer)
        for kmer, means in all_trimers.items()])]

    plot_data = [
        (kmer, kmer[-1], sig_mean, read_i)
        for kmer in kmer_levels
        for sig_mean, read_i in all_trimers[kmer]]

    trimerDat = r.DataFrame({
        'Trimer':r.FactorVector(
            r.StrVector(zip(*plot_data)[0]),
            ordered=True, levels=r.StrVector(kmer_levels)),
        'Base':r.StrVector(zip(*plot_data)[1]),
        'Signal':r.FloatVector(zip(*plot_data)[2]),
        'Read':r.StrVector(zip(*plot_data)[3])})
    # code to plot kmers as tile of colors but adds gridExtra dependency
    if False:
        kmer_plot_data = [
            (kmer_i, pos_i, base) for kmer_i, kmer in enumerate(kmer_leves)
            for pos_i, base in enumerate(kmer)]
        kmerDat = r.DataFrame({
            'Kmer':r.IntVector(zip(*kmer_plot_data)[0]),
            'Position':r.IntVector(zip(*kmer_plot_data)[1]),
            'Base':r.StrVector(zip(*kmer_plot_data)[2])})

    if read_mean:
        r.r('pdf("' + fn_base + '.read_mean.pdf", height=7, width=10)')
        plotKmerDistWReadPath(trimerDat)
        r.r('dev.off()')
    else:
        r.r('pdf("' + fn_base + '.pdf", height=7, width=10)')
        plotKmerDist(trimerDat)
        r.r('dev.off()')

    return

def main(args):
    files = [os.path.join(args.fast5_basedir, fn)
             for fn in os.listdir(args.fast5_basedir)]
    plot_kmer_dist(
        files, args.corrected_group, args.read_mean, args.kmer_length,
        args.num_trimer_threshold, args.pdf_filebase)

    return

# define function for getting parser so it can be shared in
# __main__ package script
def get_parser(with_help=True):
    import argparse
    parser = argparse.ArgumentParser(
        description='Plot distribution of signal across kmers.',
        add_help=with_help)
    parser.add_argument(
        'fast5_basedir',
        help='Directory containing fast5 files.')

    parser.add_argument(
        '--corrected-group', default='RawGenomeCorrected_000',
        help='FAST5 group to plot created by correct_raw ' +
        'script. Default: %(default)s')

    parser.add_argument(
        '--kmer-length', default=3, type=int, choices=(2,3,4),
        help='Value of K to analyze. Should be one of ' +
        '{2,3,4}. Default: %(default)d')
    parser.add_argument(
        '--num-trimer-threshold', default=4, type=int,
        help='Number of each kmer required to include ' +
        'a read in read level averages. Default: %(default)d')

    parser.add_argument(
        '--read-mean', default=False, action='store_true',
        help='Plot kmer event means across reads as opposed ' +
        'to each event.')

    parser.add_argument(
        '--pdf-filebase', default='Nanopore_kmer_distribution',
        help='Base for PDF to store plots (suffix depends on ' +
        'plot type to avoid overwriting plots). Default: %(default)s')

    return parser

def args_and_main():
    main(get_parser().parse_args())
    return

if __name__ == '__main__':
    args_and_main()
