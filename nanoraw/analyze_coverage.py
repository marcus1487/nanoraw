import sys, os

import numpy as np

from nanoraw_helper import parse_files

VERBOSE = False

def write_wiggle(files, corrected_group, wiggle_fn):
    if VERBOSE: sys.stderr.write('Parsing files.\n')
    raw_read_coverage = parse_files(files, corrected_group)

    if VERBOSE: sys.stderr.write('Calculating read coverage.\n')
    wiggle_cov = []
    for chrom, reads_data in raw_read_coverage.items():
        max_end = max(r_data.end for r_data in reads_data)
        chrom_coverage = np.zeros(max_end, dtype=np.int_)
        for r_data in reads_data:
            chrom_coverage[r_data.start:r_data.end] += 1
        wiggle_cov.append((chrom, chrom_coverage))

    if VERBOSE: sys.stderr.write('Writing wiggle.\n')
    with open(wiggle_fn, 'w') as wig_fp:
        wig_fp.write(
            'track type=wiggle_0 name={0} description={0}\n'.format(
                wiggle_fn))
        for chrm, chrm_cov in wiggle_cov:
            wig_fp.write("variableStep chrom={} span=1\n".format(
                chrm))
            wig_fp.write('\n'.join([
                str(int(pos) + 1) + " " + str(int(val))
                for pos, val in enumerate(chrm_cov) if val > 0]) +
                         '\n')

    return

def wiggle_main(args):
    global VERBOSE
    VERBOSE = not args.quiet

    files = [os.path.join(args.fast5_basedir, fn)
             for fn in os.listdir(args.fast5_basedir)]
    write_wiggle(files, args.corrected_group, args.wiggle_filename)

    return

def get_wiggle_parser():
    import argparse
    parser = argparse.ArgumentParser(
        description='Plot raw signal from from two samples where ' +
        'FAST5 files were corrected with `nanoraw correct`.',
        add_help=False)
    parser.add_argument(
        'fast5_basedir',
        help='Directory containing fast5 files.')
    parser.add_argument(
        '--corrected-group', default='RawGenomeCorrected_000',
        help='FAST5 group to plot created by correct_raw ' +
        'script. Default: %(default)s')

    parser.add_argument(
        '--wiggle-filename',
        help="Output wiggle read coverage file. Note that this will " +
        "also be the track name in the def line (only available for " +
        "single FAST5 dir currently). Default: Dont't output " +
        "coverage wiggle.")

    parser.add_argument(
        '--quiet', '-q', default=False, action='store_true',
        help="Don't print status information.")

    return parser


if __name__ == '__main__':
    raise NotImplementedError, (
        'This is a module. Run with `nanoraw plot_signal -h`')
