import sys, os

import numpy as np

from nanoraw_helper import parse_fast5s

VERBOSE = False

def write_wiggle(files, corrected_group, wiggle_fn, basecall_subgroups):
    if VERBOSE: sys.stderr.write('Parsing files.\n')
    raw_read_coverage = parse_fast5s(
        files, corrected_group, basecall_subgroups)

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

    files = [os.path.join(base_dir, fn)
             for base_dir in args.fast5_basedirs
             for fn in os.listdir(base_dir)]
    write_wiggle(files, args.corrected_group, args.wiggle_filename,
                 args.basecall_subgroups)

    return


if __name__ == '__main__':
    raise NotImplementedError, (
        'This is a module. Run with `nanoraw plot_signal -h`')
