import sys, os

import numpy as np

from nanoraw_helper import (
    parse_fast5s, parse_obs_filter, filter_reads)

VERBOSE = False

def write_wiggle(files, corrected_group, wiggle_base, basecall_subgroups,
                 obs_filter):
    if VERBOSE: sys.stderr.write('Parsing files.\n')
    raw_read_coverage = parse_fast5s(
        files, corrected_group, basecall_subgroups)
    raw_read_coverage = filter_reads(raw_read_coverage, obs_filter)

    if VERBOSE: sys.stderr.write('Calculating read coverage.\n')
    wiggle_cov = []
    for (chrom, strand), reads_data in raw_read_coverage.items():
        max_end = max(r_data.end for r_data in reads_data)
        chrom_coverage = np.zeros(max_end, dtype=np.int_)
        for r_data in reads_data:
            chrom_coverage[r_data.start:r_data.end] += 1
        wiggle_cov.append((chrom, strand, chrom_coverage))

    if VERBOSE: sys.stderr.write('Writing wiggle.\n')
    plus_wig_fp = open(wiggle_base + '.plus.wig', 'w')
    minus_wig_fp = open(wiggle_base + '.minus.wig', 'w')
    plus_wig_fp.write(
        'track type=wiggle_0 name="{0} Forward Strand" description={0}\n'.format(
            wiggle_base))
    minus_wig_fp.write(
        'track type=wiggle_0 name="{0} Reverse Strand" description={0}\n'.format(
            wiggle_base))
    for chrm, strand, chrm_cov in wiggle_cov:
        wig_fp = plus_wig_fp if strand == '+' else minus_wig_fp
        wig_fp.write("variableStep chrom={} span=1\n".format(
            chrm))
        wig_fp.write('\n'.join([
            str(int(pos) + 1) + " " + str(int(val))
            for pos, val in enumerate(chrm_cov) if val > 0]) +
                     '\n')

    plus_wig_fp.close()
    minus_wig_fp.close()

    return

def wiggle_main(args):
    global VERBOSE
    VERBOSE = not args.quiet

    files = [os.path.join(base_dir, fn)
             for base_dir in args.fast5_basedirs
             for fn in os.listdir(base_dir)]
    write_wiggle(files, args.corrected_group, args.wiggle_basename,
                 args.basecall_subgroups,
                 parse_obs_filter(args.obs_per_base_filter))

    return


if __name__ == '__main__':
    raise NotImplementedError, (
        'This is a module. See commands with `nanoraw -h`')
