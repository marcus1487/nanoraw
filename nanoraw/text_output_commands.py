import sys, os

import numpy as np

from collections import defaultdict

import nanoraw_stats as ns
import nanoraw_helper as nh

VERBOSE = False

SMALLEST_PVAL=1e-15
WIG_HEADER='track type=wiggle_0 name="{0}_{1}_{2}{3}" ' + \
    'description="{0} {1} {2}{4}"\n'
GROUP2_NAME='group2'

def write_wiggle(wig_base, group_text, data_values, type_name,
                 filter_zeros=False):
    group_w_dot = '' if group_text == '' else '.' + group_text
    group_w_us = '' if group_text == '' else '_' + group_text
    group_w_space = '' if group_text == '' else ' ' + group_text
    plus_wig_fp = open(
        wig_base + '.' + type_name + group_w_dot + '.plus.wig', 'w')
    minus_wig_fp = open(
        wig_base + '.' + type_name + group_w_dot + '.minus.wig', 'w')
    plus_wig_fp.write(WIG_HEADER.format(
        wig_base, type_name, 'fwd_strand', group_w_us, group_w_space))
    minus_wig_fp.write(WIG_HEADER.format(
        wig_base, type_name, 'rev_strand', group_w_us, group_w_space))
    for (chrm, strand), chrm_values in data_values.iteritems():
        wig_fp = plus_wig_fp if strand == '+' else minus_wig_fp
        wig_fp.write("variableStep chrom={} span=1\n".format(chrm))
        wig_fp.write('\n'.join([
            str(int(pos) + 1) + " " + str(round(val, 4))
            for pos, val in enumerate(chrm_values)
            if not (np.isnan(val) or (
                    filter_zeros and np.equal(val, 0.0)))]) + '\n')

    plus_wig_fp.close()
    minus_wig_fp.close()

    return

def write_pvals_and_qvals_wig(
        all_stats, wig_base, write_pvals, write_qvals):
    if VERBOSE: sys.stderr.write('Parsing statistics.\n')
    raw_chrm_strand_stats = defaultdict(list)
    for (pval_f, qval_f, pval, qval, pos, chrm, strand,
         cov1, cov2) in all_stats:
        raw_chrm_strand_stats[(chrm, strand)].append((pos, pval, qval))

    chrm_strand_pvals = {}
    chrm_strand_qvals = {}
    for chrm_strand, stats in raw_chrm_strand_stats.iteritems():
        chrm_poss = zip(*stats)[0]
        raw_chrm_pvals = zip(*stats)[1]
        raw_chrm_qvals = zip(*stats)[2]
        max_pos = max(chrm_poss)

        # arrange and store p-values
        chrm_pvals = np.empty(max_pos + 1)
        chrm_pvals[:] = np.nan
        np.put(chrm_pvals, chrm_poss, raw_chrm_pvals)
        chrm_strand_pvals[chrm_strand] = -np.log10(np.maximum(
            SMALLEST_PVAL, chrm_pvals))

        # arrange and store q-values
        chrm_qvals = np.empty(max_pos + 1)
        chrm_qvals[:] = np.nan
        np.put(chrm_qvals, chrm_poss, raw_chrm_qvals)
        chrm_strand_qvals[chrm_strand] = -np.log10(np.maximum(
            SMALLEST_PVAL, chrm_qvals))

    if VERBOSE: sys.stderr.write('Writing statistics wig(s).\n')
    if write_pvals:
        write_wiggle(wig_base, '', chrm_strand_pvals, 'neg_log10_pvals')
    if write_qvals:
        write_wiggle(wig_base, '', chrm_strand_qvals, 'neg_log10_qvals')

    return

def get_chrm_sizes(raw_read_coverage, raw_read_coverage2=None):
    strand_chrm_sizes = defaultdict(list)
    for (chrm, strand), cs_read_cov in \
      raw_read_coverage.iteritems():
        strand_chrm_sizes[chrm].append(max(
          r_data.end for r_data in cs_read_cov))
    if raw_read_coverage2 is not None:
        for (chrm, strand), cs_read_cov in \
          raw_read_coverage2.iteritems():
            strand_chrm_sizes[chrm].append(max(
              r_data.end for r_data in cs_read_cov))

    return dict(
        (chrm, max(strnd_sizes))
        for chrm, strnd_sizes in
        strand_chrm_sizes.iteritems())

def write_length_wig(
        raw_read_coverage, chrm_sizes, wig_base, group_name):
    if VERBOSE: sys.stderr.write('Parsing events lengths.\n')
    base_lens = nh.get_base_lengths(raw_read_coverage, chrm_sizes)

    if VERBOSE: sys.stderr.write('Writing length wig.\n')
    write_wiggle(wig_base, group_name, base_lens, 'length')

    return

def write_signal_sd_wig(
        raw_read_coverage, chrm_sizes, wig_base, group_name):
    if VERBOSE: sys.stderr.write('Parsing signal SDs.\n')
    base_sds = nh.get_base_sds(raw_read_coverage, chrm_sizes)

    if VERBOSE: sys.stderr.write('Writing signal SD wig.\n')
    write_wiggle(wig_base, group_name, base_sds, 'signalSd')

    return

def write_signal_and_diff_wigs(
        raw_read_coverage1, raw_read_coverage2, chrm_sizes,
        wig_base, group1_name, write_sig, write_diff):
    if VERBOSE: sys.stderr.write('Parsing mean base signals.\n')
    base_means1 = nh.get_base_means(raw_read_coverage1, chrm_sizes)
    if raw_read_coverage2 is not None:
        base_means2 = nh.get_base_means(raw_read_coverage2, chrm_sizes)

        if write_diff:
            if VERBOSE: sys.stderr.write(
                    'Calculating signal differences.\n')
            sig_diffs = {}
            for chrm, strand in [(c, s) for c in chrm_sizes.keys()
                                 for s in ('+', '-')]:
                # calculate difference and set no coverage
                # (nan) values to zero
                sig_diffs[(chrm, strand)] \
                    = base_means1[(chrm, strand)] - \
                    base_means2[(chrm, strand)]
            if VERBOSE: sys.stderr.write('Writing differnce wig.\n')
            write_wiggle(wig_base, '', sig_diffs, 'difference')
        if write_sig:
            if VERBOSE: sys.stderr.write('Writing signal wigs.\n')
            write_wiggle(wig_base, GROUP2_NAME, base_means2, 'signal')

    if write_sig:
        write_wiggle(wig_base, group1_name, base_means1, 'signal')

    return

def write_cov_wig(raw_read_coverage, wig_base, group_text):
    read_coverage = nh.get_coverage(raw_read_coverage)

    if VERBOSE: sys.stderr.write('Writing coverage wig.\n')
    write_wiggle(wig_base, group_text, read_coverage, 'coverage', True)

    return

def write_all_wiggles(
        files1, files2, corrected_group, basecall_subgroups, obs_filter,
        test_type, min_test_vals, stats_fn, fishers_method_offset,
        wig_base, wig_types):
    stats_file_exists = stats_fn is not None and os.path.isfile(stats_fn)
    include_stats = 'pvals' in wig_types or 'qvals' in wig_types
    if include_stats and stats_file_exists:
        if VERBOSE: sys.stderr.write('Loading statistics from file.\n')
        all_stats = ns.parse_stats(stats_fn)

    if VERBOSE: sys.stderr.write('Parsing FAST5 files.\n')
    raw_read_coverage1 = nh.parse_fast5s(
        files1, corrected_group, basecall_subgroups)
    raw_read_coverage1 = nh.filter_reads(raw_read_coverage1, obs_filter)

    group1_name = '' if files2 is None else 'group1'
    if files2 is not None:
        raw_read_coverage2 = nh.parse_fast5s(
            files2, corrected_group, basecall_subgroups)
        raw_read_coverage2 = nh.filter_reads(
            raw_read_coverage2, obs_filter)
        chrm_sizes = get_chrm_sizes(
            raw_read_coverage1, raw_read_coverage2)

        if include_stats and not stats_file_exists:
            if VERBOSE: sys.stderr.write('Calculating statistics.\n')
            all_stats = ns.get_all_significance(
                raw_read_coverage1, raw_read_coverage2, test_type,
                min_test_vals, stats_fn, fishers_method_offset)

        if VERBOSE: sys.stderr.write('Writing wiggles.\n')
        if 'coverage' in wig_types:
            write_cov_wig(raw_read_coverage2, wig_base, GROUP2_NAME)
        if 'signal_sd' in wig_types:
            write_signal_sd_wig(
                raw_read_coverage2, chrm_sizes, wig_base, GROUP2_NAME)
        if 'length' in wig_types:
            write_length_wig(raw_read_coverage2, chrm_sizes,
                             wig_base, GROUP2_NAME)

        # need to do signal and difference call once either with or
        # w/o second set of files (unlike coverage, sds and length
        if 'signal' in wig_types or 'difference' in wig_types:
            write_signal_and_diff_wigs(
                raw_read_coverage1, raw_read_coverage2, chrm_sizes,
                wig_base, group1_name, 'signal' in wig_types,
                'difference' in wig_types)
    else:
        chrm_sizes = get_chrm_sizes(raw_read_coverage1)
        if VERBOSE: sys.stderr.write('Writing wiggles.\n')
        if 'signal' in wig_types:
            write_signal_and_diff_wigs(
                raw_read_coverage1, None, chrm_sizes, wig_base,
                group1_name, 'signal' in wig_types, False)

    if 'coverage' in wig_types:
        write_cov_wig(raw_read_coverage1, wig_base, group1_name)
    if 'signal_sd' in wig_types:
        write_signal_sd_wig(raw_read_coverage1, chrm_sizes,
                            wig_base, group1_name)
    if 'length' in wig_types:
        write_length_wig(raw_read_coverage1, chrm_sizes,
                         wig_base, group1_name)
    if 'pvals' in wig_types or 'qvals' in wig_types:
        write_pvals_and_qvals_wig(
            all_stats, wig_base, 'pvals' in wig_types,
            'qvals' in wig_types)

    return

def write_most_signif(
        files1, files2, num_regions, qval_thresh, corrected_group,
        basecall_subgroups, seqs_fn, num_bases, test_type, obs_filter,
        min_test_vals, stats_fn, fasta_fn, fishers_method_offset):
    calc_stats = stats_fn is None or not os.path.isfile(stats_fn)
    if not calc_stats:
        if VERBOSE: sys.stderr.write('Loading statistics from file.\n')
        all_stats = ns.parse_stats(stats_fn)

    if calc_stats or fasta_fn is None:
        if VERBOSE: sys.stderr.write('Parsing files.\n')
        raw_read_coverage1 = nh.parse_fast5s(
            files1, corrected_group, basecall_subgroups)
        raw_read_coverage2 = nh.parse_fast5s(
            files2, corrected_group, basecall_subgroups)
        raw_read_coverage1 = nh.filter_reads(
            raw_read_coverage1, obs_filter)
        raw_read_coverage2 = nh.filter_reads(
            raw_read_coverage2, obs_filter)

    if calc_stats:
        if VERBOSE: sys.stderr.write('Calculating statistics.\n')
        all_stats = ns.get_all_significance(
            raw_read_coverage1, raw_read_coverage2, test_type,
            min_test_vals, stats_fn, fishers_method_offset)

    plot_intervals = ns.get_most_signif_regions(
        all_stats, num_bases, num_regions, qval_thresh)
    if fasta_fn is None:
        reg_seqs = get_region_sequences(
            plot_intervals, raw_read_coverage1, raw_read_coverage2,
            num_bases, corrected_group)
    else:
        fasta_records = nh.parse_fasta(fasta_fn)
        reg_seqs = [
            (p_int, fasta_records[chrm][start:start+num_bases])
            for p_int, (chrm, start, strand, reg_name)
            in plot_intervals if chrm in fasta_records]

    # get reads overlapping each region
    if VERBOSE: sys.stderr.write('Outputting region seqeuences.\n')
    with open(seqs_fn, 'w') as seqs_fp:
        for reg_i, reg_seq in reg_seqs:
            chrm, start, strand, stat = next(
                p_int for p_reg_i, p_int in plot_intervals
                if p_reg_i == reg_i)
            if strand == '-':
                reg_seq = nh.rev_comp(reg_seq)
            seqs_fp.write('>{0}::{1:d}::{2} {3}\n{4}\n'.format(
                chrm, start, strand, stat, ''.join(reg_seq)))

    return


def wiggle_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    nh.VERBOSE = VERBOSE
    ns.VERBOSE = VERBOSE

    if (any(data_type in args.wiggle_types
            for data_type in ['pvals', 'qvals']) and
        args.fast5_basedirs2 is None and
        args.statistics_filename is None):
        sys.stderr.write(
            '*' * 60 + '\nERROR: Must provide either two sets of ' +
            'FAST5s or a statistics filename to output ' +
            'pvals and/or qvals wiggle files.\n' + '*' * 60 + '\n')
        sys.exit()
    if ('difference' in args.wiggle_types and
        args.fast5_basedirs2 is None):
        sys.stderr.write(
            '*' * 60 + '\nERROR: Must provide two sets of FAST5s ' + \
            'to output difference wiggle files.\n' + '*' * 60 + '\n')
        sys.exit()

    files1, files2 = nh.get_files_lists(
        args.fast5_basedirs, args.fast5_basedirs2)
    write_all_wiggles(
        files1, files2, args.corrected_group, args.basecall_subgroups,
        nh.parse_obs_filter(args.obs_per_base_filter),
        args.test_type, args.minimum_test_reads,
        args.statistics_filename, args.fishers_method_offset,
        args.wiggle_basename, args.wiggle_types)

    return

def write_signif_diff_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    nh.VERBOSE = VERBOSE
    ns.VERBOSE = VERBOSE

    files1, files2 = nh.get_files_lists(
        args.fast5_basedirs, args.fast5_basedirs2)

    write_most_signif(
        files1, files2, args.num_regions, args.q_value_threshold,
        args.corrected_group, args.basecall_subgroups,
        args.sequences_filename, args.num_bases, args.test_type,
        nh.parse_obs_filter(args.obs_per_base_filter),
        args.minimum_test_reads, args.statistics_filename,
        args.genome_fasta, args.fishers_method_offset)

    return


if __name__ == '__main__':
    raise NotImplementedError, (
        'This is a module. See commands with `nanoraw -h`')
