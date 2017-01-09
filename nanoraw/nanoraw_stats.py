import sys, os

import numpy as np

from scipy import stats
from collections import defaultdict

def correct_multiple_testing(pvals):
    """ Use FDR Benjamini-Hochberg multiple testing correction
    """
    pvals = np.asarray(pvals)

    pvals_sortind = np.argsort(pvals)
    pvals_sorted = pvals[pvals_sortind]
    sortrevind = pvals_sortind.argsort()

    pvals_corrected_raw = pvals_sorted / (np.arange(
        1,len(pvals)+1)/float(len(pvals)))
    pvals_corrected = np.minimum.accumulate(
        pvals_corrected_raw[::-1])[::-1]
    pvals_corrected[pvals_corrected>1] = 1

    return pvals_corrected[sortrevind]

def _calc_fm_pval(pvals):
    return 1.0 - stats.chi2.cdf(
        np.sum(np.log(pvals)) * -2,
        pvals.shape[0] * 2)

def calc_fishers_method(pos_pvals, offset):
    pvals_np = np.empty(pos_pvals[-1][1] + 1)
    pvals_np[:] = np.NAN
    pvals_np[[list(zip(*pos_pvals)[1])]] = zip(*pos_pvals)[0]

    fishers_pvals = [
        (_calc_fm_pval(pvals_np[pos - offset:pos + offset + 1]),
         pos, cov1, cov2)
        for _, pos, cov1, cov2 in pos_pvals
        if pos - offset >= 0 and pos + offset + 1 <= pvals_np.shape[0]
        and not np.any(np.isnan(
            pvals_np[pos - offset:pos + offset + 1]))]

    return fishers_pvals

def mann_whitney_u_test(samp1, samp2):
    s1_len = samp1.shape[0]
    s2_len = samp2.shape[0]
    tot_len = s1_len + s2_len

    all_vals = np.concatenate([samp1, samp2])
    ranks = np.empty(tot_len, int)
    ranks[all_vals.argsort()] = np.arange(1, tot_len + 1)
    s1_ranks_sum = ranks[:s1_len].sum()
    #s2_ranks_sum = ranks[s1_len:].sum()

    u1 = s1_ranks_sum - (s1_len * (s1_len + 1)) / 2
    #u2 = s2_ranks_sum - (s2_len * (s2_len + 1)) / 2

    mu = s1_len * s2_len / 2
    rhou = np.sqrt(s1_len * s2_len * (s1_len + s2_len + 1) / 12)

    z = np.abs(u1 - mu) / rhou

    pval = stats.norm.cdf(-z) * 2.0

    return pval

def get_all_significance(
        base_events1, base_events2, test_type, min_test_vals,
        all_stats_fn, fishers_method_offset):
    if VERBOSE: sys.stderr.write(
            'Test significance of difference in base signal.\n')
    # get num_region most significantly different regions from
    # each chrm then find global most signif after
    position_pvals = []
    for chrm_strand in set(base_events1).intersection(base_events2):
        chrm, strand = chrm_strand
        if test_type == 'ttest':
            chrm_pvals = [
                (np.abs(stats.ttest_ind(
                    base_events1[chrm_strand][pos],
                    base_events2[chrm_strand][pos])[1]), pos,
                base_events1[chrm_strand][pos].shape[0],
                base_events2[chrm_strand][pos].shape[0])
                for pos in sorted(set(
                    base_events1[chrm_strand]).intersection(
                    base_events2[chrm_strand]))
                if min(base_events1[chrm_strand][pos].shape[0],
                       base_events2[chrm_strand][pos].shape[0])
                >= min_test_vals]
        elif test_type == 'mw_utest':
            # store z-scores from u-test
            chrm_pvals = [
                (mann_whitney_u_test(
                    base_events1[chrm_strand][pos],
                    base_events2[chrm_strand][pos]), pos,
                base_events1[chrm_strand][pos].shape[0],
                base_events2[chrm_strand][pos].shape[0])
                for pos in sorted(set(
                    base_events1[chrm_strand]).intersection(
                        base_events2[chrm_strand]))
                if min(base_events1[chrm_strand][pos].shape[0],
                       base_events2[chrm_strand][pos].shape[0])
                >= min_test_vals]
        else:
            raise RuntimeError, ('Invalid significance test type ' +
                                 'provided: ' + str(test_type))

        if len(chrm_pvals) == 0: continue
        if fishers_method_offset is not None:
            chrm_pvals = calc_fishers_method(
                chrm_pvals, fishers_method_offset)

        position_pvals.extend(
            (pval, pos, chrm, strand, cov1, cov2)
            for pval, pos, cov1, cov2 in chrm_pvals)

    if len(position_pvals) == 0:
        sys.stderr.write(
            '*' * 60 + '\nERROR: No regions contain minimum ' +
            'number of reads.\n' + '*' * 60 + '\n')
        sys.exit()

    position_pvals = sorted(position_pvals)
    fdr_corr_pvals = correct_multiple_testing(zip(*position_pvals)[0])
    all_stats = [(pval, qval, pos, chrm, strand, cov1, cov2)
                 for qval, (pval, pos, chrm, strand, cov1, cov2) in
                 zip(fdr_corr_pvals, position_pvals)]

    if all_stats_fn is not None:
        chrm_strand_stats = defaultdict(list)
        for pval, qval, pos, chrm, strand, cov1, cov2 in all_stats:
            chrm_strand_stats[(chrm, strand)].append((
                pos, pval, qval, cov1, cov2))
        with open(all_stats_fn, 'w') as stats_fp:
            for (chrm, strand), pos_stats in chrm_strand_stats.items():
                stats_fp.write('>>>>::' + chrm + '::' + strand + '\n')
                stats_fp.write('\n'.join([
                    '{:d}\t{:.2g}\t{:.2g}\t{:d}\t{:d}'.format(
                        pos, pval, qval, cov1, cov2)
                    for pos, pval, qval, cov1, cov2 in
                    sorted(pos_stats)]) + '\n')

    return all_stats

def parse_stats(stats_fn):
    all_stats = []
    with open(stats_fn) as stats_fp:
        curr_chrm, curr_strand = None, None
        try:
            for line in stats_fp:
                if line.startswith('>>>>'):
                    _, curr_chrm, curr_strand = line.strip().split("::")
                else:
                    if curr_chrm is None or curr_strand is None:
                        sys.stderr.write(
                            'WARNING: Incorrectly formatted ' +
                            'statistics file. No chrm or strand ' +
                            'before statistics lines\n')
                    pos, pval, qval, cov1, cov2 = line.split()
                    all_stats.append((
                        float(pval), float(qval), int(pos),
                        curr_chrm, curr_strand, int(cov1), int(cov2)))
        except ValueError:
            sys.stderr.write(
                '*' * 60  + '\nERROR: Attempt to load statistics ' +
                'file failed. May be an old version of statistics ' +
                'file. Try deleting statistics file and ' +
                'recalculating using current nanoraw version.\n' +
                '*' * 60 + '\n')
            sys.exit()

    return sorted(all_stats)
