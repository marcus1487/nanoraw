import sys, os

import numpy as np

from scipy import stats
from collections import defaultdict

from nanoraw_helper import get_reads_events

VERBOSE = False

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
        _calc_fm_pval(pvals_np[pos - offset:pos + offset + 1])
        if pos - offset >= 0 and pos + offset + 1 <= pvals_np.shape[0] and
        not np.any(np.isnan(pvals_np[pos - offset:pos + offset + 1])) else 1.0
        for _, pos, _, _ in pos_pvals]

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
        raw_read_coverage1, raw_read_coverage2, test_type, min_test_vals,
        all_stats_fn, fishers_method_offset):
    if VERBOSE: sys.stderr.write(
            'Test significance of difference in base signal.\n')
    # get num_region most significantly different regions from
    # each chrm then find global most signif after
    position_pvals = []
    for chrm_strand in set(raw_read_coverage1).intersection(
            raw_read_coverage2):
        chrm, strand = chrm_strand
        # get base events across all reads per chromosome/strand
        # so that all events aren't stored in RAM
        chrm_strand_base_events1 = get_reads_events(
            raw_read_coverage1[chrm_strand], strand == '-')
        chrm_strand_base_events2 = get_reads_events(
            raw_read_coverage2[chrm_strand], strand == '-')
        if test_type == 'ttest':
            chrm_pvals = [
                (np.abs(stats.ttest_ind(
                    chrm_strand_base_events1[pos],
                    chrm_strand_base_events2[pos])[1]), pos,
                chrm_strand_base_events1[pos].shape[0],
                chrm_strand_base_events2[pos].shape[0])
                for pos in sorted(set(
                    chrm_strand_base_events1).intersection(
                    chrm_strand_base_events2))
                if min(chrm_strand_base_events1[pos].shape[0],
                       chrm_strand_base_events2[pos].shape[0])
                >= min_test_vals]
        elif test_type == 'mw_utest':
            # store z-scores from u-test
            chrm_pvals = [
                (mann_whitney_u_test(
                    chrm_strand_base_events1[pos],
                    chrm_strand_base_events2[pos]), pos,
                chrm_strand_base_events1[pos].shape[0],
                chrm_strand_base_events2[pos].shape[0])
                for pos in sorted(set(
                    chrm_strand_base_events1).intersection(
                        chrm_strand_base_events2))
                if min(chrm_strand_base_events1[pos].shape[0],
                       chrm_strand_base_events2[pos].shape[0])
                >= min_test_vals]
        else:
            raise RuntimeError, ('Invalid significance test type ' +
                                 'provided: ' + str(test_type))

        if len(chrm_pvals) == 0: continue
        chrm_pvals_f = calc_fishers_method(
            chrm_pvals, fishers_method_offset) \
            if fishers_method_offset > 0 else zip(*chrm_pvals)[0]

        position_pvals.extend(
            (pval_f, pval, pos, chrm, strand, cov1, cov2)
            for (pval, pos, cov1, cov2), pval_f in
            zip(chrm_pvals, chrm_pvals_f))

    if len(position_pvals) == 0:
        sys.stderr.write(
            '*' * 60 + '\nERROR: No regions contain minimum ' +
            'number of reads.\n' + '*' * 60 + '\n')
        sys.exit()

    position_pvals = sorted(position_pvals)
    fdr_corr_pvals_f = correct_multiple_testing(zip(*position_pvals)[0])
    fdr_corr_pvals = correct_multiple_testing(sorted(zip(*position_pvals)[1]))
    all_stats = [(pval_f, qval_f, pval, qval, pos, chrm, strand, cov1, cov2)
                 for qval_f, qval, (pval_f, pval, pos, chrm, strand, cov1, cov2) in
                 zip(fdr_corr_pvals_f, fdr_corr_pvals, position_pvals)]

    if all_stats_fn is not None:
        chrm_strand_stats = defaultdict(list)
        for pval_f, qval_f, pval, qval, pos, chrm, strand, cov1, cov2 in all_stats:
            chrm_strand_stats[(chrm, strand)].append((
                pos, pval_f, qval_f, pval, qval, cov1, cov2))
        with open(all_stats_fn, 'w') as stats_fp:
            for (chrm, strand), pos_stats in chrm_strand_stats.items():
                stats_fp.write('>>>>::' + chrm + '::' + strand + '\n')
                stats_fp.write('\n'.join([
                    '{:d}\t{:.2g}\t{:.2g}\t{:.2g}\t{:.2g}\t{:d}\t{:d}'.format(
                        pos, pval_f, qval_f, pval, qval, cov1, cov2)
                    for pos, pval_f, qval_f, pval, qval, cov1, cov2 in
                    sorted(pos_stats)]) + '\n')

    return all_stats

def get_most_signif_regions(all_stats, num_bases, num_regions,
                            qval_thresh=None):
    # applied threshold for scores on each chromosome, so now
    # we include all here
    if qval_thresh is not None:
        num_regions = np.argmax(
            [x > qval_thresh for x in zip(*all_stats)[1]])
        if num_regions == 0:
            sys.stderr.write(
                '*' * 60 + '\nERROR: No regions identified q-value ' +
                'below thresh. Minumum q-value: {:.2g}\n'.format(
                    all_stats[0][1]) + '*' * 60 + '\n')
            sys.exit()

    plot_intervals = zip(
        ['{:03d}'.format(rn) for rn in range(num_regions)],
        [(chrm, max(pos - int(num_bases / 2.0), 0), strand,
          '(q-value:{0:.2g} p-value:{1:.2g})'.format(qval_f, pval_f))
         for pval_f, qval_f, pval, qval, pos, chrm, strand, cov1, cov2 in
         all_stats[:num_regions]])

    return plot_intervals

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
                    pos, pval_f, qval_f, pval, qval, cov1, cov2 = line.split()
                    all_stats.append((
                        float(pval_f), float(qval_f),
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


# functions for distances between event means
def sliding_window_dist(sig_diffs1, sig_diffs2, num_bases):
    return np.sqrt(min(np.sum(np.square(
        sig_diffs1[i1:i1+num_bases] - sig_diffs2[i2:i2+num_bases]))
                       for i1 in range(len(sig_diffs1) - num_bases)
                       for i2 in range(len(sig_diffs2) - num_bases)))

def euclidian_dist(sig_diffs1, sig_diffs2):
    return np.sqrt(np.sum(np.square(sig_diffs1 - sig_diffs2)))

def get_pairwise_dists(reg_sig_diffs, index_q, dists_q, num_bases):
    while not index_q.empty():
        try:
            index = index_q.get(block=False)
        except Queue.Empty:
            break

        row_dists = np.array(
            [euclidian_dist(reg_sig_diffs[index], reg_sig_diffs[j])
             for j in range(0,index+1)] +
            [0 for _ in range(index+1, len(reg_sig_diffs))])
        dists_q.put((index, row_dists))

    return


if __name__ == '__main__':
    raise NotImplementedError, (
        'This is a module. See commands with `nanoraw -h`')
