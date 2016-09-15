import os, sys

import h5py

import numpy as np

from itertools import repeat, groupby
from collections import defaultdict, namedtuple

from helper import normalize_raw_signal

# option to include reads partially overlapping regions
# may need some slight edge case fixes for this to work
# i.e. getting sequnce for region from first read
PARTIAL_OVERLAP = False

DO_PROFILE = False

readData = namedtuple('readData', (
    'start', 'end', 'segs', 'read_start_rel_to_raw',
    'strand', 'means', 'fn'))

try:
    import rpy2.robjects as r
    from rpy2.robjects.packages import importr
    ggplot = importr("ggplot2")
    # TODO: Add titles to read and separeate strands by groups when plotting
    r.r('''
    plotSingleRun <- function(dat, quantDat, BaseDat, TitleDat){
    regions <- sort(c(unique(as.character(dat$Region)),
                 unique(as.character(quantDat$Region))))
    for(reg_i in regions){
    reg_base_dat <- BaseDat[BaseDat$Region==reg_i,]
    title <- TitleDat[TitleDat$Region==reg_i,'Title']
    if(reg_i %in% dat$Region){
    reg_sig_dat <- dat[dat$Region == reg_i,]
    p <- ggplot(reg_sig_dat) +
        geom_path(aes(x=Position, y=Signal, group=Read),
                  alpha=0.3, size=0.05, show.legend=FALSE)
    } else {
    reg_quant_dat <- quantDat[quantDat$Region == reg_i,]
    p <- ggplot(reg_quant_dat) +
        geom_rect(aes(xmin=Position, xmax=Position+1,
                      ymin=Lower, ymax=Upper),
                  alpha=0.1, show.legend=FALSE) +
        ylab('Signal')
    }
    print(p + facet_grid(Strand ~ .) +
        geom_text(aes(x=Position, y=-5, label=Base, color=Base),
                  data=reg_base_dat,
                  hjust=0, vjust=0, size=3, show.legend=FALSE) +
        scale_color_manual(values=c(
            'A'='#00CC00', 'C'='#0000CC', 'G'='#FFB300', 'T'='#CC0000')) +
        geom_vline(xintercept=min(reg_base_dat$Position):(
                              max(reg_base_dat$Position) + 1), size=0.01) +
        ggtitle(title) +
        theme_bw())
}}
''')
    plotSingleRun = r.globalenv['plotSingleRun']

    r.r('''
    plotGroupComp <- function(dat, quantDat, baseDat, TitleDat){
    regions <- sort(c(unique(as.character(dat$Region)),
                        as.character(unique(quantDat$Region))))
    for(reg_i in regions){
    reg_base_dat <- baseDat[baseDat$Region==reg_i,]
    title <- TitleDat[TitleDat$Region==reg_i,'Title']
    if(reg_i %in% dat$Region){
    reg_sig_dat <- dat[dat$Region == reg_i,]
    p <- ggplot(reg_sig_dat) +
        geom_path(aes(x=Position, y=Signal, color=Group, group=Read),
                  alpha=0.3, size=0.05, show.legend=FALSE)
    } else {
    reg_quant_dat <- quantDat[quantDat$Region == reg_i,]
    p <- ggplot(reg_quant_dat) +
        geom_rect(aes(xmin=Position, xmax=Position+1,
                      ymin=Lower, ymax=Upper, color=Group),
                  alpha=0.1, show.legend=FALSE) +
        ylab('Signal')
    }
    print(p + facet_grid(Strand ~ .) +
        geom_text(aes(x=Position, y=-5, label=Base, color=Base),
                  data=reg_base_dat,
                  hjust=0, vjust=0, size=3, show.legend=FALSE) +
        scale_color_manual(values=c(
            'A'='#00CC00', 'C'='#0000CC', 'G'='#FFB300', 'T'='#CC0000',
            'Group1'='black', 'Group2'='red')) +
        geom_vline(xintercept=min(reg_base_dat$Position):(
                              max(reg_base_dat$Position) + 1), size=0.01) +
        ggtitle(title) +
        theme_bw())
}}
''')
    plotGroupComp = r.globalenv['plotGroupComp']
except:
    sys.stderr.write(
        '*' * 60 + '\nERROR: Must have rpy2, R and ' +
        'R package ggplot2 installed in order to plot.\n' +
        '*' * 60 + '\n\n')
    raise

COMP_BASES = {'A':'T', 'C':'G', 'G':'C', 'T':'A', '-':'-'}
def rev_comp(seq):
    return [COMP_BASES[b] for b in seq[::-1]]

def parse_files(files, corrected_group, norm_type='median',
                outlier_threshold=5):
    raw_read_coverage = defaultdict(list)
    for read_fn in files:
        with h5py.File(read_fn) as read_data:
            if 'Analyses/' + corrected_group not in read_data:
                continue
            align_data = dict(
                read_data['Analyses/' + corrected_group +
                          '/Alignment/'].attrs.items())
            segs = np.array(read_data['Analyses/' + corrected_group +
                                      '/template/Segments'])
            base_means = np.array(read_data[
                'Analyses/' + corrected_group +
                '/template/Events']['norm_mean'])
            read_start_rel_to_raw = read_data[
                'Analyses/' + corrected_group +
                '/template/Segments'].attrs['read_start_rel_to_raw']
            raw_read_coverage[align_data['mapped_chrom']].append(
                readData(
                    align_data['mapped_start'],
                    align_data['mapped_start'] + len(segs) - 1,
                    segs, read_start_rel_to_raw,
                    align_data['mapped_strand'], base_means, read_fn))

    return raw_read_coverage

def get_base_signal(raw_read_coverage, chrm_sizes):
    # create lists for each base to contain all signal segments
    # which overlap that base
    base_signal = dict(
        ((chrm, strand), {'base_sums':np.zeros(chrm_len),
                          'base_cov':np.zeros(chrm_len, dtype=np.int_)})
        for chrm, chrm_len in chrm_sizes.items()
        for strand in ('+', '-'))

    # calculate signal on each strand separately
    for chrm in chrm_sizes.keys():
        for r_data in raw_read_coverage[chrm]:
            strand = r_data.strand
            read_means = (r_data.means if strand == '+'
                          else r_data.means[::-1])
            base_signal[(chrm, strand)]['base_sums'][
                r_data.start:r_data.start +
                len(read_means)] += read_means
            base_signal[(chrm, strand)]['base_cov'][
                r_data.start:r_data.start +
                len(read_means)] += 1

    # take the mean over all signal overlapping each base
    old_err_settings = np.seterr(all='ignore')
    mean_base_signal = {}
    for chrm_strand, chrm_sum_cov in base_signal.items():
        mean_base_signal[chrm_strand] = np.nan_to_num(
            chrm_sum_cov['base_sums'] / chrm_sum_cov['base_cov'])
    foo = np.seterr(**old_err_settings)

    return mean_base_signal

def get_quant_data(all_reg_data, plot_signal, num_bases,
                   corrected_group, group_num='Group1',
                   pcntls=[1,10,20,30,40,49]):
    upper_pcntls = [100 - pcntl for pcntl in pcntls]
    Position, Lower, Upper, Strand, Region = [], [], [], [], []
    for reg_plot_sig, (region_i, interval_start, chrom, reg_reads) in zip(
            plot_signal, all_reg_data):
        if reg_plot_sig: continue
        def get_reg_events(r_data):
            r_means = r_data.means if r_data.strand == "+" else r_data.means[::-1]
            skipped_bases = interval_start - r_data.start

            return r_means[skipped_bases:skipped_bases + num_bases + 1]

        for strand in ('+', '-'):
            if sum(r_data.strand == strand for r_data in reg_reads) == 0:
                continue
            for pos, base_read_means in enumerate(np.column_stack(
                    [get_reg_events(r_data) for r_data in reg_reads
                     if r_data.strand == strand])):
                Position.extend(list(repeat(
                    pos + interval_start, len(pcntls))))
                Lower.extend(np.percentile(
                    base_read_means, pcntls, interpolation='nearest'))
                Upper.extend(np.percentile(
                    base_read_means, upper_pcntls, interpolation='nearest'))
                Strand.extend(list(repeat(strand, len(pcntls))))
                Region.extend(list(repeat(region_i, len(pcntls))))

    return r.DataFrame({
        'Position':r.FloatVector(Position),
        'Lower':r.FloatVector(Lower),
        'Upper':r.FloatVector(Upper),
        'Strand':r.StrVector(Strand),
        'Region':r.StrVector(Region),
        'Group':r.StrVector(list(repeat(group_num, len(Position))))})

def get_signal(read_fn, read_start_rel_to_raw, num_obs):
    with h5py.File(read_fn) as read_data:
        r_sig, shift, scale = normalize_raw_signal(
            read_data['Raw/Reads'].values()[0]['Signal'],
            read_start_rel_to_raw, num_obs, 'median', None, 5)

    return r_sig

def get_signal_data(all_reg_data, plot_signal, num_bases,
                    corrected_group, group_num='Group1'):
    Position, Signal, Read, Strand, Region = [], [], [], [], []
    for reg_plot_sig, (region_i, interval_start, chrom, reg_reads) in zip(
            plot_signal, all_reg_data):
        if not reg_plot_sig: continue
        for r_num, r_data in enumerate(reg_reads):
            r_strand = r_data.strand

            segs = np.array(r_data.segs)
            r_sig = get_signal(r_data.fn, r_data.read_start_rel_to_raw,
                               r_data.segs[-1])
            if r_strand == "-":
                segs = (r_data.segs[::-1] * -1) + r_data.segs[-1]
                r_sig = r_sig[::-1]

            skipped_bases = interval_start - r_data.start
            overlap_seg_data = segs[
                skipped_bases:skipped_bases + num_bases + 1]

            for base_i, (start, stop) in enumerate(zip(
                    overlap_seg_data[:-1], overlap_seg_data[1:])):
                Position.extend(interval_start + base_i + np.linspace(
                    0, 1, stop - start, endpoint=False))
                Signal.extend(r_sig[start:stop])
                Read.extend(list(repeat(
                    str(r_num) + '_' + group_num, stop - start)))
                Strand.extend(list(repeat(r_strand, stop - start)))
                Region.extend(list(repeat(region_i, stop - start)))

    return r.DataFrame({
        'Position':r.FloatVector(Position),
        'Signal':r.FloatVector(Signal),
        'Read':r.StrVector(Read),
        'Strand':r.StrVector(Strand),
        'Region':r.StrVector(Region),
        'Group':r.StrVector(list(repeat(group_num, len(Position))))})

def get_base_data(all_reg_data, corrected_group, num_bases):
    BaseStart, Bases, BaseRegion = [], [], []
    for region_i, interval_start, chrom, reg_reads in all_reg_data:
        # get seq data from first read FAST5 file
        with h5py.File(reg_reads[0].fn) as r0_data:
            seq = ''.join(r0_data[
                'Analyses/' + corrected_group +
                '/template/Events']['base'])
        r_base_data = seq if reg_reads[0].strand == "+" else rev_comp(seq)
        reg_base_data = r_base_data[
            interval_start - reg_reads[0].start:
            interval_start - reg_reads[0].start + num_bases]
        for i, base in enumerate(reg_base_data):
            BaseStart.append(str(i + interval_start))
            Bases.append(base)
            BaseRegion.append(region_i)

    return r.DataFrame({
        'Position':r.FloatVector(BaseStart),
        'Base':r.StrVector(Bases),
        'Region':r.StrVector(BaseRegion)})

def get_region_reads(interval_data, raw_read_coverage, num_bases):
    all_reg_data = []
    for region_i, (stat, interval_start, chrom) in interval_data:
        # get all reads that overlap this interval
        if PARTIAL_OVERLAP:
            all_reg_data.append((region_i, interval_start, chrom, [
                r_data for r_data in raw_read_coverage[chrom]
                if not (r_data.start > interval_start + num_bases or
                        r_data.end < interval_start)]))
        else:
            all_reg_data.append((region_i, interval_start, chrom, [
                r_data for r_data in raw_read_coverage[chrom]
                if r_data.start <= interval_start and
                r_data.end >= interval_start + num_bases]))

    return all_reg_data

def plot_max_diff(files1, files2, num_regions, corrected_group,
                  overplot_thresh, fn_base, num_bases=100):
    if VERBOSE: sys.stderr.write('Parsing files.\n')
    raw_read_coverage1 = parse_files(files1, corrected_group)
    raw_read_coverage2 = parse_files(files2, corrected_group)

    chrm_sizes = dict((chrm, max(
        [r_data.end for r_data in raw_read_coverage1[chrm]] +
        [r_data.end for r_data in raw_read_coverage2[chrm]]))
                       for chrm in raw_read_coverage1)

    if VERBOSE: sys.stderr.write('Getting base signal.\n')
    base_signal1 = get_base_signal(raw_read_coverage1, chrm_sizes)
    base_signal2 = get_base_signal(raw_read_coverage2, chrm_sizes)

    if VERBOSE: sys.stderr.write('Get differences between base signal.\n')
    # get num_region max diff regions from each chrm then find
    # global largest after
    largest_diff_indices = []
    for chrm, chrm_size in chrm_sizes.items():
        chrm_diffs = np.concatenate([
            np.abs(base_signal1[(chrm, '+')] -
                   base_signal2[(chrm, '+')]),
            np.abs(base_signal1[(chrm, '-')] -
                   base_signal2[(chrm, '-')])])
        chrm_max_diff_regs = np.argsort(chrm_diffs)[::-1][:num_regions]
        # mod since we looked across both strands
        largest_diff_indices.extend((
            chrm_diffs[pos],
            max(np.mod(pos, chrm_size) - int(num_bases / 2.0), 0),
            chrm) for pos in chrm_max_diff_regs)

    plot_intervals = zip(range(num_regions), sorted(
        largest_diff_indices, reverse=True)[:num_regions])

    ## get reads overlapping each region
    all_reg_data1 = get_region_reads(
        plot_intervals, raw_read_coverage1, num_bases)
    all_reg_data2 = get_region_reads(
        plot_intervals, raw_read_coverage2, num_bases)
    ## show warning for low coverage regions from either group
    if any(len(reg_data) == 0 for reg_data in
           all_reg_data1 + all_reg_data2):
        if VERBOSE: sys.stderr.write(
            '*' * 60 + '\nWarning: Some regions include only reads ' +
            'from one group. This may casue some issues in plotting. ' +
            'Probably too few reads or insufficient coverage ' +
            'supplied to script.\n' + '*' * 60 + '\n')

    ## determine whether signal or quantiles
    ## (due to overplotting) should be plotted
    strand_cov = [
        (sum(r_data.strand == '+' for r_data in reg_data1[3]),
         sum(r_data.strand == '-' for r_data in reg_data1[3]),
         sum(r_data.strand == '+' for r_data in reg_data2[3]),
         sum(r_data.strand == '-' for r_data in reg_data2[3]))
        for reg_data1, reg_data2 in zip(all_reg_data1, all_reg_data2)]
    plot_signal = [max(covs) < overplot_thresh for covs in strand_cov]
    Titles = r.DataFrame({
        'Title':r.StrVector([
            chrm + " ::: Group1 Coverage (Black): " +
            str(r_cov[0]) + " + " +
            str(r_cov[1]) + " -; Group2 Coverage (Red): " +
            str(r_cov[2]) + " + " +
            str(r_cov[3]) + " -" for chrm, r_cov in zip(
                zip(*zip(*plot_intervals)[1])[2], strand_cov)]),
        'Region':r.StrVector(zip(*plot_intervals)[0])})

    if VERBOSE: sys.stderr.write('Getting plot data.\n')
    # bases are the same from either group so only get from first
    BasesData = get_base_data(
        all_reg_data1, corrected_group, num_bases)

    # get plotting data for either quantiles of raw signal
    SignalData1 = get_signal_data(
        all_reg_data1, plot_signal, num_bases, corrected_group, 'Group1')
    SignalData2 = get_signal_data(
        all_reg_data2, plot_signal, num_bases, corrected_group, 'Group2')
    QuantData1 = get_quant_data(
        all_reg_data1, plot_signal, num_bases, corrected_group, 'Group1')
    QuantData2 = get_quant_data(
        all_reg_data2, plot_signal, num_bases, corrected_group, 'Group2')

    if VERBOSE: sys.stderr.write('Plotting.\n')
    r.r('pdf("' + fn_base + '.compare_groups.pdf", height=5, width=11)')
    plotGroupComp(r.DataFrame.rbind(SignalData1, SignalData2),
                  r.DataFrame.rbind(QuantData1, QuantData2),
                  BasesData, Titles)
    r.r('dev.off()')

    return

def plot_max_coverage(files, num_regions, corrected_group,
                      overplot_thresh, fn_base, num_bases=100):
    if VERBOSE: sys.stderr.write('Parsing files.\n')
    raw_read_coverage = parse_files(files, corrected_group)

    if VERBOSE: sys.stderr.write('Calculating read coverage.\n')
    read_coverage = []
    for chrom, reads_data in raw_read_coverage.items():
        max_end = max(r_data.end for r_data in reads_data)
        chrom_coverage = np.zeros(max_end, dtype=np.int_)
        for r_data in reads_data:
            chrom_coverage[r_data.start:r_data.end] += 1

        coverage_regions = [
            (x, len(list(y))) for x, y in groupby(chrom_coverage)]
        read_coverage.extend(zip(
            zip(*coverage_regions)[0],
            np.cumsum(np.insert(zip(*coverage_regions)[1], 0, 0)),
            repeat(chrom)))

    if VERBOSE: sys.stderr.write('Getting plot data.\n')
    plot_intervals = zip(range(num_regions),
                         sorted(read_coverage, reverse=True)[:num_regions])
    all_reg_data = get_region_reads(
        plot_intervals, raw_read_coverage, num_bases)

    strand_cov = [
        (sum(r_data.strand == '+' for r_data in reg_data[3]),
         sum(r_data.strand == '-' for r_data in reg_data[3]))
        for reg_data in all_reg_data]
    plot_signal = [max(covs) < overplot_thresh for covs in strand_cov]
    Titles = r.DataFrame({
        'Title':r.StrVector([
            chrm + " ::: Coverage: " +
            str(r_cov[0]) + " + " +
            str(r_cov[1]) + " -" for chrm, r_cov in zip(
                zip(*zip(*plot_intervals)[1])[2], strand_cov)]),
        'Region':r.StrVector(zip(*plot_intervals)[0])})

    BasesData = get_base_data(
        all_reg_data, corrected_group, num_bases)
    SignalData = get_signal_data(
        all_reg_data, plot_signal, num_bases, corrected_group)
    QuantData = get_quant_data(
        all_reg_data, plot_signal, num_bases, corrected_group)

    if VERBOSE: sys.stderr.write('Plotting.\n')
    r.r('pdf("' + fn_base + '.pdf", height=5, width=11)')
    plotSingleRun(SignalData, QuantData, BasesData, Titles)
    r.r('dev.off()')

    return

def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(
        description='Plot raw signal corrected with correct_raw.' )
    parser.add_argument('fast5_basedir',
                        help='Directory containing fast5 files.')

    parser.add_argument(
        '--fast5-basedir2',
        help='Second directory containing fast5 files. '+
        'If provided regions centered base with largest '+
        'difference in mean signal will be plotted. If ' +
        'not provide regions with max coverage will be plotted')
    parser.add_argument('--num-regions', type=int, default=10,
                        help='Number of regions to plot.')

    parser.add_argument(
        '--corrected-group', default='RawGenomeCorrected_000',
        help='FAST5 group to plot created by correct_raw ' +
        'script. Default: %(default)s')
    parser.add_argument(
        '--overplot-threshold', type=int, default=50,
        help='Number of reads to trigger plotting quantiles ' +
        'instead of raw signal due to overplotting. Default: %(default)d')

    parser.add_argument(
        '--pdf-filebase', default='Nanopore_read_coverage',
        help='Base for PDF to store plots (suffix depends ' +
        'on plot type to avoid overwriting plots). ' +
        'Default: %(default)s')

    parser.add_argument('--quiet', '-q', default=False,
                        action='store_true',
                        help="Don't print status information.")
    args = parser.parse_args()

    global VERBOSE
    VERBOSE = not args.quiet

    return (args.fast5_basedir, args.fast5_basedir2, args.num_regions,
            args.corrected_group, args.overplot_threshold,
            args.pdf_filebase)

def main():
    (filebase1, filebase2, num_regions, corrected_group,
     overplot_thresh, fn_base) = parse_arguments()

    files1 = [os.path.join(filebase1, fn)
              for fn in os.listdir(filebase1)]

    if filebase2:
        files2 = [os.path.join(filebase2, fn)
                  for fn in os.listdir(filebase2)]
        if DO_PROFILE:
            import cProfile
            cProfile.runctx(
                "plot_max_diff(files1, files2, num_regions, " +
                "corrected_group, overplot_thresh, fn_base)",
                globals(), locals(), 'profile.plot_compare.prof')
            sys.exit()
        plot_max_diff(
            files1, files2, num_regions, corrected_group,
            overplot_thresh, fn_base)
    else:
        plot_max_coverage(files1, num_regions, corrected_group,
                          overplot_thresh, fn_base)

    return

if __name__ == '__main__':
    main()
