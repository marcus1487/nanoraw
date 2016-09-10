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

DO_PROFILE = True

readData = namedtuple('readData', (
    'start', 'end', 'segs', 'read_start_rel_to_raw',
    'alignment', 'means', 'fn'))

try:
    import rpy2.robjects as r
    from rpy2.robjects.packages import importr
    ggplot = importr("ggplot2")
    # TODO: Add titles to read and separeate strands by groups when plotting
    r.r('''
plotSingleRun <- function(dat, readDat){
for(reg_df in split(dat, dat$Region)){
    num_reads <- length(unique(reg_df$Read))
    plus_reads <- length(unique(reg_df[reg_df$Strand == '+','Read']))
    minus_reads <- length(unique(reg_df[reg_df$Strand == '-','Read']))
    print(ggplot(reg_df) +
        geom_path(aes(x=Position, y=Signal, group=Read),
                  alpha=0.3, size=0.05, show.legend=FALSE) +
        facet_grid(Strand ~ .) +
        geom_text(aes(x=Position, y=-5, label=Base, color=Base),
                  data=readDat[readDat$Region==reg_df$Region[1],],
                  hjust=0, vjust=0, size=3, show.legend=FALSE) +
        scale_color_manual(values=c('#00CC00', '#0000CC', '#FFB300', '#CC0000')) +
        geom_vline(xintercept=min(reg_df$Position):(max(reg_df$Position) + 1),
                   size=0.01) +
        ggtitle(paste0(num_reads, ' Reads (', plus_reads, '+, ', minus_reads, '-)')) +
        theme_bw())
}}
''')
    plotSingleRun = r.globalenv['plotSingleRun']

    r.r('''
plotGroupComp <- function(dat, readDat){
for(reg_df in split(dat, dat$Region)){
    num_reads <- length(unique(reg_df$Read))
    plus_reads <- length(unique(reg_df[reg_df$Strand == '+','Read']))
    minus_reads <- length(unique(reg_df[reg_df$Strand == '-','Read']))
    print(ggplot(reg_df) +
        geom_path(aes(x=Position, y=Signal, color=Group, group=Read),
                  alpha=0.3, size=0.05, show.legend=FALSE) +
        facet_grid(Strand ~ .) +
        geom_text(aes(x=Position, y=-5, label=Base, color=Base),
                  data=readDat[readDat$Region==reg_df$Region[1],],
                  hjust=0, vjust=0, size=3, show.legend=FALSE) +
        scale_color_manual(values=c(
            '#00CC00', '#0000CC', '#FFB300', '#CC0000', 'black', 'red')) +
        geom_vline(xintercept=min(reg_df$Position):(max(reg_df$Position) + 1),
                   size=0.01) +
        ggtitle(paste0(num_reads, ' Reads (', plus_reads, '+, ', minus_reads, '-)')) +
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
                'Analyses/' + corrected_group + '/template/Segments'].attrs[
                    'read_start_rel_to_raw']
            raw_read_coverage[align_data['mapped_chrom']].append(readData(
                align_data['mapped_start'],
                align_data['mapped_start'] + len(segs) - 1,
                segs, read_start_rel_to_raw, align_data, base_means,
                read_fn))

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
            strand = r_data.alignment['mapped_strand']
            read_means = (r_data.means if strand == '+'
                          else r_data.means[::-1])
            base_signal[(chrm, strand)]['base_sums'][
                r_data.start:r_data.start +
                len(read_means)] += read_means
            base_signal[(chrm, strand)]['base_cov'][
                r_data.start:r_data.start +
                len(read_means)] += 1

    # TODO: setting error output doesn't seem to be working
    # take the mean over all signal overlapping each base
    old_err_settings = np.seterr(divide='ignore')
    mean_base_signal = {}
    for chrm_strand, chrm_sum_cov in base_signal.items():
        mean_base_signal[chrm_strand] = np.nan_to_num(
            chrm_sum_cov['base_sums'] / chrm_sum_cov['base_cov'])
    foo = np.seterr(**old_err_settings)

    return mean_base_signal

def get_signal(read_fn, read_start_rel_to_raw, num_obs):
    with h5py.File(read_fn) as read_data:
        r_sig = normalize_raw_signal(
            read_data['Raw/Reads'].values()[0]['Signal'],
            read_start_rel_to_raw, num_obs, 'median', None, 5)

    return r_sig

def get_plot_data(interval_data, raw_read_coverage, num_bases,
                  group_num='0'):
    Position, Signal, Read, Strand, Region = [], [], [], [], []
    BaseStart, Bases, BaseRegion = [], [], []
    for region_i, (stat, interval_start, chrom) in enumerate(
            interval_data):
        # get all reads that overlap this interval
        if PARTIAL_OVERLAP:
            reg_data = [
                r_data for r_data in raw_read_coverage[chrom]
                if not (r_data.start > interval_start + num_bases or
                        r_data.end < interval_start)]
        else:
            reg_data = [
                r_data for r_data in raw_read_coverage[chrom]
                if r_data.start <= interval_start and
                r_data.end >= interval_start + num_bases]

        # TODO: When there are greater than 100 reads should plot some
        # type of quantiles instead of all reads b/c of overplotting.
        # Possibly try to alternate plots when necessary


        # TODO: This doesn't actually work. Problem with unequal factor
        # levels when plotting if this line is encountered
        if len(reg_data) == 0:
            sys.stderr.write(
                '*' * 60 + '\nWarning: No reads in region ' +
                str(region_i) + ' overlap plotted region. ' +
                'Probably too few reads or insufficient coverage ' +
                'supplied to script.\n' + '*' * 60 + '\n')
            continue
        # get seq data from first read FAST5 file
        with h5py.File(reg_data[0].fn) as r0_data:
            seq = ''.join(r0_data[
                'Analyses/RawGenomeCorrected/template/Events']['base'])
        r_base_data = seq if reg_data[0].alignment[
            'mapped_strand'] == "+" else rev_comp(seq)
        reg_base_data = r_base_data[
            interval_start - reg_data[0].start:
            interval_start - reg_data[0].start + num_bases]
        for i, base in enumerate(reg_base_data):
            BaseStart.append(str(i + interval_start))
            Bases.append(base)
            BaseRegion.append(str(region_i))

        for r_num, r_data in enumerate(reg_data):
            r_strand = r_data.alignment['mapped_strand']

            segs = np.array(r_data.segs)
            r_sig = get_signal(r_data.fn, r_data.read_start_rel_to_raw,
                               r_data.segs[-1])
            if r_strand == "-":
                segs = (r_data.segs[::-1] * -1) + r_data.segs[-1]
                r_sig = r_sig[::-1]

            skipped_bases = interval_start - r_data.alignment[
                'mapped_start']
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

    return (BaseStart, Bases, BaseRegion,
            Position, Signal, Read, Strand, Region)

def plot_max_diff(files1, files2, num_regions, corrected_group,
                  fn_base, num_bases=100):
    sys.stderr.write('Parsing files.\n')
    raw_read_coverage1 = parse_files(files1, corrected_group)
    raw_read_coverage2 = parse_files(files2, corrected_group)

    chrm_sizes = dict((chrm, max(
        [r_data.end for r_data in raw_read_coverage1[chrm]] +
        [r_data.end for r_data in raw_read_coverage2[chrm]]))
                       for chrm in raw_read_coverage1)

    sys.stderr.write('Getting base signal.\n')
    base_signal1 = get_base_signal(raw_read_coverage1, chrm_sizes)
    base_signal2 = get_base_signal(raw_read_coverage2, chrm_sizes)

    sys.stderr.write('Get differences between base signal.\n')
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

    plot_intervals = sorted(
        largest_diff_indices, reverse=True)[:num_regions]

    sys.stderr.write('Getting plot data.\n')
    (BaseStart1, Bases1, BaseRegion1,
     Position1, Signal1, Read1, Strand1, Region1) = get_plot_data(
         plot_intervals, raw_read_coverage1, num_bases, '0')
    (BaseStart2, Bases2, BaseRegion2,
     Position2, Signal2, Read2, Strand2, Region2) = get_plot_data(
         plot_intervals, raw_read_coverage2, num_bases, '1')

    sys.stderr.write('Plotting.\n')
    rawDat = r.DataFrame({
        'Position':r.FloatVector(Position1 + Position2),
        'Signal':r.FloatVector(Signal1 + Signal2),
        'Read':r.StrVector(Read1 + Read2),
        'Strand':r.StrVector(Strand1 + Strand2),
        'Region':r.StrVector(Region1 + Region2),
        'Group':r.StrVector(list(repeat('Group1', len(Position1))) +
                            list(repeat('Group2', len(Position2))))
    })
    basesDat = r.DataFrame({
        'Position':r.FloatVector(BaseStart1),
        'Base':r.StrVector(Bases1),
        'Region':r.StrVector(BaseRegion1)
    })
    r.r('pdf("' + fn_base + '.compare_groups.pdf", height=5, width=11)')
    plotGroupComp(rawDat, basesDat)
    r.r('dev.off()')

    return

def plot_max_coverage(files, num_regions, corrected_group,
                      fn_base, num_bases=100):
    sys.stderr.write('Parsing files.\n')
    raw_read_coverage = parse_files(files, corrected_group)

    sys.stderr.write('Calculating read coverage.\n')
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

    sys.stderr.write('Getting plot data.\n')
    (BaseStart, Bases, BaseRegion,
     Position, Signal, Read, Strand, Region) = get_plot_data(
         sorted(read_coverage, reverse=True)[:num_regions],
         raw_read_coverage, num_bases)

    sys.stderr.write('Plotting.\n')
    rawDat = r.DataFrame({
        'Position':r.FloatVector(Position),
        'Signal':r.FloatVector(Signal),
        'Read':r.StrVector(Read),
        'Strand':r.StrVector(Strand),
        'Region':r.StrVector(Region)
    })
    basesDat = r.DataFrame({
        'Position':r.FloatVector(BaseStart),
        'Base':r.StrVector(Bases),
        'Region':r.StrVector(BaseRegion),
    })
    r.r('pdf("' + fn_base + '.pdf", height=5, width=11)')
    plotSingleRun(rawDat, basesDat)
    r.r('dev.off()')

    return

def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(
        description='Plot raw signal corrected with correct_raw.' )
    parser.add_argument('fast5_basedir',
                        help='Directory containing fast5 files.')

    parser.add_argument('--fast5-basedir2',
                        help='Second directory containing fast5 files. '+
                        'If provided regions centered base with largest '+
                        'difference in mean signal will be plotted. If ' +
                        'not provide regions with max coverage will be plotted')
    parser.add_argument('--num-regions', type=int,
                        help='Number of regions to plot.')

    parser.add_argument('--corrected-group', default='RawGenomeCorrected_000',
                        help='FAST5 group to plot created by correct_raw ' +
                        'script. Default: %(default)s')

    parser.add_argument('--pdf-filebase', default='Nanopore_read_coverage',
                        help='Base for PDF to store plots (suffix depends ' +
                        'on plot type to avoid overwriting plots). ' +
                        'Default: %(default)s')

    parser.add_argument('--verbose', '-v', default=False,
                        action='store_true',
                        help='Whether or not to print status ' +
                        'information.')
    args = parser.parse_args()

    return args.fast5_basedir, args.fast5_basedir2, args.num_regions, \
        args.corrected_group, args.pdf_filebase

def main():
    filebase1, filebase2, num_regions, corrected_group, \
        fn_base = parse_arguments()

    files1 = [os.path.join(filebase1, fn) for fn in os.listdir(filebase1)]

    if filebase2:
        files2 = [os.path.join(filebase2, fn) for fn in os.listdir(filebase2)]
        if DO_PROFILE:
            import cProfile
            cProfile.runctx(
                "plot_max_diff(files1, files2, num_regions, corrected_group, fn_base)",
                globals(), locals(), 'profile.plot_compare.prof')
            sys.exit()
        plot_max_diff(files1, files2, num_regions, corrected_group, fn_base)
    else:
        plot_max_coverage(files1, num_regions, corrected_group, fn_base)

    return

if __name__ == '__main__':
    main()
