import os, sys

import h5py
import Queue

import numpy as np
import multiprocessing as mp

from time import sleep
from collections import defaultdict
from itertools import repeat, groupby
from pkg_resources import resource_string

# import nanoraw functions
import nanoraw_stats as ns
import nanoraw_helper as nh

VERBOSE = False

SMALLEST_PVAL = 1e-15

# quantiles and especially violin plots at leat 3 values
# to be meaningful
QUANT_MIN = 3

# plotting names for strands
FWD_STRAND = 'Forward Strand'
REV_STRAND = 'Reverse Strand'


###################################
#### ggplot via rpy2 functions ####
###################################

GG_LOAD_ERROR='*' * 60 + '\nERROR: Must have rpy2, R and ' + \
    'R package ggplot2 installed in order to plot.\n' + \
    '*' * 60 + '\n\n'
try:
    import rpy2.robjects as r
    from rpy2.robjects.packages import importr
except:
    # pass here and raise error when main functions are actually called
    pass

############################################
#### Kmer signal distribution functions ####
############################################

def plot_kmer_dist(files, corrected_group, basecall_subgroups,
                   read_mean, upstrm_bases, dnstrm_bases,
                   kmer_thresh, num_reads, pdf_fn):
    if VERBOSE: sys.stderr.write(
            'Parsing files and tabulating k-mers.\n')
    kmer_len = upstrm_bases + dnstrm_bases + 1
    reads_added = 0
    all_kmers = defaultdict(list)
    # randomly pick files instead of ordered from listing
    np.random.shuffle(files)
    for fn, basecall_subgroup in [(fn, bc_grp) for fn in files
                                  for bc_grp in basecall_subgroups]:
        try:
            read_data = h5py.File(fn, 'r')
        except IOError:
            # probably truncated file
            continue
        if ('/Analyses/' + corrected_group + '/' +
            basecall_subgroup + '/Events') not in read_data:
            continue
        event_data = read_data[
            '/Analyses/' + corrected_group + '/' + basecall_subgroup +
            '/Events'].value
        seq = event_data['base']
        means = event_data['norm_mean']
        read_kmers = defaultdict(list)
        for kmer, event_mean in zip(
                [''.join(bs) for bs in zip(*[
                    seq[i:] for i in range(kmer_len)])],
                means[kmer_len - 1 - dnstrm_bases:]):
            read_kmers[kmer].append(event_mean)
        if min(len(x) for x in read_kmers.values()) > kmer_thresh:
            reads_added += 1
            for kmer, kmer_means in read_kmers.items():
                if read_mean:
                    all_kmers[kmer].append((
                        np.mean(kmer_means), reads_added))
                else:
                    all_kmers[kmer].extend(
                        zip(kmer_means, repeat(reads_added)))

        if reads_added >= num_reads:
            break

    if VERBOSE: sys.stderr.write('Preparing plot data.\n')
    kmer_levels = [kmer for means, kmer in sorted([
        (np.mean(zip(*means)[0]), kmer)
        for kmer, means in all_kmers.items()])]

    plot_data = [
        (kmer, kmer[upstrm_bases], sig_mean, read_i)
        for kmer in kmer_levels
        for sig_mean, read_i in all_kmers[kmer]]

    kmerDat = r.DataFrame({
        'Kmer':r.FactorVector(
            r.StrVector(zip(*plot_data)[0]),
            ordered=True, levels=r.StrVector(kmer_levels)),
        'Base':r.StrVector(zip(*plot_data)[1]),
        'Signal':r.FloatVector(zip(*plot_data)[2]),
        'Read':r.StrVector(zip(*plot_data)[3])})
    # df to plot kmers as tile of colors requires cowplot R package
    try:
        cowplot = importr("cowplot")
        baseDat = r.DataFrame({
            'Kmer':r.FactorVector(
                r.StrVector([kmer for kmer in kmer_levels
                             for _ in range(kmer_len)]),
                ordered=True, levels=r.StrVector(kmer_levels)),
            'Base':r.StrVector([kmer[i] for kmer in kmer_levels
                                for i in range(kmer_len)]),
            'Position':r.IntVector([
                i - upstrm_bases for kmer in kmer_levels
                for i in range(kmer_len)])})
    except:
        sys.stderr.write(
            '********* WARNING: Install R package `cowplot` for ' +
            'visual kmer display. Using text kmer display. ********\n')
        baseDat = r.NA_Character

    if VERBOSE: sys.stderr.write('Plotting.\n')
    r.r(resource_string(__name__, 'R_scripts/plotKmerDist.R'))
    r.r('pdf("' + pdf_fn + '", height=7, width=10)')
    if read_mean:
        r.globalenv['plotKmerDistWReadPath'](kmerDat, baseDat)
    else:
        r.globalenv['plotKmerDist'](kmerDat, baseDat)
    r.r('dev.off()')

    return



########################################
#### General data parsing functions ####
########################################

def get_read_correction_data(
        filename, reg_type, reg_width, corr_basecall_group,
        region_name=None, start_at_zero=False):
    fast5_data = h5py.File(filename, 'r')
    raw_data = fast5_data['/Raw/Reads'].values()[0]
    if ( '/Analyses/' + corr_basecall_group) not in fast5_data:
        return None, None, None, None
    corr_data = fast5_data['/Analyses/' + corr_basecall_group]

    read_id = raw_data.attrs['read_id']

    # if a random region should be selected
    events_end = corr_data['Events']['start'][-1] + \
                 corr_data['Events']['length'][-1]
    if reg_type == 'start':
        reg_start = 0
    elif reg_type == 'end':
        reg_start = events_end - reg_width
    elif reg_type == 'random':
        reg_start = np.random.randint(0, events_end - reg_width)
    else:
        # reg_type should be an integer which is the raw start position
        assert isinstance(reg_type, int)
        reg_start = reg_type

    raw_offset = corr_data['Events'].attrs['read_start_rel_to_raw']
    shift = corr_data.attrs['shift']
    scale = corr_data.attrs['scale']
    lower_lim = corr_data.attrs['lower_lim']
    upper_lim = corr_data.attrs['upper_lim']
    norm_reg_signal, scale_values = nh.normalize_raw_signal(
        raw_data['Signal'].value, raw_offset + reg_start, reg_width,
        shift=shift, scale=scale, lower_lim=lower_lim,
        upper_lim=upper_lim)

    # calculate running difference
    min_seg_len = 4
    sig_cs = np.cumsum(np.insert(norm_reg_signal, 0, 0))
    running_diffs = np.abs((2 * sig_cs[min_seg_len:-min_seg_len]) -
                           sig_cs[:-2*min_seg_len] -
                           sig_cs[2*min_seg_len:])

    # note that I need to check that both new and old segments are
    # in the region as the start of the same genomic position can
    # shift in raw space (i.e. only the old or new position could be
    # in the region of interest)
    old_segs = corr_data['Alignment/read_segments'].value
    old_segs_in_reg = np.where(np.logical_and(
            reg_start <= old_segs, old_segs < reg_start + reg_width))[0]
    old_reg_segs = old_segs[old_segs_in_reg]
    old_align_vals = corr_data['Alignment/read_alignment'].value
    new_segs = np.concatenate([corr_data['Events']['start'],
                               [events_end,]])
    new_segs_in_reg = np.where(np.logical_and(
            reg_start <= new_segs, new_segs < reg_start + reg_width))[0]
    new_reg_segs = new_segs[new_segs_in_reg]
    new_align_vals = corr_data['Alignment/genome_alignment'].value

    i_old_segs = iter(old_segs)
    i_new_segs = iter(new_segs)
    align_vals = [((old_b, next(i_old_segs) if old_b != '-' else -1),
                   (new_b, next(i_new_segs) if new_b != '-' else -1))
                  for old_b, new_b in zip(
                          old_align_vals, new_align_vals)]
    reg_align_vals = [
        ((old_b, old_pos, old_pos in old_reg_segs),
         (new_b, new_pos, new_pos in new_reg_segs))
        for (old_b, old_pos), (new_b, new_pos) in align_vals
        if old_pos in old_reg_segs or new_pos in new_reg_segs]

    # summarize alignment for old and new segments
    old_is_del, old_is_mismatch, new_is_ins = [], [], []
    last_was_del = False
    for (old_b, old_pos, old_in_reg), (
            new_b, new_pos, new_in_reg) in reg_align_vals:
        if old_b == '-' and new_in_reg:
            new_is_ins.append(True)
        elif new_b == '-' and old_in_reg:
            old_is_del.append(True)
            old_is_mismatch.append(False)
            last_was_del = True
        else:
            if new_in_reg:
                new_is_ins.append(False)
            if old_in_reg:
                if last_was_del:
                    old_is_del.append(True)
                    last_was_del = False
                else:
                    old_is_del.append(False)
                old_is_mismatch.append(old_b != new_b)

    old_bases, old_reg_segs = zip(*[
        (b, pos) for b, pos, in_reg in zip(*reg_align_vals)[0]
        if in_reg]) if len(reg_align_vals) > 0 else ([], [])
    new_bases, new_reg_segs = zip(*[
        (b, pos) for b, pos, in_reg in zip(*reg_align_vals)[1]
        if in_reg]) if len(reg_align_vals) > 0 else ([], [])

    # bring positions to zero start if aligning multiple sequences
    sig_range = range(reg_start, reg_start + reg_width)
    if start_at_zero:
        old_reg_segs = [
            old_seg_pos - reg_start for old_seg_pos in old_reg_segs]
        new_reg_segs = [
            new_seg_pos - reg_start for new_seg_pos in new_reg_segs]
        sig_range = range(0, reg_width)

    old_dat = {
        'Position':r.FloatVector(old_reg_segs),
        'Base':r.StrVector(old_bases),
        'IsDel':r.BoolVector(old_is_del),
        'IsMismatch':r.BoolVector(old_is_mismatch),
        'Read':r.StrVector([read_id for _ in range(len(old_bases))])}
    new_dat = {
        'Position':r.FloatVector(new_reg_segs),
        'Base':r.StrVector(new_bases),
        'IsIns':r.BoolVector(new_is_ins),
        'Read':r.StrVector([read_id for _ in range(len(new_bases))])}
    sig_dat = {
        'Signal':r.FloatVector(norm_reg_signal),
        'Position':r.FloatVector(sig_range),
        'Read':r.StrVector([
            read_id for _ in range(len(norm_reg_signal))])}
    diff_dat = {
        'Signal':r.FloatVector(running_diffs),
        'Position':r.FloatVector(sig_range[
            min_seg_len:len(running_diffs) + min_seg_len]),
        'Read':r.StrVector([
            read_id for _ in range(len(running_diffs))])}
    # add region is applicable
    if region_name is not None:
        old_dat['Region'] = r.StrVector([
            region_name for _ in range(len(old_bases))])
        new_dat['Region'] = r.StrVector([
            region_name for _ in range(len(new_bases))])
        sig_dat['Region'] = r.StrVector([
            region_name for _ in range(len(norm_reg_signal))])
        diff_dat['Region'] = r.StrVector([
            region_name for _ in range(len(running_diffs))])

    old_dat = r.DataFrame(old_dat)
    new_dat = r.DataFrame(new_dat)
    sig_dat = r.DataFrame(sig_dat)
    diff_dat = r.DataFrame(diff_dat)

    return old_dat, new_dat, sig_dat, diff_dat

def get_reg_events(r_data, interval_start, num_bases):
    if r_data.means is None:
        with h5py.File(r_data.fn) as read_data:
            r_means = read_data[
                'Analyses/' + r_data.corr_group + '/Events']['norm_mean']
    else:
        r_means = r_data.means
    r_means = r_means if (
        r_data.strand == "+") else r_means[::-1]
    if r_data.start > interval_start:
        # handle reads that start in middle of region
        start_overlap = interval_start + num_bases - r_data.start
        # create region with nan values
        region_means = np.empty(num_bases)
        region_means[:] = np.NAN
        region_means[-start_overlap:] = r_means[:start_overlap]
    elif r_data.end < interval_start + num_bases:
        # handle reads that end inside region
        end_overlap = r_data.end - interval_start
        # create region with nan values
        region_means = np.empty(num_bases)
        region_means[:] = np.NAN
        region_means[:end_overlap] = r_means[-end_overlap:]
    else:
        skipped_bases = interval_start - r_data.start
        region_means = r_means[
            skipped_bases:skipped_bases + num_bases]

    return region_means

def get_event_data(
        all_reg_data, plot_types, num_bases, corrected_group,
        overplot_thresh, group_num='Group1'):
    Position, Signal, Strand, Region = [], [], [], []
    for reg_plot_sig, (
            region_i, interval_start, chrom, reg_reads) in zip(
                plot_types, all_reg_data):
        if reg_plot_sig != 'Violin': continue

        for strand in ('+', '-'):
            if sum(r_data.strand == strand
                   for r_data in reg_reads) == 0:
                continue
            reg_events = [
                get_reg_events(r_data, interval_start, num_bases)
                for r_data in reg_reads if r_data.strand == strand]
            for pos, base_read_means in enumerate(
                    np.column_stack(reg_events)):
                # skip bases with no coverage
                if sum(~np.isnan(base_read_means)) == 0:
                    continue
                # remove nan  regions of reads from partial overlaps
                base_read_means = base_read_means[
                    ~np.isnan(base_read_means)]
                Position.extend(repeat(
                    pos + interval_start, base_read_means.shape[0]))
                Signal.extend(base_read_means)
                Strand.extend(repeat(
                    FWD_STRAND if strand == '+' else REV_STRAND,
                    base_read_means.shape[0]))
                Region.extend(repeat(
                    region_i, base_read_means.shape[0]))

    return r.DataFrame({
        'Position':r.IntVector(Position),
        'Signal':r.FloatVector(Signal),
        'Strand':r.StrVector(Strand),
        'Region':r.StrVector(Region),
        'Group':r.StrVector(list(repeat(group_num, len(Position))))})

def get_boxplot_data(
        all_reg_data, plot_types, num_bases, corrected_group,
        overplot_thresh, group_num='Group1'):
    (Position, SigMin, Sig25, SigMed, Sig75, SigMax, Strand, Region) = (
        [], [], [], [], [], [], [], [])
    for reg_plot_sig, (
            region_i, interval_start, chrom, reg_reads) in zip(
                plot_types, all_reg_data):
        if reg_plot_sig != 'Boxplot': continue

        for strand in ('+', '-'):
            if sum(r_data.strand == strand
                   for r_data in reg_reads) == 0:
                continue
            reg_events = [
                get_reg_events(r_data, interval_start, num_bases)
                for r_data in reg_reads if r_data.strand == strand]
            for pos, base_read_means in enumerate(
                    np.column_stack(reg_events)):
                # skip regions with no coverage
                if sum(~np.isnan(base_read_means)) == 0:
                    continue
                # remove nan  regions of reads from partial overlaps
                base_read_means = base_read_means[
                    ~np.isnan(base_read_means)]
                Position.append(pos + interval_start)
                SigMin.append(np.percentile(base_read_means, 0))
                Sig25.append(np.percentile(base_read_means, 25))
                SigMed.append(np.percentile(base_read_means, 50))
                Sig75.append(np.percentile(base_read_means, 75))
                SigMax.append(np.percentile(base_read_means, 100))
                Strand.append(
                    FWD_STRAND if strand == '+' else REV_STRAND)
                Region.append(region_i)

    return r.DataFrame({
        'Position':r.IntVector(Position),
        'SigMin':r.FloatVector(SigMin),
        'Sig25':r.FloatVector(Sig25),
        'SigMed':r.FloatVector(SigMed),
        'Sig75':r.FloatVector(Sig75),
        'SigMax':r.FloatVector(SigMax),
        'Strand':r.StrVector(Strand),
        'Region':r.StrVector(Region),
        'Group':r.StrVector(list(repeat(group_num, len(Position))))})

def get_quant_data(
        all_reg_data, plot_types, num_bases, corrected_group,
        overplot_thresh, group_num='Group1', pos_offest=0,
        pcntls=[1,10,20,30,40,49]):
    upper_pcntls = [100 - pcntl for pcntl in pcntls]
    Position, Lower, Upper, Strand, Region = [], [], [], [], []
    for reg_plot_sig, (
            region_i, interval_start, chrom, reg_reads) in zip(
                plot_types, all_reg_data):
        if reg_plot_sig != 'Quantile': continue

        for strand in ('+', '-'):
            if sum(r_data.strand == strand
                   for r_data in reg_reads) == 0:
                continue
            reg_events = [
                get_reg_events(r_data, interval_start, num_bases)
                for r_data in reg_reads if r_data.strand == strand]
            for pos, base_read_means in enumerate(
                    np.column_stack(reg_events)):
                # skip regions with no coverage
                if sum(~np.isnan(base_read_means)) == 0:
                    continue
                # remove nan  regions of reads from partial overlaps
                base_read_means = base_read_means[
                    ~np.isnan(base_read_means)]
                Position.extend(list(repeat(
                    pos + interval_start + pos_offest, len(pcntls))))
                Lower.extend(np.percentile(
                    base_read_means, pcntls, interpolation='nearest'))
                Upper.extend(np.percentile(
                    base_read_means, upper_pcntls,
                    interpolation='nearest'))
                Strand.extend(
                    list(repeat(FWD_STRAND if strand == '+' else
                                REV_STRAND, len(pcntls))))
                Region.extend(list(repeat(region_i, len(pcntls))))

    return r.DataFrame({
        'Position':r.FloatVector(Position),
        'Lower':r.FloatVector(Lower),
        'Upper':r.FloatVector(Upper),
        'Strand':r.StrVector(Strand),
        'Region':r.StrVector(Region),
        'Group':r.StrVector(list(repeat(group_num, len(Position))))})

def get_signal(read_fn, read_start_rel_to_raw, num_obs, corrected_group):
    with h5py.File(read_fn) as fast5_data:
        # retrieve shift and scale computed in correction script
        corr_subgrp = fast5_data['/Analyses/' + corrected_group]
        shift = corr_subgrp.attrs['shift']
        scale = corr_subgrp.attrs['scale']
        lower_lim = corr_subgrp.attrs['lower_lim']
        upper_lim = corr_subgrp.attrs['upper_lim']
        r_sig, scale_values = nh.normalize_raw_signal(
            fast5_data['/Raw/Reads'].values()[0]['Signal'],
            read_start_rel_to_raw, num_obs, shift=shift, scale=scale,
            lower_lim=lower_lim, upper_lim=upper_lim)

    return r_sig

def get_signal_data(
        all_reg_data, plot_types, num_bases, corrected_group,
        overplot_thresh, group_num='Group1'):
    Position, Signal, Read, Strand, Region = [], [], [], [], []
    for reg_plot_sig, (
            region_i, interval_start, chrom, reg_reads) in zip(
                plot_types, all_reg_data):
        if not reg_plot_sig in ('Signal', 'Downsample'): continue
        if reg_plot_sig == 'Downsample':
            plus_reads = [r_data for r_data in reg_reads
                          if r_data.strand == '+']
            minus_reads = [r_data for r_data in reg_reads
                           if r_data.strand == '-']
            # randomly select reads to plot if too many
            if len(plus_reads) > overplot_thresh:
                np.random.shuffle(plus_reads)
                plus_reads = plus_reads[:overplot_thresh]
            if len(minus_reads) > overplot_thresh:
                np.random.shuffle(minus_reads)
                minus_reads = minus_reads[:overplot_thresh]
            reg_reads = plus_reads + minus_reads
        for r_num, r_data in enumerate(reg_reads):
            r_strand = r_data.strand

            segs = r_data.segs
            if r_strand == "-":
                segs = (segs[::-1] * -1) + segs[-1]

            if interval_start < r_data.start:
                # handle reads that start in the middle of the interval
                start_offset = r_data.start - interval_start
                overlap_seg_data = segs[:num_bases - start_offset + 1]
            else:
                start_offset = 0
                skipped_bases = interval_start - r_data.start
                overlap_seg_data = segs[
                    skipped_bases:skipped_bases + num_bases + 1]

            num_reg_obs = overlap_seg_data[-1] - overlap_seg_data[0]
            if r_strand == "+":
                reg_start_rel_raw = (r_data.read_start_rel_to_raw +
                                     overlap_seg_data[0])
                r_sig = get_signal(
                    r_data.fn, reg_start_rel_raw, num_reg_obs,
                    r_data.corr_group)
            else:
                reg_start_rel_raw = (r_data.read_start_rel_to_raw +
                                     segs[-1] - overlap_seg_data[-1])
                r_sig = get_signal(
                    r_data.fn, reg_start_rel_raw, num_reg_obs,
                    r_data.corr_group)
                r_sig = r_sig[::-1]

            for base_i, (start, stop) in enumerate(zip(
                    overlap_seg_data[:-1], overlap_seg_data[1:])):
                Position.extend(
                    interval_start + base_i + start_offset +
                    np.linspace(0, 1, stop - start, endpoint=False))
                Signal.extend(r_sig[start-overlap_seg_data[0]:
                                    stop-overlap_seg_data[0]])
                Read.extend(list(repeat(
                    str(r_num) + '_' + group_num, stop - start)))
                Strand.extend(list(repeat(
                    FWD_STRAND if r_strand == '+' else
                    REV_STRAND, stop - start)))
                Region.extend(list(repeat(region_i, stop - start)))

    return r.DataFrame({
        'Position':r.FloatVector(Position),
        'Signal':r.FloatVector(Signal),
        'Read':r.StrVector(Read),
        'Strand':r.StrVector(Strand),
        'Region':r.StrVector(Region),
        'Group':r.StrVector(list(repeat(group_num, len(Position))))})

def get_plot_types_data(plot_args, quant_offset=0):
    SignalData = get_signal_data(*plot_args)
    QuantData = get_quant_data(*plot_args, pos_offest=quant_offset)
    BoxData = get_boxplot_data(*plot_args)
    EventData = get_event_data(*plot_args)

    return SignalData, QuantData, BoxData, EventData

def get_reg_base_data(all_reg_data, corrected_group, num_bases):
    all_reg_base_data = []
    for region_i, interval_start, chrom, reg_reads in all_reg_data:
        # try to find first read to overlap whole region
        try:
            full_cov_read = next(
                read_data for read_data in reg_reads
                if read_data.start <= interval_start and
                    read_data.end >= interval_start + num_bases)
            # get seq data from first read FAST5 file
            with h5py.File(full_cov_read.fn) as r_data:
                seq = ''.join(r_data[
                    'Analyses/' + full_cov_read.corr_group +
                    '/Events']['base'])
            r_base_data = (seq if full_cov_read.strand == "+"
                           else nh.rev_comp(seq))
            reg_base_data = r_base_data[
                interval_start - full_cov_read.start:
                interval_start - full_cov_read.start + num_bases]
        except StopIteration:
            # handle case where no read overlaps whole region
            # let each read contibute its sequence and fill the rest
            # with dashes
            reg_base_data = ['-',] * num_bases
            for read_data in reg_reads:
                with h5py.File(read_data.fn) as r_data:
                    seq = ''.join(r_data[
                        'Analyses/' + read_data.corr_group +
                        '/Events']['base'])
                if read_data.strand == "-":
                    seq = nh.rev_comp(seq)
                if read_data.start > interval_start:
                    # handle reads that start in the middle of a region
                    start_overlap = (interval_start + num_bases -
                                    read_data.start)
                    reg_base_data[-start_overlap:] = seq[:start_overlap]
                else:
                    # get the number of bases from end of read that
                    # overlap the region
                    end_overlap = read_data.end - interval_start
                    reg_base_data[:end_overlap] = seq[-end_overlap:]

            reg_base_data = ''.join(reg_base_data)

        all_reg_base_data.append(reg_base_data)

    return all_reg_base_data

def get_base_r_data(all_reg_data, all_reg_base_data):
    BaseStart, Bases, BaseRegion = [], [], []
    for (region_i, interval_start, chrom, reg_reads
    ), reg_base_data in zip(
        all_reg_data, all_reg_base_data):
        for i, base in enumerate(reg_base_data):
            BaseStart.append(str(i + interval_start))
            Bases.append(base)
            BaseRegion.append(region_i)

    return r.DataFrame({
        'Position':r.FloatVector(BaseStart),
        'Base':r.StrVector(Bases),
        'Region':r.StrVector(BaseRegion)})

def get_region_reads(
        plot_intervals, raw_read_coverage, num_bases,
        filter_no_cov=True):
    all_reg_data = []
    for region_i, (chrm, interval_start, strand, stat) in plot_intervals:
        # get all reads that overlap this interval
        # note that this includes partial overlaps as these contribute
        # to coverage and other statistics so can't really restrict to
        # full coverage as previous versions of code did
        all_reg_data.append((
            region_i, interval_start, chrm, [
                r_data for r_data in raw_read_coverage[(chrm, strand)]
                if not (r_data.start >= interval_start + num_bases or
                        r_data.end < interval_start + 1)]
            if (chrm, strand) in raw_read_coverage else []))

    no_cov_regions = [
        (len(r_data) == 0, chrm + ':' + str(start))
        for reg_i, start, chrm, r_data in all_reg_data]
    if not filter_no_cov:
        return all_reg_data, no_cov_regions

    # filter out no coverage regions
    plot_intervals = [
        p_int for p_int, no_cov in zip(plot_intervals, no_cov_regions)
        if not no_cov[0]]
    all_reg_data = [
        reg_data for reg_data in all_reg_data if len(reg_data[3]) > 0]

    if any(no_cov[0] for no_cov in no_cov_regions):
        sys.stderr.write(
            '**** WARNING **** No coverage in regions: ' +
            '; '.join([reg for no_cov, reg in no_cov_regions
                       if no_cov]) + '\n')

    return all_reg_data, plot_intervals



########################################
#### Base plotting linker functions ####
########################################

def plot_corrections(
        plot_intervals, reg_width, num_reads,
        corrected_group, basecall_subgroup, pdf_fn):
    if VERBOSE: sys.stderr.write('Preparing plot data.\n')
    OldSegDat, NewSegDat, SigDat, DiffDat = [], [], [], []
    for read_fn, reg_type in plot_intervals:
        try:
            old_dat, new_dat, signal_dat, diff_dat \
                = get_read_correction_data(
                    read_fn, reg_type, reg_width, corrected_group + '/' +
                    basecall_subgroup)
        # some FAST5 files give an error:
        #     "Can't read data (Inflate() failed)"
        except IOError:
            continue
        if old_dat is None:
            # skip reads that don't have correction slots b/c they
            # couldn't be corrected
            continue
        OldSegDat.append(old_dat)
        NewSegDat.append(new_dat)
        SigDat.append(signal_dat)
        DiffDat.append(diff_dat)
        if len(OldSegDat) >= num_reads:
            break
    if VERBOSE and len(OldSegDat) < num_reads:
        sys.stderr.write(
            'WARNING: Fewer reads than requested were able to ' +
            'be processed. Likely too few reads provided or ' +
            'those provided were not corrected.\n')
    OldSegDat = r.DataFrame.rbind(*OldSegDat)
    NewSegDat = r.DataFrame.rbind(*NewSegDat)
    SigDat = r.DataFrame.rbind(*SigDat)
    DiffDat = r.DataFrame.rbind(*DiffDat)

    if VERBOSE: sys.stderr.write('Plotting.\n')
    r.r(resource_string(__name__, 'R_scripts/plotReadCorr.R'))
    r.r('pdf("' + pdf_fn + '", height=7, width=11)')
    r.globalenv['plotReadCorr'](OldSegDat, NewSegDat, SigDat, DiffDat)
    r.r('dev.off()')

    return

def plot_multi_corrections(
        files, num_reads_per_plot, num_regions, reg_width,
        corrected_group, basecall_subgroups, pdf_fn, include_orig_bcs,
        genome_locations):
    raw_read_coverage = nh.parse_fast5s(
        files, corrected_group, basecall_subgroups)
    read_coverage = nh.get_coverage(raw_read_coverage)
    coverage_regions = []
    for chrom_strand, chrom_coverage in read_coverage.items():
        chrm_coverage_regions = [
            (x, len(list(y))) for x, y in groupby(chrom_coverage)]
        chrm_reg_starts = np.cumsum(np.insert(
            zip(*chrm_coverage_regions)[1], 0, 0))
        coverage_regions.extend(zip(
            zip(*chrm_coverage_regions)[0],
            [start + (reg_len / 2) for start, reg_len in
             zip(chrm_reg_starts, zip(*chrm_coverage_regions)[1])],
            repeat(chrom_strand[0]), repeat(chrom_strand[1])))

    if genome_locations is None:
        # randomly select regions with at least num_reads_to_plot regions
        coverage_regions = [
            (chrm, reg_center, strand)
            for stat, reg_center, chrm, strand in
            coverage_regions if stat >= num_reads_per_plot]
        np.random.shuffle(coverage_regions)
        plot_intervals = zip(
            ['{:03d}'.format(rn) for rn in range(num_regions)],
            coverage_regions[:num_regions])
        if len(plot_intervals) < num_regions:
            sys.stderr.write(
                '*' * 60 + '\nWarning: Fewer regions contain minimum ' +
                'number of reads than requested.\n' + '*' * 60 + '\n')
    else:
        if VERBOSE: sys.stderr.write('Parsing genome locations.\n')
        parsed_locations = []
        for chrm_pos_strand in genome_locations:
            split_vals = chrm_pos_strand.replace('"', '').replace(
                "'", "").split(':')[:3]
            # default to plus strand if not specified
            if len(split_vals) == 2:
                parsed_locations.append((
                    split_vals[0], split_vals[1], '+'))
            else:
                parsed_locations.append(split_vals)
        plot_intervals = [
            ('{:03d}'.format(i), (chrm, int(pos) - 1, strand))
            for i, (chrm, pos, strand) in enumerate(parsed_locations)]
        # filter regions with no coverage
        plot_intervals = [
            (reg_i, (chrm, start, strand))
            for (reg_i, (chrm, start, strand)) in plot_intervals
            if read_coverage[(chrm, strand)][start] > 0]
        if len(plot_intervals) < len(parsed_locations):
            sys.stderr.write(
                '*' * 60 + '\nWarning: Some regions did not contain ' +
                'read coverage.\n' + '*' * 60 + '\n')

    if len(plot_intervals) == 0:
        sys.stderr.write(
            '*' * 60 + '\nERROR: No regions contain minimum ' +
            'number of reads.\n' + '*' * 60 + '\n')
        sys.exit()

    if VERBOSE: sys.stderr.write('Preparing plot data.\n')
    OldSegDat, NewSegDat, SigDat = [], [], []
    for reg_i, (chrm, reg_center, strand) in plot_intervals:
        reg_num_reads = 0
        ## get num_reads_per_region reads from this region
        reg_reads = [
            r_data for r_data in raw_read_coverage[(chrm, strand)]
            if r_data.start <= reg_center - (reg_width / 2.0) and
            r_data.end > reg_center + (reg_width / 2.0) and
            r_data.strand == strand]
        for r_data in reg_reads:
            # calculate raw position start
            if strand == '+':
                raw_start = int(r_data.segs[reg_center - r_data.start]
                                - (reg_width / 2) - 1)
            else:
                raw_start = int((r_data.segs[
                    len(r_data.segs) - (reg_center - r_data.start) - 1]
                                 - reg_width) + (reg_width / 2))
            try:
                old_dat, new_dat, signal_dat, diff_dat \
                    = get_read_correction_data(
                        r_data.fn, raw_start, reg_width,
                        r_data.corr_group, reg_i, True)
            # some FAST5 files give an error:
            #     "Can't read data (Inflate() failed)"
            except IOError:
                continue
            if old_dat is None:
                # skip reads that don't have correction slots b/c they
                # couldn't be corrected
                continue
            OldSegDat.append(old_dat)
            NewSegDat.append(new_dat)
            SigDat.append(signal_dat)
            reg_num_reads += 1
            if reg_num_reads >= num_reads_per_plot:
                break
        if reg_num_reads < num_reads_per_plot:
            # TODO: figure out if we should warn here
            pass
    OldSegDat = r.DataFrame.rbind(*OldSegDat)
    NewSegDat = r.DataFrame.rbind(*NewSegDat)
    SigDat = r.DataFrame.rbind(*SigDat)

    if VERBOSE: sys.stderr.write('Plotting.\n')
    r.r(resource_string(__name__, 'R_scripts/plotMultiReadCorr.R'))
    r.r('pdf("' + pdf_fn + '", height=5, width=11)')
    if include_orig_bcs:
        r.globalenv['plotMultiReadCorr'](OldSegDat, NewSegDat, SigDat)
    else:
        r.globalenv['plotMultiReadCorrNoOrig'](NewSegDat, SigDat)
    r.r('dev.off()')

    return

def plot_single_sample(
        plot_intervals, raw_read_coverage, num_bases, overplot_thresh,
        overplot_type, corrected_group, pdf_fn):
    if VERBOSE: sys.stderr.write('Preparing plot data.\n')
    all_reg_data, plot_intervals = get_region_reads(
        plot_intervals, raw_read_coverage, num_bases)
    if len(plot_intervals) == 0:
        sys.stderr.write(
            '*' * 60 + '\nERROR: No reads in any selected regions.\n'
            + '*' * 60 + '\n')
        sys.exit()

    strand_cov = [
        (sum(r_data.strand == '+' for r_data in reg_data[3]),
         sum(r_data.strand == '-' for r_data in reg_data[3]))
        for reg_data in all_reg_data]
    if overplot_type == "Downsample":
        plot_types = ["Downsample" for covs in strand_cov]
        dnspl_stars = [
            ['*' if grp_r_ovr_cov > overplot_thresh else ''
             for grp_r_ovr_cov in r_cov] for r_cov in strand_cov]
    else:
        plot_types = [
            'Signal' if (max(covs) < overplot_thresh or
                         min(covs) < QUANT_MIN)
            else overplot_type for covs in strand_cov]
        dnspl_stars = [['' for _ in r_cov] for r_cov in strand_cov]
    Titles = r.DataFrame({
        'Title':r.StrVector([
            chrm + ":" + strand + ' ' + stat +
            " ::: Coverage: " + str(r_cov[0]) + r_ovp[0] + " + " +
            str(r_cov[1]) + r_ovp[1] + " -"
            for (chrm, i_start, strand, stat), r_cov, r_ovp in zip(
                    zip(*plot_intervals)[1], strand_cov,
                    dnspl_stars)]),
        'Region':r.StrVector(zip(*plot_intervals)[0])})

    all_reg_base_data = get_reg_base_data(
        all_reg_data, corrected_group, num_bases)
    BasesData = get_base_r_data(all_reg_data, all_reg_base_data)
    SignalData, QuantData, BoxData, EventData = get_plot_types_data(
        (all_reg_data, plot_types, num_bases, corrected_group,
         overplot_thresh, 'Group1'))

    if VERBOSE: sys.stderr.write('Plotting.\n')
    r.r(resource_string(__name__, 'R_scripts/plotSingleRun.R'))
    r.r('pdf("' + pdf_fn + '", height=5, width=11)')
    r.globalenv['plotSingleRun'](SignalData, QuantData, BoxData,
                                 EventData, BasesData, Titles)
    r.r('dev.off()')

    return

def filter_group_regs(plot_intervals, grps_reg_data, grps_no_cov):
    both_no_cov = [
        all(zip(*reg_no_covs)[0])
        for reg_no_covs in zip(*grps_no_cov)]
    both_no_cov_regs = [
        reg for reg_both_no_cov, reg in zip(
            both_no_cov, zip(*grps_no_cov[0])[1]) if reg_both_no_cov]
    if any(both_no_cov) and VERBOSE:
        sys.stderr.write(
            '*' * 60 + '\nWarning: Some regions include no reads: ' +
            '\t'.join(both_no_cov_regs) + '\n' + '*' * 60 + '\n')

    # filter plot intervals and region read data
    plot_intervals = [plot_data for plot_data, no_cov in
                      zip(plot_intervals, both_no_cov) if not no_cov]
    if len(plot_intervals) == 0:
        sys.stderr.write(
            '*' * 60 + '\nERROR: No reads in any selected regions.\n'
            + '*' * 60 + '\n')
        sys.exit()

    grps_reg_data = zip(*[reg_data for reg_data, no_cov in
                          zip(zip(*grps_reg_data), both_no_cov)
                          if not no_cov])

    return plot_intervals, grps_reg_data

def plot_two_samples(
        plot_intervals, raw_read_coverage1, raw_read_coverage2,
        num_bases, overplot_thresh, overplot_type, corrected_group,
        pdf_fn, seqs_fn=None):
    if VERBOSE: sys.stderr.write('Preparing plot data.\n')
    # get reads overlapping each region
    all_reg_data1, no_cov_regs1 = get_region_reads(
        plot_intervals, raw_read_coverage1, num_bases,
        filter_no_cov=False)
    all_reg_data2, no_cov_regs2 = get_region_reads(
        plot_intervals, raw_read_coverage2, num_bases,
        filter_no_cov=False)

    # filter regions with no coverage in either read group
    plot_intervals, (all_reg_data1, all_reg_data2) = filter_group_regs(
        plot_intervals, (all_reg_data1, all_reg_data2),
        (no_cov_regs1, no_cov_regs2))

    ## determine whether signal or quantiles
    ## (due to overplotting) should be plotted
    strand_cov = [
        (sum(r_data.strand == '+' for r_data in reg_data1[3]),
         sum(r_data.strand == '-' for r_data in reg_data1[3]),
         sum(r_data.strand == '+' for r_data in reg_data2[3]),
         sum(r_data.strand == '-' for r_data in reg_data2[3]))
        for reg_data1, reg_data2 in zip(all_reg_data1, all_reg_data2)]
    # downsample can be plotted with any number of reads on either strand
    if overplot_type == "Downsample":
        plot_types = ["Downsample" for covs in strand_cov]
        dnspl_stars = [
            ['*' if grp_r_ovr_cov > overplot_thresh else ''
             for grp_r_ovr_cov in r_cov] for r_cov in strand_cov]
    else:
        plot_types = [
            'Signal' if (max(covs) < overplot_thresh or
                         min(covs) < QUANT_MIN)
            else overplot_type for covs in strand_cov]
        dnspl_stars = [['' for _ in r_cov] for r_cov in strand_cov]
    Titles = r.DataFrame({
        'Title':r.StrVector([
            chrm + ":" + strand + ' ' + stat +
            " ::: Coverage: Group1 (Blue): " +
            str(r_cov[0]) + r_dnspl[0] + " + " +
            str(r_cov[1]) + r_dnspl[1] + " -; Group2 (Red): " +
            str(r_cov[2]) + r_dnspl[2] + " + " +
            str(r_cov[3]) + r_dnspl[3] + " -"
            for (chrm, i_start, strand, stat), r_cov, r_dnspl in zip(
                    zip(*plot_intervals)[1], strand_cov,
                    dnspl_stars)]),
        'Region':r.StrVector(zip(*plot_intervals)[0])})

    # bases are the same from either group so only get from first
    merged_reg_data = [
        (reg_id, start, chrm, reg_data1 + reg_data2)
        for (reg_id, start, chrm, reg_data1),
        (_, _, _, reg_data2) in  zip(all_reg_data1, all_reg_data2)]
    all_reg_base_data = get_reg_base_data(
        merged_reg_data, corrected_group, num_bases)
    BasesData = get_base_r_data(merged_reg_data, all_reg_base_data)

    # get plotting data for either quantiles of raw signal
    SignalData1, QuantData1, BoxData1, EventData1 = get_plot_types_data(
        (all_reg_data1, plot_types, num_bases, corrected_group,
         overplot_thresh, 'Group1'), 0.1)
    SignalData2, QuantData2, BoxData2, EventData2 = get_plot_types_data(
        (all_reg_data2, plot_types, num_bases, corrected_group,
         overplot_thresh, 'Group2'), 0.5)

    if VERBOSE: sys.stderr.write('Plotting.\n')
    r.r(resource_string(__name__, 'R_scripts/plotGroupComp.R'))
    r.r('pdf("' + pdf_fn + '", height=5, width=11)')
    r.globalenv['plotGroupComp'](
        r.DataFrame.rbind(SignalData1, SignalData2),
        r.DataFrame.rbind(QuantData1, QuantData2),
        r.DataFrame.rbind(BoxData1, BoxData2),
        r.DataFrame.rbind(EventData1, EventData2),
        BasesData, Titles, 0.4)
    r.r('dev.off()')

    if seqs_fn is not None:
        if VERBOSE: sys.stderr.write('Outputting region seqeuences.\n')
        reg_seqs = zip(zip(*merged_reg_data)[0], all_reg_base_data)
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

def plot_motif_centered_with_stats(
        raw_read_coverage1, raw_read_coverage2, plot_intervals,
        pval_locs, num_bases, corrected_group, overplot_thresh,
        overplot_type, pdf_fn):
    if VERBOSE: sys.stderr.write('Preparing plot data.\n')
    all_reg_data1, no_cov_regs1 = get_region_reads(
        plot_intervals, raw_read_coverage1, num_bases,
        filter_no_cov=False)
    all_reg_data2, no_cov_regs2 = get_region_reads(
        plot_intervals, raw_read_coverage2, num_bases,
        filter_no_cov=False)

    # TODO: Get reads regardless of strand. Needs to go back
    # to region selection
    # filter for only plus strand reads
    all_reg_data1 = [
        (reg_i, i_start, chrm, [r_data for r_data in reads
                                if r_data.strand == '+'])
        for reg_i, i_start, chrm, reads in all_reg_data1]
    all_reg_data2 = [
        (reg_i, i_start, chrm, [r_data for r_data in reads
                                if r_data.strand == '+'])
        for reg_i, i_start, chrm, reads in all_reg_data2]

    merged_reg_data = [
        (reg_id, start, chrm, reg_data1 + reg_data2)
        for (reg_id, start, chrm, reg_data1),
        (_, _, _, reg_data2) in  zip(all_reg_data1, all_reg_data2)]

    plot_types = [overplot_type for _ in merged_reg_data]

    # sigDat lists
    SignalData1, QuantData1, BoxData1, EventData1 = get_plot_types_data(
        (all_reg_data1, plot_types, num_bases, corrected_group,
         overplot_thresh, 'Group1'), 0.1)
    SignalData2, QuantData2, BoxData2, EventData2 = get_plot_types_data(
        (all_reg_data2, plot_types, num_bases, corrected_group,
         overplot_thresh, 'Group2'), 0.5)
    if overplot_type in ('Signal', 'Downsample'):
        PlotData = r.DataFrame.rbind(SignalData1, SignalData2)
    elif overplot_type == 'Quantile':
        PlotData = r.DataFrame.rbind(QuantData1, QuantData2)
    elif overplot_type == 'Boxplot':
        PlotData = r.DataFrame.rbind(BoxData1, BoxData2)
    elif overplot_type == 'Violin':
        PlotData = r.DataFrame.rbind(EventData1, EventData2)

    all_reg_base_data = get_reg_base_data(
        merged_reg_data, corrected_group, num_bases)
    BasesData = get_base_r_data(merged_reg_data, all_reg_base_data)

    # stat lists
    StatsData = r.DataFrame({
        'Position':r.FloatVector(zip(*pval_locs)[0]),
        'NegLogPValue':r.FloatVector(zip(*pval_locs)[1])})

    if VERBOSE: sys.stderr.write('Plotting.\n')
    r.r(resource_string(__name__, 'R_scripts/plotMotifStats.R'))
    r.r('pdf("' + pdf_fn + '", height=3, width=5)')
    r.globalenv['plotMotifStats'](PlotData, BasesData,
                                 StatsData, overplot_type, 0.4)
    r.r('dev.off()')

    return


#################################
#### Plot processing methods ####
#################################

def plot_max_coverage(
        files, files2, num_regions, corrected_group, basecall_subgroups,
        overplot_thresh, pdf_fn, num_bases, overplot_type, obs_filter):
    if VERBOSE: sys.stderr.write('Parsing files.\n')
    raw_read_coverage = nh.parse_fast5s(
        files, corrected_group, basecall_subgroups)
    raw_read_coverage = nh.filter_reads(raw_read_coverage, obs_filter)
    read_coverage = nh.get_coverage(raw_read_coverage)

    if files2 is None:
        coverage_regions = []
        for (chrom, strand), chrom_coverage in read_coverage.items():
            chrm_coverage_regions = [
                (x, len(list(y))) for x, y in groupby(chrom_coverage)]
            coverage_regions.extend(zip(
                zip(*chrm_coverage_regions)[0],
                np.cumsum(np.insert(
                    zip(*chrm_coverage_regions)[1], 0, 0)),
                repeat(chrom), repeat(strand)))

            plot_intervals = zip(
                ['{:03d}'.format(rn) for rn in range(num_regions)],
                [(chrm, start, strand, '')
                 for stat, start, chrm, strand in
                 sorted(coverage_regions, reverse=True)[:num_regions]])

        plot_single_sample(
            plot_intervals, raw_read_coverage, num_bases,
            overplot_thresh, overplot_type, corrected_group, pdf_fn)
    else:
        raw_read_coverage2 = nh.parse_fast5s(
            files2, corrected_group, basecall_subgroups)
        raw_read_coverage2 = nh.filter_reads(
            raw_read_coverage2, obs_filter)
        read_coverage2 = nh.get_coverage(raw_read_coverage2)
        coverage_regions = []
        # only process chromosomes in both read groups
        for (chrom, strand) in set(read_coverage).intersection(
                read_coverage2):
            chrom_coverage = read_coverage[(chrom, strand)]
            chrom_coverage2 = read_coverage2[(chrom, strand)]
            if chrom_coverage.shape[0] >= chrom_coverage2.shape[0]:
                merged_chrom_cov = np.pad(
                    chrom_coverage2, (0,chrom_coverage.shape[0] -
                                      chrom_coverage2.shape[0]),
                    'constant', constant_values=0) + chrom_coverage
            else:
                merged_chrom_cov = np.pad(
                    chrom_coverage, (0, chrom_coverage2.shape[0] -
                                     chrom_coverage.shape[0]),
                    'constant', constant_values=0) + chrom_coverage2

            chrm_coverage_regions = [
                (x, len(list(y))) for x, y in groupby(merged_chrom_cov)]
            coverage_regions.extend(zip(
                zip(*chrm_coverage_regions)[0],
                np.cumsum(np.insert(
                    zip(*chrm_coverage_regions)[1], 0, 0)),
                repeat(chrom), repeat(strand)))

            plot_intervals = zip(
                ['{:03d}'.format(rn) for rn in range(num_regions)],
                [(chrm, start, strand, '')
                 for stat, start, chrm, strand in
                 sorted(coverage_regions, reverse=True)[:num_regions]])

        plot_two_samples(
            plot_intervals, raw_read_coverage, raw_read_coverage2,
            num_bases, overplot_thresh, overplot_type, corrected_group,
            pdf_fn)

    return

def plot_genome_locations(
        files, files2, corrected_group, basecall_subgroups,
        overplot_thresh, pdf_fn, num_bases, overplot_type,
        genome_locations, obs_filter):
    if VERBOSE: sys.stderr.write('Parsing genome locations.\n')
    # ignore strand for genome location plotting
    genome_locations = [
        chrm_pos.replace('"', '').replace("'", "").split(':')[:2]
        for chrm_pos in genome_locations]
    # minus one here as all python internal coords are 0-based, but
    # genome is generally 1-based
    plot_intervals = [
        ('{:03d}'.format(i), (
            chrm, max(0, int(int(pos) -
                             np.floor(num_bases / 2.0) - 1)),
            '+', '')) for i, (chrm, pos) in enumerate(
                genome_locations)]

    if VERBOSE: sys.stderr.write('Parsing files.\n')
    raw_read_coverage = nh.parse_fast5s(
        files, corrected_group, basecall_subgroups)
    raw_read_coverage = nh.filter_reads(raw_read_coverage, obs_filter)

    if files2 is not None:
        raw_read_coverage2 = nh.parse_fast5s(
            files2, corrected_group, basecall_subgroups)
        raw_read_coverage2 = nh.filter_reads(
            raw_read_coverage2, obs_filter)
        plot_two_samples(
            plot_intervals, raw_read_coverage, raw_read_coverage2,
            num_bases, overplot_thresh, overplot_type, corrected_group,
            pdf_fn)
    else:
        plot_single_sample(
            plot_intervals, raw_read_coverage, num_bases,
            overplot_thresh, overplot_type, corrected_group, pdf_fn)

    return

def plot_motif_centered(
        files, files2, num_regions, corrected_group, basecall_subgroups,
        overplot_thresh, pdf_fn, num_bases, overplot_type, motif,
        fasta_fn, deepest_coverage, obs_filter):
    if VERBOSE: sys.stderr.write(
            'Identifying genomic k-mer locations.\n')
    fasta_records = nh.parse_fasta(fasta_fn)
    motif_pat = nh.parse_motif(motif)
    motif_len = len(motif)
    def get_motif_locs(covered_chrms):
        # TODO: search over negative strand as well
        motif_locs = []
        for chrm, seq in fasta_records.iteritems():
            if chrm not in covered_chrms: continue
            for motif_loc in motif_pat.finditer(seq):
                motif_locs.append((chrm, motif_loc.start()))

        if len(motif_locs) == 0:
            sys.stderr.write(
                'Motif (' + motif_pat.pattern +
                ') not found in genome.\n')
            sys.exit()
        elif len(motif_locs) < num_regions:
            sys.stderr.write(
                'WARNING: Motif (' + motif_pat.pattern +
                ') only found ' + str(len(motif_locs)) +
                'times in genome.\n')
            num_region = len(motif_locs)
        np.random.shuffle(motif_locs)

        return motif_locs

    if VERBOSE: sys.stderr.write('Parsing files.\n')
    raw_read_coverage = nh.parse_fast5s(
        files, corrected_group, basecall_subgroups)
    raw_read_coverage = nh.filter_reads(raw_read_coverage, obs_filter) 
    if deepest_coverage:
        read_coverage = nh.get_coverage(raw_read_coverage)
    if files2 is not None:
        raw_read_coverage2 = nh.parse_fast5s(
            files2, corrected_group, basecall_subgroups)
        raw_read_coverage2 = nh.filter_reads(
            raw_read_coverage2, obs_filter)

        covered_chrms = set(zip(*raw_read_coverage)[0]).intersection(
            zip(*raw_read_coverage2)[0])
        # filter out motif_locs to chromosomes not covered
        motif_locs = get_motif_locs(covered_chrms)

        if deepest_coverage:
            read_coverage2 = nh.get_coverage(raw_read_coverage2)
        if deepest_coverage:
            if VERBOSE: sys.stderr.write(
                    'Finding deepest coverage regions.\n')
            def get_cov(chrm, pos):
                try:
                    return min(read_coverage[(chrm, '+')][pos],
                               read_coverage[(chrm, '-')][pos],
                               read_coverage2[(chrm, '+')][pos],
                               read_coverage2[(chrm, '-')][pos])
                except (IndexError, KeyError):
                    return 0
            motif_locs_cov = sorted([
                (get_cov(chrm, pos), chrm, pos)
                for chrm, pos in motif_locs], reverse=True)
            if motif_locs_cov[0][0] == 0:
                sys.stderr.write(
                    '*' * 60 + '\nk-mer not covered ' +
                    'by both groups at any positions.\n'
                    + '*' * 60 + '\n')
                sys.exit()

            plot_intervals = zip(
                ['{:03d}'.format(rn) for rn in range(num_regions)],
                ((chrm, max(pos - int(
                    (num_bases - motif_len + 1) / 2.0), 0), '+', '')
                 for cov, chrm, pos in motif_locs_cov))
        else:
            # zip over iterator of regions that have at least a
            # read overlapping so we don't have to check all reads
            plot_intervals = zip(
                ['{:03d}'.format(rn) for rn in range(num_regions)],
                ((chrm, max(pos - int(
                    (num_bases - motif_len + 1) / 2.0), 0), strand, '')
                 for chrm, pos in motif_locs
                 for strand in ('+', '-')
                 if (any(r_data.start < pos < r_data.end
                         for r_data in raw_read_coverage[
                                 (chrm, strand)]) and
                     any(r_data2.start < pos < r_data2.end
                         for r_data2 in raw_read_coverage2[
                                 (chrm, strand)]))))

        plot_two_samples(
            plot_intervals, raw_read_coverage, raw_read_coverage2,
            num_bases, overplot_thresh, overplot_type, corrected_group,
            pdf_fn)
    else:
        covered_chrms = set(zip(*raw_read_coverage)[0])
        # filter out motif_locs to chromosomes not covered
        motif_locs = get_motif_locs(covered_chrms)

        if deepest_coverage:
            if VERBOSE: sys.stderr.write(
                    'Finding deepest coverage regions.\n')
            def get_cov(chrm, pos):
                try:
                    return min(read_coverage[(chrm, '+')][pos],
                               read_coverage[(chrm, '-')][pos])
                except IndexError, KeyError:
                    return 0
            motif_locs_cov = [
                (get_cov(chrm, pos), chrm, pos)
                for chrm, pos in motif_locs]

            plot_intervals = zip(
                ['{:03d}'.format(rn) for rn in range(num_regions)],
                ((chrm, max(pos - int(
                    (num_bases - motif_len + 1) / 2.0), 0), '+', '')
                 for cov, chrm, pos in sorted(
                         motif_locs_cov, reverse=True)))
        else:
            # zip over iterator of regions that have at least a
            # read overlapping so we don't have to check all reads
            plot_intervals = zip(
                ['{:03d}'.format(rn) for rn in range(num_regions)],
                ((chrm, max(pos - int(
                    (num_bases - motif_len + 1) / 2.0), 0), '+', '')
                 for chrm, pos in motif_locs
                 if any(r_data.start < pos < r_data.end
                        for r_data in raw_read_coverage[(chrm, '+')] +
                        raw_read_coverage[(chrm, '-')])))

        plot_single_sample(
            plot_intervals, raw_read_coverage, num_bases,
            overplot_thresh, overplot_type, corrected_group, pdf_fn)

    return

def plot_max_diff(
        files1, files2, num_regions, corrected_group, basecall_subgroups,
        overplot_thresh, pdf_fn, seqs_fn, num_bases, overplot_type,
        obs_filter):
    if VERBOSE: sys.stderr.write('Parsing files.\n')
    raw_read_coverage1 = nh.parse_fast5s(
        files1, corrected_group, basecall_subgroups)
    raw_read_coverage2 = nh.parse_fast5s(
        files2, corrected_group, basecall_subgroups)
    raw_read_coverage1 = nh.filter_reads(raw_read_coverage1, obs_filter)
    raw_read_coverage2 = nh.filter_reads(raw_read_coverage2, obs_filter)

    chrm_sizes = dict(
        (chrm, max(
            [r_data.end for r_data in raw_read_coverage1[(chrm, '+')] +
             raw_read_coverage1[(chrm, '-')]] +
            [r_data.end for r_data in raw_read_coverage2[(chrm, '+')] +
             raw_read_coverage2[(chrm, '-')]]))
        for chrm in zip(*set(raw_read_coverage1).intersection(
                raw_read_coverage2))[0])

    if VERBOSE: sys.stderr.write('Getting base signal.\n')
    base_means1 = nh.get_base_means(raw_read_coverage1, chrm_sizes)
    base_means2 = nh.get_base_means(raw_read_coverage2, chrm_sizes)

    if VERBOSE: sys.stderr.write(
            'Get differences between base signal.\n')
    # get num_region max diff regions from each chrm then find
    # global largest after
    largest_diff_indices = []
    for chrm, chrm_size in chrm_sizes.items():
        for strand in ('+', '-'):
            # calculate difference and set no coverage (nan) values
            # to zero
            chrm_diffs = np.nan_to_num(
                np.abs(base_means1[(chrm, strand)] -
                       base_means2[(chrm, strand)]))
            chrm_max_diff_regs = np.argsort(
                chrm_diffs)[::-1][:num_regions]
            largest_diff_indices.extend((
                chrm_diffs[pos], max(pos - int(num_bases / 2.0), 0),
                chrm, strand) for pos in chrm_max_diff_regs)

    plot_intervals = zip(
        ['{:03d}'.format(rn) for rn in range(num_regions)],
        [(chrm, start, strand, '(Mean diff: {:.2f})'.format(stat))
         for stat, start, chrm, strand in
         sorted(largest_diff_indices, reverse=True)[:num_regions]])

    plot_two_samples(
        plot_intervals, raw_read_coverage1, raw_read_coverage2,
        num_bases, overplot_thresh, overplot_type, corrected_group,
        pdf_fn, seqs_fn)

    return

def plot_most_signif(
        files1, files2, num_regions, corrected_group, basecall_subgroups,
        overplot_thresh, pdf_fn, seqs_fn, num_bases, overplot_type,
        test_type, obs_filter, qval_thresh, min_test_vals, stats_fn,
        fishers_method_offset):
    calc_stats = stats_fn is None or not os.path.isfile(stats_fn)
    # if we are loading stats from file. Try to load the plotting
    # regions before calculating raw read coverage
    if not calc_stats:
        if VERBOSE: sys.stderr.write('Loading statistics from file.\n')
        all_stats = ns.parse_stats(stats_fn)

    if VERBOSE: sys.stderr.write('Parsing files.\n')
    raw_read_coverage1 = nh.parse_fast5s(
        files1, corrected_group, basecall_subgroups)
    raw_read_coverage2 = nh.parse_fast5s(
        files2, corrected_group, basecall_subgroups)
    raw_read_coverage1 = nh.filter_reads(raw_read_coverage1, obs_filter)
    raw_read_coverage2 = nh.filter_reads(raw_read_coverage2, obs_filter)

    if calc_stats:
        if VERBOSE: sys.stderr.write('Calculating statistics.\n')
        all_stats = ns.get_all_significance(
            raw_read_coverage1, raw_read_coverage2, test_type,
            min_test_vals, stats_fn, fishers_method_offset)

    plot_intervals = ns.get_most_signif_regions(
        all_stats, num_bases, num_regions, qval_thresh)
    plot_two_samples(
        plot_intervals, raw_read_coverage1, raw_read_coverage2,
        num_bases, overplot_thresh, overplot_type, corrected_group,
        pdf_fn, seqs_fn)

    return

def get_region_sequences(
        plot_intervals, raw_read_coverage1, raw_read_coverage2,
        num_bases, corrected_group):
    all_reg_data1, no_cov_regs1 = get_region_reads(
        plot_intervals, raw_read_coverage1, num_bases,
        filter_no_cov=False)
    all_reg_data2, no_cov_regs2 = get_region_reads(
        plot_intervals, raw_read_coverage2, num_bases,
        filter_no_cov=False)
    merged_reg_data = [
        (reg_id, start, chrm, reg_data1 + reg_data2)
        for (reg_id, start, chrm, reg_data1),
        (_, _, _, reg_data2) in  zip(all_reg_data1, all_reg_data2)]
    all_reg_base_data = get_reg_base_data(
        merged_reg_data, corrected_group, num_bases)

    return zip(zip(*merged_reg_data)[0], all_reg_base_data)

def plot_motif_centered_signif(
        files1, files2, num_regions, corrected_group, basecall_subgroups,
        overplot_thresh, overplot_type, pdf_fn, motif, context_width,
        test_type, obs_filter, min_test_vals, num_stat_values, stats_fn,
        fishers_method_offset):
    motif_pat = nh.parse_motif(motif)
    motif_len = len(motif)

    calc_stats = stats_fn is None or not os.path.isfile(stats_fn)
    if not calc_stats:
        if VERBOSE: sys.stderr.write('Loading statistics from file.\n')
        all_stats = ns.parse_stats(stats_fn)

    if VERBOSE: sys.stderr.write('Parsing files.\n')
    raw_read_coverage1 = nh.parse_fast5s(
        files1, corrected_group, basecall_subgroups)
    raw_read_coverage2 = nh.parse_fast5s(
        files2, corrected_group, basecall_subgroups)
    raw_read_coverage1 = nh.filter_reads(raw_read_coverage1, obs_filter)
    raw_read_coverage2 = nh.filter_reads(raw_read_coverage2, obs_filter)

    if calc_stats:
        if VERBOSE: sys.stderr.write('Calculating statistics.\n')
        all_stats = ns.get_all_significance(
            raw_read_coverage1, raw_read_coverage2, test_type,
            min_test_vals, stats_fn, fishers_method_offset)

    all_stats_dict = dict(
        ((chrm, strand, pos), (pval, qval))
        for pval, qval, pos, chrm, strand, cov1, cov2 in all_stats)
    covered_poss = defaultdict(set)
    for _, _, pos, chrm, strand, _, _ in all_stats:
        covered_poss[(chrm, strand)].add(pos)

    if VERBOSE: sys.stderr.write(
            'Finding signficant regions with motif.\n')
    motif_regions_data = []
    for pval, qval, pos, chrm, strand, cov1, cov2 in all_stats:
        reg_seq = get_region_sequences(
            [('0', (chrm, pos - motif_len + 1, strand, pval)),],
            raw_read_coverage1, raw_read_coverage2,
            (motif_len * 2) - 1, corrected_group)[0][1]

        reg_match = motif_pat.search(reg_seq)
        if reg_match:
            motif_regions_data.append((
                pval, qval, pos, chrm, strand, reg_match.start()))
        if len(motif_regions_data) >= num_stat_values:
            break

    # get plot intervals for all stat regions then trim to
    # num_regions after getting all p-values for plotting
    plot_intervals = zip(
        ['{:03d}'.format(rn) for rn in range(num_stat_values)],
        [(chrm, pos - motif_len + offset - context_width + 1, strand, '')
         for pval, qval, pos, chrm, strand, offset in
         motif_regions_data])
    plot_width = motif_len + (context_width * 2)
    # note check for key in stats dict as significant position
    # may lie next to region with coverage below the threshold
    pval_locs = [(pos - start, -np.log10(
        max(SMALLEST_PVAL, all_stats_dict[(chrm, strand, pos)][0]
            if (chrm, strand, pos) in all_stats_dict else 1)))
                 for chrm, start, strand, _ in zip(*plot_intervals)[1]
                 for pos in range(start, start + plot_width)]
    # TODO: Fix so that negative strand reads are plotted too.
    # requires adding "don't reverse signal" option in getting plot
    # data
    plot_intervals = [interval_data for interval_data in plot_intervals
                      if interval_data[1][2] == '+']
    # unique genomic regions filter
    uniq_p_intervals = []
    used_intervals = defaultdict(set)
    for p_int, (chrm, start, strand, reg_name) in plot_intervals:
        # could have significant region immediately next to
        # beginning/end of reads
        interval_poss = range(start, start + plot_width)
        if start not in used_intervals[(chrm, strand)] and all(
                pos in covered_poss[(chrm, strand)]
                for pos in interval_poss):
            uniq_p_intervals.append(
                (p_int, (chrm, start, strand, reg_name)))
            used_intervals[(chrm, strand)].update(interval_poss)

    uniq_p_intervals = uniq_p_intervals[:num_regions]

    plot_motif_centered_with_stats(
        raw_read_coverage1, raw_read_coverage2, uniq_p_intervals,
        pval_locs, plot_width, corrected_group, overplot_thresh,
        overplot_type, pdf_fn)

    return

def cluster_most_signif(
        files1, files2, num_regions, qval_thresh, corrected_group,
        basecall_subgroups, pdf_fn, num_bases, test_type, obs_filter,
        min_test_vals, r_struct_fn, num_processes, fasta_fn, stats_fn,
        fishers_method_offset):
    calc_stats = stats_fn is None or not os.path.isfile(stats_fn)
    if not calc_stats:
        if VERBOSE: sys.stderr.write('Loading statistics from file.\n')
        all_stats = ns.parse_stats(stats_fn)

    if VERBOSE: sys.stderr.write('Parsing files.\n')
    raw_read_coverage1 = nh.parse_fast5s(
        files1, corrected_group, basecall_subgroups)
    raw_read_coverage2 = nh.parse_fast5s(
        files2, corrected_group, basecall_subgroups)
    raw_read_coverage1 = nh.filter_reads(raw_read_coverage1, obs_filter)
    raw_read_coverage2 = nh.filter_reads(raw_read_coverage2, obs_filter)

    # calculate positions covered by at least one read in both sets
    read_coverage1 = nh.get_coverage(raw_read_coverage1)
    read_coverage2 = nh.get_coverage(raw_read_coverage2)
    covered_poss = dict(
        (chrm_strand, set(
            np.where(read_coverage1[chrm_strand] > 0)[0]).intersection(
                np.where(read_coverage2[chrm_strand] > 0)[0]))
        for chrm_strand in set(read_coverage1).intersection(
                read_coverage2))

    if calc_stats:
        if VERBOSE: sys.stderr.write('Calculating statistics.\n')
        all_stats = ns.get_all_significance(
            raw_read_coverage1, raw_read_coverage2, test_type,
            min_test_vals, stats_fn, fishers_method_offset)

    plot_intervals = ns.get_most_signif_regions(
        all_stats, num_bases, num_regions, qval_thresh)

    # unique genomic regions filter
    uniq_p_intervals = []
    used_intervals = defaultdict(set)
    for p_int, (chrm, start, strand, reg_name) in plot_intervals:
        # could have significant region immediately next to
        # beginning/end of reads
        interval_poss = range(start, start + num_bases)
        if start not in used_intervals[(chrm, strand)] and all(
                pos in covered_poss[(chrm, strand)]
                for pos in interval_poss):
            uniq_p_intervals.append(
                (p_int, (chrm, start, strand, reg_name)))
            used_intervals[(chrm, strand)].update(interval_poss)

    # get region data if outputting R data structure
    if r_struct_fn is not None:
        if VERBOSE: sys.stderr.write('Getting sequences.\n')
        if fasta_fn is None:
            # add region sequences to column names for saved dist matrix
            reg_seqs = zip(*get_region_sequences(
                uniq_p_intervals, raw_read_coverage1, raw_read_coverage2,
                num_bases, corrected_group))[1]
        else:
            fasta_records = nh.parse_fasta(fasta_fn)
            reg_seqs = [
                fasta_records[chrm][start:start+num_bases]
                for p_int, (chrm, start, strand, reg_name)
                in uniq_p_intervals]

    if VERBOSE: sys.stderr.write('Getting base signal.\n')
    chrm_sizes = dict(
        (chrm, max(
            [r_data.end for r_data in raw_read_coverage1[(chrm, '+')] +
             raw_read_coverage1[(chrm, '-')]] +
            [r_data.end for r_data in raw_read_coverage2[(chrm, '+')] +
             raw_read_coverage2[(chrm, '-')]]))
        for chrm in zip(*set(raw_read_coverage1).intersection(
                raw_read_coverage2))[0])

    base_means1 = nh.get_base_means(raw_read_coverage1, chrm_sizes)
    base_means2 = nh.get_base_means(raw_read_coverage2, chrm_sizes)

    if VERBOSE: sys.stderr.write('Getting region signal difference.\n')
    reg_sig_diffs = [
        np.nan_to_num(
            base_means1[(chrm, strand)][start:start+num_bases] -
            base_means2[(chrm, strand)][start:start+num_bases])
        for chrm, start, strand, _ in zip(*uniq_p_intervals)[1]]

    # some code to output reg signal for discovery plotting
    """sys.stdout.write('\n'.join(
        '\t'.join(('\t'.join(map(str, reg_diff)), reg_seq,
                   chrm, str(start), strand))
        for reg_seq, (chrm, start, strand, _), reg_diff in
        zip(reg_seqs, zip(*uniq_p_intervals)[1], reg_sig_diffs)) + '\n')
    sys.exit()"""

    if VERBOSE: sys.stderr.write('Getting distance between signals.\n')
    manager = mp.Manager()
    index_q = manager.Queue()
    dists_q = manager.Queue()

    for i in range(len(reg_sig_diffs)):
        index_q.put(i)

    args = (reg_sig_diffs, index_q, dists_q, num_bases)
    processes = []
    for p_id in xrange(num_processes):
        p = mp.Process(target=ns.get_pairwise_dists,
                       args=args)
        p.start()
        processes.append(p)

    reg_sig_diff_dists = []
    while any(p.is_alive() for p in processes):
        try:
            row_dists = dists_q.get(block=False)
            reg_sig_diff_dists.append(row_dists)
        except Queue.Empty:
            sleep(1)
            continue
    # empty any entries left in queue after processes have finished
    while not dists_q.empty():
        row_dists = dists_q.get(block=False)
        reg_sig_diff_dists.append(row_dists)

    reg_sig_diff_dists = zip(*sorted(reg_sig_diff_dists))[1]

    reg_sig_diff_dists = r.r.matrix(
        r.FloatVector(np.concatenate(reg_sig_diff_dists)),
        ncol=len(reg_sig_diffs), byrow=True)

    if r_struct_fn is not None:
        reg_sig_diff_dists.colnames = r.StrVector(
            ['::'.join((seq, chrm, strand, str(start))) for seq, (
                region_i, (chrm, interval_start, strand, stat)) in
             zip(reg_seqs, uniq_p_intervals)])
        r_struct_fn = r.StrVector([r_struct_fn,])
    else:
        r_struct_fn = r.NA_Character

    if VERBOSE: sys.stderr.write('Plotting (and saving data).\n')
    r.r(resource_string(__name__, 'R_scripts/plotSigMDS.R'))
    r.r('pdf("' + pdf_fn + '", height=7, width=7)')
    r.globalenv['plotSigMDS'](reg_sig_diff_dists, r_struct_fn)
    r.r('dev.off()')

    return



#########################################
#### Non-genome based plotting mains ####
#########################################

def plot_correction_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    nh.VERBOSE = VERBOSE

    try:
        ggplot = importr("ggplot2")
    except:
        sys.stderr.write(GG_LOAD_ERROR)
        sys.exit()

    files = nh.get_files_list(args.fast5_basedirs)
    plot_intervals = zip(files, repeat(args.region_type))
    plot_corrections(
        plot_intervals, args.num_obs, args.num_reads,
        args.corrected_group, args.basecall_subgroup, args.pdf_filename)

    return

def plot_multi_correction_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    nh.VERBOSE = VERBOSE

    try:
        ggplot = importr("ggplot2")
    except:
        sys.stderr.write(GG_LOAD_ERROR)
        sys.exit()

    # make sure the signal is an odd for region centering
    num_regions = args.num_regions if args.num_regions % 2 == 0 else \
                  args.num_regions + 1
    files = nh.get_files_list(args.fast5_basedirs)
    plot_multi_corrections(
        files, args.num_reads_per_plot, num_regions, args.num_obs,
        args.corrected_group, args.basecall_subgroups, args.pdf_filename,
        args.include_original_basecalls, args.genome_locations)

    return

def kmer_dist_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    nh.VERBOSE = VERBOSE

    try:
        ggplot = importr("ggplot2")
    except:
        sys.stderr.write(GG_LOAD_ERROR)
        sys.exit()

    files = nh.get_files_list(args.fast5_basedirs)
    plot_kmer_dist(
        files, args.corrected_group, args.basecall_subgroups,
        args.read_mean, args.upstream_bases, args.downstream_bases,
        args.num_kmer_threshold, args.num_reads, args.pdf_filename)

    return



##############################################
#### Signal plot parsers and main methods ####
##############################################

def max_cov_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    nh.VERBOSE = VERBOSE

    try:
        ggplot = importr("ggplot2")
    except:
        sys.stderr.write(GG_LOAD_ERROR)
        sys.exit()

    files1, files2 = nh.get_files_lists(
        args.fast5_basedirs, args.fast5_basedirs2)

    plot_max_coverage(
        files1, files2, args.num_regions, args.corrected_group,
        args.basecall_subgroups, args.overplot_threshold,
        args.pdf_filename, args.num_bases, args.overplot_type,
        nh.parse_obs_filter(args.obs_per_base_filter))

    return

def genome_loc_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    nh.VERBOSE = VERBOSE

    try:
        ggplot = importr("ggplot2")
    except:
        sys.stderr.write(GG_LOAD_ERROR)
        sys.exit()

    files1, files2 = nh.get_files_lists(
        args.fast5_basedirs, args.fast5_basedirs2)

    plot_genome_locations(
        files1, files2, args.corrected_group, args.basecall_subgroups,
        args.overplot_threshold, args.pdf_filename,
        args.num_bases, args.overplot_type, args.genome_locations,
        nh.parse_obs_filter(args.obs_per_base_filter))

    return

def motif_loc_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    nh.VERBOSE = VERBOSE

    try:
        ggplot = importr("ggplot2")
    except:
        sys.stderr.write(GG_LOAD_ERROR)
        sys.exit()

    files1, files2 = nh.get_files_lists(
        args.fast5_basedirs, args.fast5_basedirs2)

    plot_motif_centered(
        files1, files2, args.num_regions, args.corrected_group,
        args.basecall_subgroups, args.overplot_threshold,
        args.pdf_filename, args.num_bases, args.overplot_type,
        args.motif, args.genome_fasta, args.deepest_coverage,
        nh.parse_obs_filter(args.obs_per_base_filter))

    return

def max_diff_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    nh.VERBOSE = VERBOSE

    try:
        ggplot = importr("ggplot2")
    except:
        sys.stderr.write(GG_LOAD_ERROR)
        sys.exit()

    files1, files2 = nh.get_files_lists(
        args.fast5_basedirs, args.fast5_basedirs2)

    plot_max_diff(
        files1, files2, args.num_regions, args.corrected_group,
        args.basecall_subgroups, args.overplot_threshold,
        args.pdf_filename, args.sequences_filename, args.num_bases,
        args.overplot_type,
        nh.parse_obs_filter(args.obs_per_base_filter))

    return

def signif_diff_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    nh.VERBOSE = VERBOSE
    ns.VERBOSE = VERBOSE

    try:
        ggplot = importr("ggplot2")
    except:
        sys.stderr.write(GG_LOAD_ERROR)
        sys.exit()

    files1, files2 = nh.get_files_lists(
        args.fast5_basedirs, args.fast5_basedirs2)

    plot_most_signif(
        files1, files2, args.num_regions, args.corrected_group,
        args.basecall_subgroups, args.overplot_threshold,
        args.pdf_filename, args.sequences_filename, args.num_bases,
        args.overplot_type, args.test_type,
        nh.parse_obs_filter(args.obs_per_base_filter),
        args.q_value_threshold, args.minimum_test_reads,
        args.statistics_filename, args.fishers_method_offset)

    return

def motif_signif_diff_main(args):
    try:
        ggplot = importr("ggplot2")
    except:
        sys.stderr.write(GG_LOAD_ERROR)
        sys.exit()
    try:
        cowplot = importr("cowplot")
    except:
        sys.stderr.write(
            '*' * 60 + '\nERROR: Must have R packge `cowplot` ' +
            'installed in order to create motif centered plots ' +
            '(install via `install.packages(cowplot)` from ' +
            'an R prompt).\n' + '*' * 60 + '\n\n')
        sys.exit()

    global VERBOSE
    VERBOSE = not args.quiet
    nh.VERBOSE = VERBOSE
    ns.VERBOSE = VERBOSE

    files1, files2 = nh.get_files_lists(
        args.fast5_basedirs, args.fast5_basedirs2)

    plot_motif_centered_signif(
        files1, files2, args.num_regions, args.corrected_group,
        args.basecall_subgroups, args.overplot_threshold,
        args.overplot_type, args.pdf_filename, args.motif,
        args.num_context, args.test_type,
        nh.parse_obs_filter(args.obs_per_base_filter),
        args.minimum_test_reads, args.num_statistics,
        args.statistics_filename, args.fishers_method_offset)

    return

def cluster_signif_diff_main(args):
    global VERBOSE
    VERBOSE = not args.quiet
    nh.VERBOSE = VERBOSE
    ns.VERBOSE = VERBOSE

    try:
        ggplot = importr("ggplot2")
    except:
        sys.stderr.write(GG_LOAD_ERROR)
        sys.exit()

    files1, files2 = nh.get_files_lists(
        args.fast5_basedirs, args.fast5_basedirs2)

    cluster_most_signif(
        files1, files2, args.num_regions, args.q_value_threshold,
        args.corrected_group, args.basecall_subgroups,
        args.pdf_filename, args.num_bases, args.test_type,
        nh.parse_obs_filter(args.obs_per_base_filter),
        args.minimum_test_reads, args.r_data_filename, args.processes,
        args.genome_fasta, args.statistics_filename,
        args.fishers_method_offset)

    return

if __name__ == '__main__':
    raise NotImplementedError, (
        'This is a module. See commands with `nanoraw -h`')
