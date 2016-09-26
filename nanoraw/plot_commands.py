import os, sys

import re
import h5py

import numpy as np

from copy import copy
from collections import defaultdict
from itertools import repeat, groupby

from nanoraw_helper import (
    normalize_raw_signal, parse_fast5s, parse_fasta)

DO_PROFILE = False
VERBOSE = False

# quantiles and especially violin plots at leat 3 values
# to be meaningful
QUANT_MIN = 3


COMP_BASES = {'A':'T', 'C':'G', 'G':'C', 'T':'A', '-':'-'}
def rev_comp(seq):
    return [COMP_BASES[b] for b in seq[::-1]]



###################################
#### ggplot via rpy2 functions ####
###################################

try:
    import rpy2.robjects as r
    from rpy2.robjects.packages import importr
    ggplot = importr("ggplot2")
    r.r('''
    plotSingleRun <- function(
        dat, quantDat, boxDat, eventDat, BaseDat, TitleDat){
    regions <- sort(c(unique(as.character(dat$Region)),
                 unique(as.character(quantDat$Region)),
                 unique(as.character(boxDat$Region)),
                 unique(as.character(eventDat$Region))))
    for(reg_i in regions){
    reg_base_dat <- BaseDat[BaseDat$Region==reg_i,]
    title <- TitleDat[TitleDat$Region==reg_i,'Title']
    if(reg_i %in% dat$Region){
    reg_sig_dat <- dat[dat$Region == reg_i,]
    p <- ggplot(reg_sig_dat) +
        geom_path(aes(x=Position, y=Signal, group=Read),
                  alpha=0.3, size=0.05, show.legend=FALSE)
    } else if(reg_i %in% quantDat$Region) {
    reg_quant_dat <- quantDat[quantDat$Region == reg_i,]
    p <- ggplot(reg_quant_dat) +
        geom_rect(aes(xmin=Position, xmax=Position+1,
                      ymin=Lower, ymax=Upper),
                  alpha=0.1, show.legend=FALSE) +
        ylab('Signal')
    } else if(reg_i %in% boxDat$Region) {
    reg_box_dat <- boxDat[boxDat$Region == reg_i,]
    p <- ggplot(reg_box_dat) +
        geom_boxplot(aes(Position + 0.5, ymin=SigMin, lower=Sig25,
                         middle=SigMed, upper=Sig75, ymax=SigMax),
                     size=0.2, alpha=0.5,
                     stat="identity", show.legend=FALSE) +
        ylab('Signal') + xlab('Position')
    } else {
    reg_event_dat <- eventDat[eventDat$Region == reg_i,]
    p <- ggplot(reg_event_dat) +
        geom_violin(aes(x=Position + 0.5, y=Signal, group=Position),
                    size=0, show.legend=FALSE) +
        ylab('Signal') + xlab('Position')
    }
    print(p + facet_grid(Strand ~ .) +
        geom_text(aes(x=Position+0.5, y=-5, label=Base, color=Base),
                  data=reg_base_dat,
                  hjust=0.5, vjust=0, size=3, show.legend=FALSE) +
        scale_color_manual(values=c(
            'A'='#00CC00', 'C'='#0000CC', 'G'='#FFB300', 'T'='#CC0000',
            '-'='black')) +
        geom_vline(xintercept=min(reg_base_dat$Position):(
                              max(reg_base_dat$Position) + 1),
                   size=0.01) +
        ggtitle(title) +
        theme_bw() + theme(axis.text.x=element_text(hjust=0)))
}}
''')
    plotSingleRun = r.globalenv['plotSingleRun']

    r.r('''
    plotGroupComp <- function(dat, quantDat, boxDat, eventDat,
                              baseDat, TitleDat, QuantWidth){
    regions <- sort(c(unique(as.character(dat$Region)),
                      unique(as.character(quantDat$Region)),
                      unique(as.character(boxDat$Region)),
                      unique(as.character(eventDat$Region))))
    for(reg_i in regions){
    reg_base_dat <- baseDat[baseDat$Region==reg_i,]
    title <- TitleDat[TitleDat$Region==reg_i,'Title']
    if(reg_i %in% dat$Region){
    reg_sig_dat <- dat[dat$Region == reg_i,]
    p <- ggplot(reg_sig_dat) +
        geom_path(aes(x=Position, y=Signal, color=Group, group=Read),
                  alpha=0.3, size=0.05, show.legend=FALSE)
    } else if(reg_i %in% quantDat$Region) {
    reg_quant_dat <- quantDat[quantDat$Region == reg_i,]
    p <- ggplot(reg_quant_dat) +
        geom_rect(aes(xmin=Position, xmax=Position + QuantWidth,
                      ymin=Lower, ymax=Upper, fill=Group),
                  alpha=0.1, show.legend=FALSE) +
        ylab('Signal')
    } else if (reg_i %in% boxDat$Region) {
    reg_box_dat <- boxDat[boxDat$Region == reg_i,]
    p <- ggplot(reg_box_dat) +
        geom_boxplot(aes(Position + 0.5, ymin=SigMin, lower=Sig25,
                         middle=SigMed, upper=Sig75, ymax=SigMax,
                         fill=Group), size=0.2, alpha=0.3,
                     stat="identity", show.legend=FALSE) +
        ylab('Signal') + xlab('Position')
    } else {
    reg_event_dat <- eventDat[eventDat$Region == reg_i,]
    p <- ggplot(reg_event_dat) +
        geom_violin(aes(x=Position + 0.5, y=Signal, fill=Group,
                       group=paste0(Group, Position)),
                     size=0, show.legend=FALSE) +
        ylab('Signal') + xlab('Position')
    }
    print(p + facet_grid(Strand ~ .) +
        geom_text(aes(x=Position+0.5, y=-5, label=Base, color=Base),
                  data=reg_base_dat,
                  hjust=0.5, vjust=0, size=3, show.legend=FALSE) +
        scale_color_manual(values=c(
            'A'='#00CC00', 'C'='#0000CC', 'G'='#FFB300', 'T'='#CC0000',
            '-'='black', 'Group1'='blue', 'Group2'='red')) +
        scale_fill_manual(values=c(
            'Group1'='blue', 'Group2'='red')) +
        geom_vline(xintercept=min(reg_base_dat$Position):(
                              max(reg_base_dat$Position) + 1),
                   size=0.01) +
        ggtitle(title) +
        theme_bw() + theme(axis.text.x=element_text(hjust=0)))
}}
''')
    plotGroupComp = r.globalenv['plotGroupComp']

    r.r('''
    plotKmerDist <- function(dat){
    print(ggplot(dat) +
        geom_boxplot(aes(x=Trimer, y=Signal, color=Base)) +
        theme_bw() + theme(
            axis.text.x=element_text(angle=60, hjust=1, size=8)) +
        scale_color_manual(
            values=c('#00CC00', '#0000CC', '#FFB300', '#CC0000')))
}
''')
    plotKmerDist = r.globalenv['plotKmerDist']
    r.r('''
    plotKmerDistWReadPath <- function(dat){
    print(ggplot(dat) +
        geom_boxplot(aes(x=Trimer, y=Signal, color=Base)) +
        theme_bw() +
        theme(axis.text.x=element_text(angle=60, hjust=1, size=8)) +
        scale_color_manual(
            values=c('#00CC00', '#0000CC', '#FFB300', '#CC0000')))
    print(ggplot(dat) +
        geom_path(aes(x=Trimer, y=Signal, group=Read), alpha=0.1) +
        theme_bw() +
        theme(axis.text.x=element_text(angle=60, hjust=1, size=8)) +
        scale_color_manual(
        values=c('#00CC00', '#0000CC', '#FFB300', '#CC0000')))
}
''')
    plotKmerDistWReadPath = r.globalenv['plotKmerDistWReadPath']
except:
    sys.stderr.write(
        '*' * 60 + '\nERROR: Must have rpy2, R and ' +
        'R package ggplot2 installed in order to plot.\n' +
        '*' * 60 + '\n\n')
    raise



############################################
#### Kmer signal distribution functions ####
############################################

def plot_kmer_dist(files, corrected_group, read_mean, kmer_len,
                   kmer_thresh, pdf_fn):
    if VERBOSE: sys.stderr.write('Parsing files.\n')
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

    if VERBOSE: sys.stderr.write('Tabulating k-mers.\n')
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
                    all_trimers[trimer].append((
                        np.mean(trimer_means), read_i))
                else:
                    all_trimers[trimer].extend(
                        zip(trimer_means, repeat(read_i)))

    if VERBOSE: sys.stderr.write('Preparing plot data.\n')
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
            (kmer_i, pos_i, base)
            for kmer_i, kmer in enumerate(kmer_leves)
            for pos_i, base in enumerate(kmer)]
        kmerDat = r.DataFrame({
            'Kmer':r.IntVector(zip(*kmer_plot_data)[0]),
            'Position':r.IntVector(zip(*kmer_plot_data)[1]),
            'Base':r.StrVector(zip(*kmer_plot_data)[2])})

    if VERBOSE: sys.stderr.write('Plotting.\n')
    if read_mean:
        r.r('pdf("' + pdf_fn + '", height=7, width=10)')
        plotKmerDistWReadPath(trimerDat)
        r.r('dev.off()')
    else:
        r.r('pdf("' + pdf_fn + '", height=7, width=10)')
        plotKmerDist(trimerDat)
        r.r('dev.off()')

    return

def kmer_main(args):
    global VERBOSE
    VERBOSE = not args.quiet

    files = [os.path.join(args.fast5_basedir, fn)
             for fn in os.listdir(args.fast5_basedir)]
    plot_kmer_dist(
        files, args.corrected_group, args.read_mean, args.kmer_length,
        args.num_trimer_threshold, args.pdf_filename)

    return

def get_kmer_parser():
    import argparse
    parser = argparse.ArgumentParser(
        description='Plot distribution of signal across kmers.',
        add_help=False)
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
        '--pdf-filename', default='Nanopore_kmer_distribution.pdf',
        help='PDF filename to store plot(s). Default: %(default)s')

    parser.add_argument(
        '--quiet', '-q', default=False, action='store_true',
        help="Don't print status information.")

    return parser



########################################
#### General data parsing functions ####
########################################

def get_reg_events(r_data, interval_start, num_bases):
    if r_data.means is None:
        with h5py.File(r_data.fn) as read_data:
            r_means = read_data[
                'Analyses/' + r_data.corr_group +
                '/template/Events']['norm_mean']
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
        group_num='Group1'):
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
                    strand, base_read_means.shape[0]))
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
        group_num='Group1'):
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
                Strand.append(strand)
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
        group_num='Group1', pos_offest=0,
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

def get_signal_data(all_reg_data, plot_types, num_bases,
                    corrected_group, group_num='Group1'):
    Position, Signal, Read, Strand, Region = [], [], [], [], []
    for reg_plot_sig, (
            region_i, interval_start, chrom, reg_reads) in zip(
                plot_types, all_reg_data):
        if reg_plot_sig != 'Signal': continue
        for r_num, r_data in enumerate(reg_reads):
            r_strand = r_data.strand

            segs = r_data.segs
            r_sig = get_signal(
                r_data.fn, r_data.read_start_rel_to_raw, segs[-1])
            if r_strand == "-":
                segs = (segs[::-1] * -1) + segs[-1]
                r_sig = r_sig[::-1]

            if interval_start < r_data.start:
                # handle reads that start in the middle of the interval
                start_offset = r_data.start - interval_start
                overlap_seg_data = segs[:num_bases - start_offset + 1]
            else:
                start_offset = 0
                skipped_bases = interval_start - r_data.start
                overlap_seg_data = segs[
                    skipped_bases:skipped_bases + num_bases + 1]

            for base_i, (start, stop) in enumerate(zip(
                    overlap_seg_data[:-1], overlap_seg_data[1:])):
                Position.extend(
                    interval_start + base_i + start_offset +
                    np.linspace(0, 1, stop - start, endpoint=False))
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

def get_plot_types_data(plot_args, quant_offset=0):
    SignalData = get_signal_data(*plot_args)
    QuantData = get_quant_data(*plot_args, pos_offest=quant_offset)
    BoxData = get_boxplot_data(*plot_args)
    EventData = get_event_data(*plot_args)

    return SignalData, QuantData, BoxData, EventData

def get_base_data(all_reg_data, corrected_group, num_bases):
    BaseStart, Bases, BaseRegion = [], [], []
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
                    'Analyses/' + corrected_group +
                    '/template/Events']['base'])
            r_base_data = (seq if full_cov_read.strand == "+"
                           else rev_comp(seq))
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
                        'Analyses/' + corrected_group +
                        '/template/Events']['base'])
                if read_data.strand == "-":
                    seq = rev_comp(seq)
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
    for region_i, (chrm, interval_start, strand) in plot_intervals:
        # get all reads that overlap this interval
        # note that this includes partial overlaps as these contribute
        # to coverage and other statistics so can't really restrict to
        # full coverage as previous versions of code did
        all_reg_data.append((region_i, interval_start, chrm, [
            r_data for r_data in raw_read_coverage[chrm]
            if not (r_data.start > interval_start + num_bases or
                    r_data.end < interval_start)]))

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
                r_data.start:r_data.start + len(read_means)] += 1

    # take the mean over all signal overlapping each base
    old_err_settings = np.seterr(all='ignore')
    mean_base_signal = {}
    for chrm_strand, chrm_sum_cov in base_signal.items():
        mean_base_signal[chrm_strand] = np.nan_to_num(
            chrm_sum_cov['base_sums'] / chrm_sum_cov['base_cov'])
    foo = np.seterr(**old_err_settings)

    return mean_base_signal



########################################
#### Base plotting linker functions ####
########################################

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
    plot_types = [
        'Signal' if (max(covs) < overplot_thresh or
                     min(covs) < QUANT_MIN)
        else overplot_type for covs in strand_cov]
    Titles = r.DataFrame({
        'Title':r.StrVector([
            chrm + " ::: Coverage: " +
            str(r_cov[0]) + " + " +
            str(r_cov[1]) + " -" for chrm, r_cov in zip(
                zip(*zip(*plot_intervals)[1])[0], strand_cov)]),
        'Region':r.StrVector(zip(*plot_intervals)[0])})

    BasesData = get_base_data(
        all_reg_data, corrected_group, num_bases)
    SignalData, QuantData, BoxData, EventData = get_plot_types_data(
        (all_reg_data, plot_types, num_bases, corrected_group, 'Group1'))

    if VERBOSE: sys.stderr.write('Plotting.\n')
    r.r('pdf("' + pdf_fn + '", height=5, width=11)')
    plotSingleRun(SignalData, QuantData, BoxData, EventData,
                  BasesData, Titles)
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
        pdf_fn):
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
    plot_types = [
        'Signal' if (max(covs) < overplot_thresh or
                     min(covs) < QUANT_MIN)
        else overplot_type for covs in strand_cov]
    Titles = r.DataFrame({
        'Title':r.StrVector([
            chrm + " ::: Group1 Coverage (Blue): " +
            str(r_cov[0]) + " + " +
            str(r_cov[1]) + " -; Group2 Coverage (Red): " +
            str(r_cov[2]) + " + " +
            str(r_cov[3]) + " -" for chrm, r_cov in zip(
                zip(*zip(*plot_intervals)[1])[0], strand_cov)]),
        'Region':r.StrVector(zip(*plot_intervals)[0])})

    if VERBOSE: sys.stderr.write('Getting plot data.\n')
    # bases are the same from either group so only get from first
    merged_reg_data = [
        (reg_id, start, chrm, reg_data1 + reg_data2)
        for (reg_id, start, chrm, reg_data1),
        (_, _, _, reg_data2) in  zip(all_reg_data1, all_reg_data2)]
    BasesData = get_base_data(
        merged_reg_data, corrected_group, num_bases)

    # get plotting data for either quantiles of raw signal
    SignalData1, QuantData1, BoxData1, EventData1 = get_plot_types_data(
        (all_reg_data1, plot_types, num_bases, corrected_group,
         'Group1'), 0.1)
    SignalData2, QuantData2, BoxData2, EventData2 = get_plot_types_data(
        (all_reg_data2, plot_types, num_bases, corrected_group,
         'Group2'), 0.5)

    if VERBOSE: sys.stderr.write('Plotting.\n')
    r.r('pdf("' + pdf_fn + '", height=5, width=11)')
    plotGroupComp(r.DataFrame.rbind(SignalData1, SignalData2),
                  r.DataFrame.rbind(QuantData1, QuantData2),
                  r.DataFrame.rbind(BoxData1, BoxData2),
                  r.DataFrame.rbind(EventData1, EventData2),
                  BasesData, Titles, 0.4)
    r.r('dev.off()')

    return



############################################
#### One or two sample plotting methods ####
###########################################

def plot_genome_locations(
        files, files2, corrected_group, overplot_thresh, pdf_fn,
        num_bases, overplot_type, genome_locations):
    genome_locations = [chrm_pos.split(':')
                        for chrm_pos in genome_locations]
    plot_intervals = [
        ('{:03d}'.format(i), (
            chrm, max(0, int(int(pos) - np.floor(num_bases / 2.0))),
            None)) for i, (chrm, pos) in enumerate(genome_locations)]

    raw_read_coverage = parse_fast5s(files, corrected_group)

    if files2 is not None:
        raw_read_coverage2 = parse_fast5s(files2, corrected_group)
        plot_two_samples(
            plot_intervals, raw_read_coverage, raw_read_coverage2,
            num_bases, overplot_thresh, overplot_type, corrected_group,
            pdf_fn)
    else:
        plot_single_sample(
            plot_intervals, raw_read_coverage, num_bases,
            overplot_thresh, overplot_type, corrected_group, pdf_fn)
    return

def plot_kmer_centered(
        files, files2, num_regions, corrected_group, overplot_thresh,
        pdf_fn, num_bases, overplot_type, kmer, fasta_fn,
        deepest_coverage):
    with open(fasta_fn) as fasta_fp:
        fasta_records = parse_fasta(fasta_fp)

    kmer_locs = []
    for chrm, seq in fasta_records.items():
        for kmer_loc in re.finditer(kmer, seq):
            kmer_locs.append((chrm, kmer_loc.start()))

    if len(kmer_locs) == 0:
        sys.stderr.write('Kmer (' + kmer + ') not found in genome.\n')
        sys.exit()
    elif len(kmer_locs) < num_regions:
        sys.stderr.write(
            'WARNING: Kmer (' + kmer + ') only found ' +
            str(len(kmer_locs)) + 'times in genome.\n')
        num_region = len(kmer_locs)
    np.random.shuffle(kmer_locs)

    raw_read_coverage = parse_fast5s(files, corrected_group)
    if files2 is not None:
        raw_read_coverage2 = parse_fast5s(files2, corrected_group)
        if deepest_coverage:
            raise NotImplementedError, 'Not currenly a working option.'
        else:
            # zip over iterator of regions that have at least a
            # read overlapping so we don't have to check all reads
            plot_intervals = zip(
                ['{:03d}'.format(rn) for rn in range(num_regions)],
                ((chrm, max(pos - int(
                    (num_bases - len(kmer) + 1) / 2.0), 0), None)
                 for chrm, pos in kmer_locs
                 if (any(r_data.start < pos < r_data.end
                         for r_data in raw_read_coverage[chrm]) and
                     any(r_data2.start < pos < r_data2.end
                         for r_data2 in raw_read_coverage2[chrm]))))

        plot_two_samples(
            plot_intervals, raw_read_coverage, raw_read_coverage2,
            num_bases, overplot_thresh, overplot_type, corrected_group,
            pdf_fn)
    else:
        if deepest_coverage:
            raise NotImplementedError, 'Not currenly a working option.'
        else:
            # zip over iterator of regions that have at least a
            # read overlapping so we don't have to check all reads
            plot_intervals = zip(
                ['{:03d}'.format(rn) for rn in range(num_regions)],
                ((chrm, max(pos - int(
                    (num_bases - len(kmer) + 1) / 2.0), 0), None)
                 for chrm, pos in kmer_locs
                 if any(r_data.start < pos < r_data.end
                        for r_data in raw_read_coverage[chrm])))

        plot_single_sample(
            plot_intervals, raw_read_coverage, num_bases,
            overplot_thresh, overplot_type, corrected_group, pdf_fn)

    return



########################################
#### Single sample plotting methods ####
########################################

def plot_max_coverage(files, num_regions, corrected_group,
                      overplot_thresh, pdf_fn, num_bases,
                      overplot_type):
    if VERBOSE: sys.stderr.write('Parsing files.\n')
    raw_read_coverage = parse_fast5s(files, corrected_group)

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
            repeat(chrom), repeat(None)))

    plot_intervals = zip(
        ['{:03d}'.format(rn) for rn in range(num_regions)],
        [(chrm, start, strand) for stat, start, chrm, strand in
         sorted(read_coverage, reverse=True)[:num_regions]])

    plot_single_sample(
        plot_intervals, raw_read_coverage, num_bases,
        overplot_thresh, overplot_type, corrected_group, pdf_fn)

    return

def single_sample_main(args):
    global VERBOSE
    VERBOSE = not args.quiet

    files = [os.path.join(args.fast5_basedir, fn)
             for fn in os.listdir(args.fast5_basedir)]
    if DO_PROFILE:
        import cProfile
        cProfile.runctx(
            "plot_max_coverage(" +
            "files, args.num_regions, args.corrected_group," +
            "args.overplot_threshold, args.pdf_filename, " +
            "args.num_bases, args.overplot_type)",
            globals(), locals(), 'profile.plot_compare.prof')
        sys.exit()
    if args.genome_locations is not None:
        plot_genome_locations(
            files, None, args.corrected_group, args.overplot_threshold,
            args.pdf_filename, args.num_bases, args.overplot_type,
            args.genome_locations)
    elif args.kmer:
        plot_kmer_centered(
            files, None, args.num_regions, args.corrected_group,
            args.overplot_threshold, args.pdf_filename,
            args.num_bases, args.overplot_type, args.kmer,
            args.fasta_filename, args.deepest_coverage)
    else:
        plot_max_coverage(
            files, args.num_regions, args.corrected_group,
            args.overplot_threshold, args.pdf_filename, args.num_bases,
            args.overplot_type)

    return

def get_single_sample_parser():
    import argparse
    parser = argparse.ArgumentParser(
        description='Plot raw signal from from two samples where ' +
        'FAST5 files were corrected with `nanoraw correct`.',
        add_help=False)
    parser.add_argument(
        'fast5_basedir',
        help='Directory containing fast5 files.')

    parser.add_argument(
        '--pdf-filename',
        default='Nanopore_read_coverage.pdf',
        help='PDF filename to store plots. Default: %(default)s')

    return parser



#####################################
#### Two sample plotting methods ####
#####################################

def plot_max_diff(files1, files2, num_regions, corrected_group,
                  overplot_thresh, pdf_fn, num_bases, overplot_type):
    if VERBOSE: sys.stderr.write('Parsing files.\n')
    raw_read_coverage1 = parse_fast5s(files1, corrected_group, True)
    raw_read_coverage2 = parse_fast5s(files2, corrected_group, True)

    chrm_sizes = dict((chrm, max(
        [r_data.end for r_data in raw_read_coverage1[chrm]] +
        [r_data.end for r_data in raw_read_coverage2[chrm]]))
                       for chrm in raw_read_coverage1)

    if VERBOSE: sys.stderr.write('Getting base signal.\n')
    base_signal1 = get_base_signal(raw_read_coverage1, chrm_sizes)
    base_signal2 = get_base_signal(raw_read_coverage2, chrm_sizes)

    if VERBOSE: sys.stderr.write(
            'Get differences between base signal.\n')
    # get num_region max diff regions from each chrm then find
    # global largest after
    largest_diff_indices = []
    for chrm, chrm_size in chrm_sizes.items():
        for strand in ('+', '-'):
            chrm_diffs = np.concatenate([
                np.abs(base_signal1[(chrm, strand)] -
                       base_signal2[(chrm, strand)])])
            chrm_max_diff_regs = np.argsort(
                chrm_diffs)[::-1][:num_regions]
            largest_diff_indices.extend((
                chrm_diffs[pos], max(pos - int(num_bases / 2.0), 0),
                chrm, strand) for pos in chrm_max_diff_regs)

    plot_intervals = zip(
        ['{:03d}'.format(rn) for rn in range(num_regions)],
        [(chrm, start, strand) for stat, start, chrm, strand in
         sorted(largest_diff_indices, reverse=True)[:num_regions]])

    plot_two_samples(
        plot_intervals, raw_read_coverage1, raw_read_coverage2,
        num_bases, overplot_thresh, overplot_type, corrected_group,
        pdf_fn)

    return

def compare_main(args):
    global VERBOSE
    VERBOSE = not args.quiet

    files1 = [os.path.join(args.fast5_basedir, fn)
              for fn in os.listdir(args.fast5_basedir)]
    files2 = [os.path.join(args.fast5_basedir2, fn)
              for fn in os.listdir(args.fast5_basedir2)]

    if DO_PROFILE:
        import cProfile
        cProfile.runctx(
            "plot_max_diff("
            "files1, files2, args.num_regions, " +
            "args.corrected_group, args.overplot_threshold, " +
            "args.pdf_filename, args.num_bases, " +
            "args.overplot_type)",
            globals(), locals(), 'profile.plot_compare.prof')
        sys.exit()
    if args.genome_locations is not None:
        plot_genome_locations(
            files1, files2, args.corrected_group,
            args.overplot_threshold, args.pdf_filename,
            args.num_bases, args.overplot_type, args.genome_locations)
    elif args.kmer:
        plot_kmer_centered(
            files1, files2, args.num_regions, args.corrected_group,
            args.overplot_threshold, args.pdf_filename,
            args.num_bases, args.overplot_type, args.kmer,
            args.fasta_filename, args.deepest_coverage)
    elif args.t_test_locations:
        # TODO test locations instead of difference in means
        pass
    else:
        plot_max_diff(
            files1, files2, args.num_regions, args.corrected_group,
            args.overplot_threshold, args.pdf_filename, args.num_bases,
            args.overplot_type)

    return

def get_compare_parser():
    import argparse
    parser = argparse.ArgumentParser(
        description='Plot raw signal from from two samples where ' +
        'FAST5 files were corrected with `nanoraw correct`.',
        add_help=False)
    parser.add_argument(
        'fast5_basedir',
        help='Directory containing fast5 files.')
    parser.add_argument(
        'fast5_basedir2',
        help='Second directory containing fast5 files.')

    parser.add_argument(
        '--pdf-filename',
        default='Nanopore_read_coverage.compare_groups.pdf',
        help='PDF filename to store plots. Default: %(default)s')

    parser.add_argument(
        '--t-test-locations', default=False, action='store_true',
        help="Select locations to plot based on t-test of " +
        "event means instead of difference in means.")

    return parser



########################################
#### General plotting option parser ####
########################################

def get_plot_parser():
    import argparse
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        '--num-regions', type=int, default=10,
        help='Number of regions to plot. Default: %(default)d')
    parser.add_argument(
        '--genome-locations', nargs='*',
        help='Plot genomic locations instead of regions selected ' +
        'by criterion (default is max coverage regions). Format ' +
        'locations as "chrm:position [chrm2:position2 ...]"')
    parser.add_argument(
        '--kmer',
        help='Plot regions centered on a particular kmer of ' +
        'interest. This option requires a FASTA also be provided. ' +
        'Default is to randomly select [--num_regions] of regions ' +
        'with this kmer. Set --deepest-coverage flag to select ' +
        'applicable regions for plotting')
    parser.add_argument(
        '--fasta-filename',
        help='FASTA file used to map reads with "correct" ' +
        'command (required for --kmer option).')
    parser.add_argument(
        '--deepest-coverage')

    parser.add_argument(
        '--num-bases', type=int, default=51,
        help='Number of bases to plot from region. Default: %(default)d')

    parser.add_argument(
        '--corrected-group', default='RawGenomeCorrected_000',
        help='FAST5 group to plot created by correct_raw ' +
        'script. Default: %(default)s')
    parser.add_argument(
        '--overplot-threshold', type=int, default=50,
        help='Number of reads to trigger alternative plot type ' +
        'instead of raw signal due to overplotting. ' +
        'Default: %(default)d')
    parser.add_argument(
        '--overplot-type', default='Boxplot',
        choices=['Boxplot', 'Quantile', 'Violin'],
        help='Plot type for regions with higher coverage. ' +
        'Choices: Boxplot (default), Quantile, Violin')

    parser.add_argument(
        '--quiet', '-q', default=False, action='store_true',
        help="Don't print status information.")

    return parser


if __name__ == '__main__':
    raise NotImplementedError, (
        'This is a module. Run with `nanoraw plot_signal -h`')
