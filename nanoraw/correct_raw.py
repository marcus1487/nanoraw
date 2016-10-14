import os, sys

import h5py
import Queue

import numpy as np
import multiprocessing as mp

from glob import glob
from time import sleep, time
from subprocess import call, STDOUT
from itertools import groupby, izip
from tempfile import NamedTemporaryFile
from collections import defaultdict, namedtuple

from nanoraw_helper import normalize_raw_signal

fd = sys.stderr.fileno()
def _redirect_stderr(to):
    sys.stderr.close()
    os.dup2(to.fileno(), fd)
    sys.stderr = os.fdopen(fd, 'w')

with os.fdopen(os.dup(fd), 'w') as old_stderr:
    with open(os.devnull, 'w') as fp:
        _redirect_stderr(fp)
    try:
        # load changepoint from R since pythons isn't very stable
        import rpy2.robjects as r
        from rpy2.robjects.packages import importr
        rChpt = importr("changepoint")
        USE_R_CPTS = True
    except:
        USE_R_CPTS = False
    finally:
        _redirect_stderr(old_stderr)


NANORAW_VERSION = '0.1'
DO_PROFILE = False
VERBOSE = False

COMP_BASES = {'A':'T', 'C':'G', 'G':'C', 'T':'A', '-':'-'}
def rev_comp(seq):
    return ''.join(COMP_BASES[b] for b in seq[::-1])

BASES = ['A', 'C', 'G', 'T']
baseP = ['p_A', 'p_C', 'p_G', 'p_T']

indelStats = namedtuple('indelStats',
                        ('start', 'end', 'diff'))
indelGroupStats = namedtuple('indelGroupStats',
                             ('start', 'end', 'cpts', 'indels'))
channelInfo = namedtuple(
    'channelInfo',
    ('offset', 'range', 'digitisation', 'number', 'sampling_rate'))
readInfo = namedtuple(
    'readInfo',
    ('ID', 'Start', 'TrimStart', 'TrimEnd',
     'Insertions', 'Deletions', 'Matches', 'Mismatches'))
genomeLoc = namedtuple(
    'genomeLoc', ('Start', 'Strand', 'Chrom'))

GRAPHMAP_FIELDS = (
    'qName', 'qLength', 'qStart', 'qEnd', 'qStrand',
    'tName', 'tLength', 'tStart', 'tEnd', 'tStrand',
    'score', 'numMatch', 'numMismatch', 'numIns', 'numDel',
    'mapQV', 'qAlignedSeq', 'matchPattern', 'tAlignedSeq')


def write_new_fast5_group(
        filename, genome_location, read_info,
        read_start_rel_to_raw, new_segs, align_seq, alignVals,
        norm_signal, shift, scale, basecall_group, corrected_group):
    # save new events as new hdf5 Group
    read_data = h5py.File(filename, 'r+')
    analyses_grp = read_data['Analyses']
    # already checked if group existed and not overwrite option
    # before running correction function, so shouldn't need to
    # check here (may be an issue if script is run in
    # multiple simultaneous instances)
    if corrected_group in analyses_grp:
        del analyses_grp[corrected_group]
    corr_grp = analyses_grp.create_group(corrected_group)
    corr_grp.attrs['version'] = NANORAW_VERSION
    corr_grp.attrs['basecall_group'] = basecall_group

    # store alignment statistics
    corr_alignment = corr_grp.create_group('Alignment')
    corr_alignment.attrs['mapped_start'] = genome_location.Start
    corr_alignment.attrs['mapped_strand'] = genome_location.Strand
    corr_alignment.attrs['mapped_chrom'] = genome_location.Chrom

    corr_alignment.attrs['trimmed_obs_start'] = read_info.TrimStart
    corr_alignment.attrs['trimmed_obs_end'] = read_info.TrimEnd
    corr_alignment.attrs['num_insertions'] = read_info.Insertions
    corr_alignment.attrs['num_deletions'] = read_info.Deletions
    corr_alignment.attrs['num_matches'] = read_info.Matches
    corr_alignment.attrs['num_mismatches'] = read_info.Mismatches

    np_read_align = np.chararray(len(alignVals))
    np_read_align[:] = zip(*alignVals)[0]
    corr_alignment.create_dataset('read_alignment', data=np_read_align)
    np_genome_align = np.chararray(len(alignVals))
    np_genome_align[:] = zip(*alignVals)[1]
    corr_alignment.create_dataset(
        'genome_alignment', data=np_genome_align)

    # store new aligned segments and sequence
    corr_temp = corr_grp.create_group('template')
    corr_temp.attrs['shift'] = shift
    corr_temp.attrs['scale'] = scale
    corr_events = corr_temp.create_dataset('Segments', data=new_segs)
    corr_events.attrs['read_start_rel_to_raw'] = read_start_rel_to_raw

    # Add Events to data frame with event means, SDs and lengths
    """raw_mean_sd = [(base_sig.mean(), base_sig.std()) for base_sig in
                   np.split(raw_signal, new_segs[1:-1])]"""
    norm_mean_sd = [(base_sig.mean(), base_sig.std()) for base_sig in
                    np.split(norm_signal, new_segs[1:-1])]

    event_data = np.array(
        zip(zip(*norm_mean_sd)[0], zip(*norm_mean_sd)[1],
            new_segs[:-1], np.diff(new_segs), list(align_seq)),
        dtype=[('norm_mean', '<f8'), ('norm_stdev', '<f8'),
               ('start', '<f8'), ('length', '<f8'), ('base', 'S1')])
    corr_events = corr_temp.create_dataset('Events', data=event_data)

    read_data.flush()
    read_data.close()

    return

def get_indel_groups(
        align_dat, align_segs, raw_signal, min_seg_len, timeout,
        num_cpts_limit):
    def get_deletion_stats():
        read_seq = ''.join(zip(*align_dat['aligned'])[0])
        del_stats = []
        # group contiguous deleted bases
        for k, g in groupby(enumerate(
                [j for j, b in enumerate(read_seq) if b == "-"]),
                            lambda (i, x): i-x):
            g = list(g)
            # store deletion start position, stop position and
            # total bases effect on final sequence
            del_stats.append(indelStats(
                max(0 ,g[0][1] - 1),
                min(g[0][1] + len(g) + 1, len(align_segs) - 1),
                -1 * len(g)))
        return del_stats
    def get_insertion_stats():
        # similar to deletion stats store start, stop and bases effect
        return [
            indelStats(max(0, event_num - 1),
                       min(event_num + 1, len(align_segs)), len(bs))
            for event_num, bs in align_dat['insertion']]
    def extend_group(indel_group):
        group_start = min(indel.start for indel in indel_group)
        group_stop = max(indel.end for indel in indel_group)
        num_cpts = sum(indel.diff for indel in indel_group
                       ) + group_stop - group_start - 1
        # check that there are enough points to split
        # add an extra set of values to ensure no zero changepoint
        while align_segs[group_stop] - align_segs[group_start] < (
                num_cpts + 2) * min_seg_len:
            num_cpts += int(group_start > 0) + int(
                group_stop < len(align_segs) - 1)
            group_start = max(0, group_start - 1)
            group_stop = min(len(align_segs) - 1, group_stop + 1)
        return group_start, group_stop, num_cpts
    def extend_and_join(indel_group):
        group_start, group_stop, num_cpts = extend_group(indel_group)
        # check if the extension hits the previous group
        while (len(indel_groups) > 0) and (
                group_start <= indel_groups[-1].end):
            indel_group = indel_groups[-1].indels + indel_group
            del indel_groups[-1]
            group_start, group_stop, num_cpts = extend_group(indel_group)
        return group_start, group_stop, num_cpts, indel_group
    def get_cpts_alt(group_start, group_stop, num_cpts):
        """
        Get changepoints where the raw difference between min_seg_len
        obs to the left and min_seg_len obs to the right is largest
        while maintaining the min_seg_len between changepoints.
        Still need to test this function for off by one bugs etc.
        """
        if num_cpts_limit is not None and num_cpts > num_cpts_limit:
            raise RuntimeError, ('Reached maximum number of ' +
                                 'changepoints for a single indel')
        sig_cs = raw_signal[align_segs[group_start]:
                            align_segs[group_stop]]
        sig_cs = np.cumsum(np.insert(sig_cs, 0, 0))
        # get difference between all neighboring min_seg_len regions
        running_diffs = np.abs((2 * sig_cs[min_seg_len:-min_seg_len]) -
                               sig_cs[:-2*min_seg_len] -
                               sig_cs[2*min_seg_len:])
        cpts = []
        blacklist_pos = set()
        for pos in running_diffs.argsort()[::-1]:
            if pos not in blacklist_pos:
                cpts.append(pos)
                blacklist_pos.update(range(
                    pos-min_seg_len+1, pos+min_seg_len+1))
            if len(cpts) == num_cpts:
                break
        if len(cpts) < num_cpts:
            return None
        return sorted([cpt + min_seg_len for cpt in cpts])
    def get_cpts_R(group_start, group_stop, num_cpts):
        if num_cpts_limit is not None and num_cpts > num_cpts_limit:
            raise RuntimeError, ('Reached maximum number of ' +
                                 'changepoints for a single indel')
        cpts = rChpt.cpt_mean(
            r.FloatVector(raw_signal[
                align_segs[group_start]:align_segs[group_stop]]),
                method="BinSeg", Q=num_cpts, penalty="None",
            minseglen=min_seg_len).do_slot('cpts')
        return [int(x) for x in cpts[:-1]]
    get_cpts = get_cpts_R if USE_R_CPTS else get_cpts_alt
    def extend_for_cpts(group_start, group_stop, num_cpts, indel_group):
        cpts = get_cpts(group_start, group_stop, num_cpts)
        # There is a bug in the changepoint package that allows a zero
        # width first segment. If one is present extend the region and
        # find cpts again
        # from rpy2 version
        #while cpts[0] == 0:
        while cpts is None:
            num_cpts += int(group_start > 0) + int(
                group_stop < len(align_segs) - 1)
            group_start = max(0, group_start - 1)
            group_stop = min(len(align_segs) - 1, group_stop + 1)
            if (len(indel_groups) > 0) and (
                    group_start <= indel_groups[-1].end):
                indel_group = indel_groups[-1].indels + indel_group
                del indel_groups[-1]
                group_start, group_stop, num_cpts = extend_group(
                    indel_group)
            cpts = get_cpts(group_start, group_stop, num_cpts)

        return [x + align_segs[group_start]
                for x in cpts], group_start, group_stop, indel_group

    if timeout is not None:
        timeout_start = time()

    # sort indels in order of start positions
    all_indels = sorted(get_insertion_stats() + get_deletion_stats())
    indel_groups = []
    curr_group = [all_indels[0],]
    for indel in all_indels[1:]:
        if timeout is not None and time() - timeout_start > timeout:
            raise RuntimeError, 'Read took too long to re-segment.'
        # check if indel hits current group
        if max(indel.end for indel in curr_group) >= indel.start:
            curr_group.append(indel)
        else:
            (curr_start, curr_stop, num_cpts,
             curr_group) = extend_and_join(curr_group)
            cpts, curr_start, curr_stop, curr_group = extend_for_cpts(
                curr_start, curr_stop, num_cpts, curr_group)
            # if the indel group still reaches the next indel
            if curr_stop >= indel.start:
                curr_group.append(indel)
            else:
                indel_groups.append(indelGroupStats(
                    curr_start, curr_stop, cpts, curr_group))
                curr_group = [indel,]

    # handle the last indel group if it is not yet included
    if indel_groups[-1].indels[-1] != all_indels[-1]:
        curr_start, curr_stop, num_cpts, curr_group = extend_and_join(
            curr_group)
        cpts, curr_start, curr_stop, curr_group = extend_for_cpts(
            curr_start, curr_stop, num_cpts, curr_group)
        indel_groups.append(indelGroupStats(
            curr_start, curr_stop, cpts, curr_group))

    return indel_groups

def get_aligned_seq(align_dat):
    # get read sequence with deletions and then sequence
    # after fixing deletions and insertions
    read_seq = ''.join(zip(*align_dat['aligned'])[0])
    # if there are insertions add those to the sequence
    if len(align_dat['insertion']) > 0:
        read_seq = ''.join([
            read_seq[i:j] + bs for i, (j, bs) in zip(
                    [0,] + sorted(
                        zip(*align_dat['insertion'])[0]),
                    sorted(align_dat['insertion']) +
                    [(len(read_seq), ''),])])
    return read_seq.replace("-", "")

def align_to_genome(basecalls, genome_filename, graphmap_path):
    ## align to genome
    read_fp = NamedTemporaryFile(
        delete=False, suffix='.fasta')
    read_fp.write(">" + ''.join(basecalls[:5]) + '\n' +
                  ''.join(basecalls) + '\n')
    read_fp.close()
    out_fp = NamedTemporaryFile(delete=False)
    try:
        # suppress output from graphmap with devnull sink
        with open(os.devnull, 'w') as FNULL:
            exitStatus = call([
                graphmap_path, 'align',
                '-r', genome_filename,
                '-d', read_fp.name,
                '-o', out_fp.name,
                '-L', 'm5'], stdout=FNULL, stderr=STDOUT)

        alignment = dict(zip(GRAPHMAP_FIELDS, out_fp.readline().split()))
        out_fp.close()
    except:
        raise OSError, ('Problem running/parsing graphmap. Ensure ' +
                        'you have a compatible version installed.')

    ## flip read and genome to match raw signal if neg strand mapped
    if len(alignment) != len(GRAPHMAP_FIELDS):
        raise NotImplementedError, (
            'Graphmap did not produce alignment.')

    if alignment['tStrand'] != '+':
        raise NotImplementedError, (
            'Graphmap indicates negative strand reference mapping.')

    if alignment['qStrand'] == "+":
        alignVals = zip(alignment['qAlignedSeq'],
                        alignment['tAlignedSeq'])
    else:
        alignVals = zip(rev_comp(alignment['qAlignedSeq']),
                        rev_comp(alignment['tAlignedSeq']))

    return alignVals, genomeLoc(
        int(alignment['tStart']), alignment['qStrand'],
        alignment['tName'])

def trim_alignment(alignVals, starts_rel_to_read,
                   read_start_rel_to_raw, abs_event_start,
                   genome_location):
    # trim read to first matching bases
    start_trim_align_pos = 0
    start_trim_read_pos = 0
    start_trim_genome_pos = 0
    r_base, g_base = alignVals[0]
    while r_base == '-' or g_base == '-':
        start_trim_read_pos += int(r_base != '-')
        start_trim_genome_pos += int(g_base != '-')
        start_trim_align_pos += 1
        r_base, g_base = alignVals[start_trim_align_pos]

    end_trim_align_pos = 0
    end_trim_read_pos = 0
    end_trim_genome_pos = 0
    r_base, g_base = alignVals[-1]
    while r_base == '-' or g_base == '-':
        end_trim_read_pos += int(r_base != '-')
        end_trim_genome_pos += int(g_base != '-')
        end_trim_align_pos += 1
        r_base, g_base = alignVals[-1 * (end_trim_align_pos + 1)]

    alignVals = alignVals[start_trim_align_pos:]
    if end_trim_align_pos > 0:
        alignVals = alignVals[:-1*end_trim_align_pos]

    if start_trim_read_pos > 0:
        start_trim_obs = starts_rel_to_read[start_trim_read_pos]
        starts_rel_to_read = starts_rel_to_read[
            start_trim_read_pos:] - start_trim_obs
        read_start_rel_to_raw += start_trim_obs
        abs_event_start += start_trim_obs

    if end_trim_read_pos > 0:
        starts_rel_to_read = starts_rel_to_read[:-1 * end_trim_read_pos]

    if genome_location.Strand == '+' and start_trim_genome_pos > 0:
        genome_location = genomeLoc(
            genome_location.Start + start_trim_genome_pos, '+',
            genome_location.Chrom)
    elif genome_location.Strand == '-' and end_trim_genome_pos > 0:
        genome_location = genomeLoc(
            genome_location.Start + end_trim_genome_pos, '-',
            genome_location.Chrom)

    return (alignVals, starts_rel_to_read, read_start_rel_to_raw,
            abs_event_start, start_trim_read_pos, end_trim_read_pos,
            genome_location)

def fix_stay_states(
        called_dat, starts_rel_to_read, basecalls,
        read_start_rel_to_raw, abs_event_start):
    change_states = (called_dat['model_state'][:-1] !=
                     called_dat['model_state'][1:])
    start_trim = 0
    event_change_state = change_states[0]
    while not event_change_state:
        if start_trim >= len(change_states) - 2:
            raise RuntimeError, (
                'Read is composed entirely of stay model ' +
                'states and cannot be processed')
        start_trim += 1
        event_change_state = change_states[start_trim]
    end_trim = 0
    event_change_state = change_states[-1]
    while not event_change_state:
        end_trim += 1
        event_change_state = change_states[-(end_trim+1)]

    # trim all applicable data structures
    change_states = change_states[start_trim:]
    starts_rel_to_read = starts_rel_to_read[start_trim:]
    basecalls = basecalls[start_trim:]
    if end_trim > 0:
        change_states = change_states[:-end_trim]
        starts_rel_to_read = starts_rel_to_read[:-end_trim]
        basecalls = basecalls[:-end_trim]
    if start_trim > 0:
        start_trim_obs = starts_rel_to_read[0]
        starts_rel_to_read = starts_rel_to_read - start_trim_obs
        read_start_rel_to_raw += start_trim_obs
        abs_event_start += start_trim_obs

    # now actually remove internal stay states
    change_states = np.insert(
        change_states, (0, len(change_states) - 1), True)
    starts_rel_to_read = starts_rel_to_read[change_states]
    basecalls = basecalls[change_states[:-1]]

    return (starts_rel_to_read, basecalls, read_start_rel_to_raw,
            abs_event_start, start_trim, end_trim)

def get_read_data(filename, rmStayStates, basecall_group):
    try:
        fast5_data = h5py.File(filename, 'r')
    except IOError:
        raise IOError, 'Error opening file. Likely a corrupted file.'

    try:
        channel_info = fast5_data['UniqueGlobalKey/channel_id'].attrs
    except:
        raise RuntimeError, ("No channel_id group in HDF5 file. " +
                             "Probably mux scan HDF5 file.")
    channel_info = channelInfo(
        channel_info['offset'], channel_info['range'],
        channel_info['digitisation'], channel_info['channel_number'],
        channel_info['sampling_rate'].astype('int_'))

    # TODO: support complement and template event correction
    try:
        called_dat = fast5_data[
            'Analyses/' + basecall_group + '/BaseCalled_template/Events']
    except:
        raise RuntimeError, ('No events in file. Likely a ' +
                             'segmentation error.')
    raw_dat = fast5_data['Raw/Reads/'].values()[0]
    all_raw_signal = raw_dat['Signal'].value
    read_id = raw_dat.attrs['read_id']

    abs_event_start = int(called_dat.attrs['start_time'] *
                          channel_info.sampling_rate)
    read_start_rel_to_raw = int(
        abs_event_start - raw_dat.attrs['start_time'])

    last_event = called_dat[-1]
    starts_rel_to_read = np.append(
        called_dat['start'], last_event['start'] + last_event['length'])
    starts_rel_to_read = (
        starts_rel_to_read *
        channel_info.sampling_rate).astype('int_') - abs_event_start
    basecalls = np.array([
        BASES[eventBasePs.argmax()]
        for eventBasePs in np.column_stack(
            called_dat[basePVal] for basePVal in baseP)])

    if any(len(vals) <= 1 for vals in (
        all_raw_signal, starts_rel_to_read, basecalls,
        called_dat['model_state'])):
        raise NotImplementedError, (
            'One or no segments or signal present in read.')

    start_trim, end_trim = 0, 0
    if rmStayStates:
        (starts_rel_to_read, basecalls, read_start_rel_to_raw,
         abs_event_start, start_trim, end_trim) = fix_stay_states(
             called_dat, starts_rel_to_read, basecalls,
             read_start_rel_to_raw, abs_event_start)

    fast5_data.close()

    return (all_raw_signal, read_start_rel_to_raw, starts_rel_to_read,
            basecalls, channel_info, abs_event_start, read_id,
            start_trim, end_trim)

def correct_raw_data(
        filename, genome_filename, graphmap_path, basecall_group,
        corrected_group, rmStayStates=True, outlier_threshold=5,
        timeout=None, min_event_obs=4, num_cpts_limit=None,
        overwrite=True, in_place=True):
    if not overwrite and in_place:
        try:
            read_data = h5py.File(filename, 'r')
            if 'Analyses/' + corrected_group in read_data:
                raise RuntimeError, (
                    "Raw genome corrected data exists for " +
                    "this read and --overwrite is not set.")
            read_data.close()
        except IOError:
            raise IOError, 'Error opening file. Likely a corrupted file.'

        raise NotImplementedError, (
            'Not currently implimenting new hdf5 file writing.')
    else:
        # check that the file is writeable before trying to correct
        if not os.access(filename, os.W_OK):
            raise IOError, 'FAST5 file is not writable.'

    (all_raw_signal, read_start_rel_to_raw, starts_rel_to_read,
     basecalls, channel_info, abs_event_start, read_id,
     rm_stay_start_trim, rm_stay_end_trim) = get_read_data(
         filename, rmStayStates, basecall_group)

    alignVals, genome_location = align_to_genome(
        basecalls, genome_filename, graphmap_path)

    (alignVals, starts_rel_to_read, read_start_rel_to_raw,
     abs_event_start, startTrim, endTrim,
     genome_location) = trim_alignment(
         alignVals, starts_rel_to_read, read_start_rel_to_raw,
         abs_event_start, genome_location)

    read_info = readInfo(
        read_id, abs_event_start / float(channel_info.sampling_rate),
        startTrim, endTrim,
        sum(b == "-" for b in zip(*alignVals)[0]),
        sum(b == "-" for b in zip(*alignVals)[1]),
        sum(rb == gb for rb, gb in alignVals),
        sum(rb != gb for rb, gb in alignVals
            if rb != '-' and gb != '-'))

    # normalize signal
    norm_signal, shift, scale = normalize_raw_signal(
        all_raw_signal, read_start_rel_to_raw, starts_rel_to_read[-1],
        'median', channel_info, outlier_threshold)

    # parse out aligned (including deleted) and inserted bases
    align_dat = {'aligned':[], 'insertion':defaultdict(str)}
    for r_base, g_base in alignVals:
        if r_base == "-":
            align_dat['insertion'][
                len(align_dat['aligned'])] += g_base
        else:
            align_dat['aligned'].append((
                g_base, g_base == r_base))
    align_dat['insertion'] = align_dat['insertion'].items()

    # group indels that are adjacent for re-segmentation
    indel_groups = get_indel_groups(
        align_dat, starts_rel_to_read, norm_signal, min_event_obs,
        timeout, num_cpts_limit)

    new_segs = []
    prev_stop = 0
    for group_start, group_stop, cpts, group_indels in indel_groups:
        ## add segments from last indel to this one and new segments
        new_segs.append(
            np.append(starts_rel_to_read[prev_stop:group_start+1], cpts))
        prev_stop = group_stop
    # handle end of read
    new_segs.append(starts_rel_to_read[prev_stop:])
    new_segs = np.concatenate(new_segs)

    align_seq = get_aligned_seq(align_dat)
    if new_segs.shape[0] != len(align_seq) + 1:
        raise ValueError, ('Aligned sequence does not match number ' +
                           'of segments produced.')

    if in_place:
        write_new_fast5_group(
            filename, genome_location, read_info,
            read_start_rel_to_raw, new_segs, align_seq, alignVals,
            norm_signal, shift, scale, basecall_group, corrected_group)
    else:
        # create new hdf5 file to hold new read signal
        raise NotImplementedError, (
            'Not currently implementing new hdf5 file writing.')

    return

def get_aligned_seq_worker(
        filenames_q, failed_reads_q, genome_filename, graphmap_path,
        basecall_group, corrected_group, timeout, num_cpts_limit,
        overwrite):
    if DO_PROFILE:
        import cProfile
        cProfile.runctx(
            """while not filenames_q.empty():
        try:
            fn = filenames_q.get(block=False)
        except Queue.Empty:
            break

        try:
            correct_raw_data(
                fn, genome_filename, graphmap_path, basecall_group,
                corrected_group, timeout=timeout,
                num_cpts_limit=num_cpts_limit)
        except Exception as e:
            failed_reads_q.put((str(e), fn))""",
            globals(), locals(), 'profile.correct_reads.prof')
        sys.exit()
    num_processed = 0
    while not filenames_q.empty():
        try:
            fn = filenames_q.get(block=False)
        except Queue.Empty:
            break

        num_processed += 1
        if VERBOSE and num_processed % 100 == 0:
            sys.stderr.write('.')
            sys.stderr.flush()

        try:
            correct_raw_data(
                fn, genome_filename, graphmap_path, basecall_group,
                corrected_group, timeout=timeout,
                num_cpts_limit=num_cpts_limit, overwrite=overwrite)
        except Exception as e:
            failed_reads_q.put((str(e), fn))

    return

def get_all_reads_data(fast5_fns, genome_filename, graphmap_path,
                       basecall_group, corrected_group,
                       num_processes=1, timeout=None,
                       num_cpts_limit=None, overwrite=True):
    manager = mp.Manager()
    fast5_q = manager.Queue()
    failed_reads_q = manager.Queue()
    num_reads = 0
    for fast5_fn in fast5_fns:
        num_reads += 1
        fast5_q.put(fast5_fn)

    args = (fast5_q, failed_reads_q, genome_filename,
            graphmap_path, basecall_group, corrected_group, timeout,
            num_cpts_limit, overwrite)
    processes = []
    for p_id in xrange(num_processes):
        p = mp.Process(target=get_aligned_seq_worker, args=args)
        p.start()
        processes.append(p)

    if VERBOSE: sys.stderr.write(
            'Correcting ' + str(num_reads) + ' reads (Will ' +
            'print a dot for each 100 reads completed).\n')
    failed_reads = defaultdict(list)
    while any(p.is_alive() for p in processes):
        try:
            errorType, fn = failed_reads_q.get(block=False)
            failed_reads[errorType].append(fn)
        except Queue.Empty:
            sleep(1)
            continue

    # empty any entries left in queue after processes have finished
    while not failed_reads_q.empty():
        errorType, fn = failed_reads_q.get(block=False)
        failed_reads[errorType].append(fn)

    # print newline after read progress dots
    if VERBOSE: sys.stderr.write('\n')

    return dict(failed_reads)

def main(args):
    global VERBOSE
    VERBOSE = not args.quiet

    global USE_R_CPTS
    if not USE_R_CPTS and args.use_r_cpts:
        sys.stderr.write(
            'WARNING: Could not load rpy2 (inlucding package ' +
            '"changepoint"). Using python segmentation method.')
    USE_R_CPTS = args.use_r_cpts and USE_R_CPTS

    num_p = 1 if DO_PROFILE else args.processes

    if VERBOSE: sys.stderr.write('Getting file list.\n')
    if args.fast5_pattern:
        files = glob(os.path.join(
            args.fast5_basedir, args.fast5_pattern))
    else:
        files = [os.path.join(args.fast5_basedir, fn)
                 for fn in os.listdir(args.fast5_basedir)]

    failed_reads = get_all_reads_data(
        files, args.genome_fasta, args.graphmap_path,
        args.basecall_group, args.corrected_group, num_p, args.timeout,
        args.cpts_limit, args.overwrite)
    sys.stderr.write('Failed reads summary:\n' + '\n'.join(
        "\t" + err + " :\t" + str(len(fns))
        for err, fns in failed_reads.items()) + '\n')
    if args.failed_reads_filename is not None:
        with open(args.failed_reads_filename, 'w') as fp:
            fp.write('\n'.join((
                err + '\t' + ','.join(fns)
                for err, fns in failed_reads.items())) + '\n')

    return

# define function for getting parser so it can be shared in
# __main__ package script
def get_parser():
    import argparse
    parser = argparse.ArgumentParser(
        description='Parse raw data from oxford nanopore R9 FAST5 ' +
        'files and correct segmentation to match genomic regions.',
        add_help=False)

    parser.add_argument(
        'fast5_basedir',
        help='Directory containing fast5 files.')
    parser.add_argument(
        'graphmap_path',
        help='Full path to built graphmap version.')
    parser.add_argument(
        'genome_fasta',
        help='Path to fasta file for mapping.')

    parser.add_argument(
        '--fast5-pattern',
        help='A pattern to search for a subset of files within ' +
        'fast5-basedir. Note that on the unix command line patterns ' +
        'may be expanded so it is best practice to quote patterns.')

    parser.add_argument(
        '--processes', default=1, type=int,
        help='Number of processes. Default: %(default)d')
    parser.add_argument(
        '--timeout', default=None, type=int,
        help='Timeout in seconds for the processing of a single ' +
        'read.  Default: No timeout.')
    parser.add_argument(
        '--cpts-limit', default=None, type=int,
        help='Maximum number of changepoints to find within a single ' +
        'indel group. (Not setting this option can cause a process ' +
        'to stall and cannot be controlled by the timeout option). ' +
        'Default: No limit.')

    parser.add_argument(
        '--basecall-group', default='Basecall_1D_000',
        help='FAST5 group to use for obtaining original ' +
        'basecalls (under Analyses group). Default: %(default)s')
    parser.add_argument(
        '--corrected-group', default='RawGenomeCorrected_000',
        help='FAST5 group to use for storing corrected ' +
        'events produced by this script. Default: %(default)s')

    parser.add_argument(
        '--use-r-cpts', default=False, action='store_true',
        help='Use R changepoint package to determine new event ' +
        'delimiters. (requires rpy2, R and R package ' +
        '"changepoint" to be installed)')
    parser.add_argument(
        '--failed-reads-filename',
        help='Output failed read filenames into a this file ' +
        'with assoicated error for each read. Default: ' +
        'Do not store failed reads.')

    parser.add_argument(
        '--overwrite', default=False, action='store_true',
        help='Overwrite previous corrected group in FAST5/HDF5 ' +
        'file. (Note this only effects the group defined ' +
        'by --corrected-group).')

    parser.add_argument(
        '--quiet', '-q', default=False, action='store_true',
        help="Don't print status information.")

    return parser

def args_and_main():
    main(get_parser().parse_args())
    return

if __name__ == '__main__':
    args_and_main()
