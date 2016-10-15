import sys

import correct_raw
import analyze_coverage
import plot_commands

from collections import OrderedDict

def main(args=None):
    """The main routine."""
    if args is None:
        args = sys.argv[1:]

    commands = OrderedDict(
        [('Main comands:',OrderedDict([
            ('correct  ','Correct annotation of raw signal with ' +
            'genomic aignement of existing basecalls'),
            ('write_wiggle','Write wiggle file of genome coverage.')
        ])),
        ('Plotting comands:',OrderedDict([
            ('plot_signal','Plot signal across selected ' +
             'genomic regions.'),
            ('plot_comparison','Plot comparison of two read groups.'),
            ('plot_correction','Plot segmentation before and after correction.'),
            ('plot_kmer','Plot signal quantiles acorss kmers.'),
        ]))])
    raw_cmds = [cmd.strip() for grp_cmds in commands.values()
                for cmd in grp_cmds]
    desc = '\n\n'.join([
        grp + '\n' + '\n'.join([
            '\t' + cmd + "\t" + cmd_help
            for cmd, cmd_help in cmds.items()])
        for grp, cmds in commands.items()])

    import argparse
    parser = argparse.ArgumentParser(
        prog='nanoraw',
        description='nanoraw is a python toolset to analyze and ' +
        'visualize raw nanopore sequencing data.',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(
        title='commands', description=desc,
        help='Additional help available for subcommands.')

    # create the parser for the "correct" command
    correct_parser = correct_raw.get_parser()
    subparser_correct = subparsers.add_parser(
        raw_cmds[0], parents=[correct_parser,])
    subparser_correct.set_defaults(func=correct_raw.main)

    # create the parser for the "write_wiggle" command
    wiggle_parser = analyze_coverage.get_wiggle_parser()
    subparser_wiggle = subparsers.add_parser(
        raw_cmds[1], parents=[wiggle_parser,])
    subparser_wiggle.set_defaults(func=analyze_coverage.wiggle_main)


    # plotting functions
    plot_parser = plot_commands.get_plot_parser()

    # create the parser for the "plot_signal" command
    signal_parser = plot_commands.get_single_sample_parser()
    subparser_signal = subparsers.add_parser(
        raw_cmds[2], parents=[plot_parser, signal_parser])
    subparser_signal.set_defaults(func=plot_commands.single_sample_main)

    compare_parser = plot_commands.get_compare_parser()
    subparser_compare = subparsers.add_parser(
        raw_cmds[3], parents=[plot_parser, compare_parser])
    subparser_compare.set_defaults(func=plot_commands.compare_main)

    plot_correct_parser = plot_commands.get_correction_parser()
    subparser_plot_correct = subparsers.add_parser(
        raw_cmds[4], parents=[plot_correct_parser,])
    subparser_plot_correct.set_defaults(
        func=plot_commands.plot_correction_main)

    # create the parser for the "plot_kmer" command
    kmer_parser = plot_commands.get_kmer_parser()
    subparser_kmer = subparsers.add_parser(
        raw_cmds[5], parents=[kmer_parser,])
    subparser_kmer.set_defaults(func=plot_commands.kmer_main)

    args = parser.parse_args(args)

    args.func(args)

    return

if __name__ == "__main__":
    main()
