import sys

import correct_raw
import plot_signal
import plot_kmer_quantiles

def main(args=None):
    """The main routine."""
    if args is None:
        args = sys.argv[1:]

    import argparse

    parser = argparse.ArgumentParser(
        prog='nanoraw',
        description='nanoraw is a python toolset to analyze and ' +
        'visualize raw nanopore sequencing data.')
    subparsers = parser.add_subparsers()

    # create the parser for the "correct" command
    correct_parser = correct_raw.get_parser(False)
    subparser_correct = subparsers.add_parser(
        'correct', parents=[correct_parser,],
        help='Correct annotation of raw signal with genomic ' +
        'aignement of existing basecalls.')
    subparser_correct.set_defaults(func=correct_raw.main)

    # create the parser for the "plot_signal" command
    signal_parser = plot_signal.get_parser(False)
    subparser_signal = subparsers.add_parser(
        'plot_signal', parents=[signal_parser,],
        help='Plot signal across selected genomic regions.')
    subparser_signal.set_defaults(func=plot_signal.main)

    # create the parser for the "plot_kmer" command
    kmer_parser = plot_kmer_quantiles.get_parser(False)
    subparser_kmer = subparsers.add_parser(
        'plot_kmer', parents=[kmer_parser,],
        help='Plot signal quantiles acorss kmers.')
    subparser_kmer.set_defaults(func=plot_kmer_quantiles.main)

    args = parser.parse_args(args)

    args.func(args)

    return

if __name__ == "__main__":
    main()
