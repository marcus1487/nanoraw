Auxiliary Commands
******************

`write_most_significant_fasta`
------------------------------

Write sequences in FASTA format where signal differs the most significantly between two groups. This subcommand is suggested for running downstream motif discovery pipelines. Options are similar or equivalent to significance plotting commands.

`write_wiggles`
---------------

Write wiggle files for attributes from a group of FAST5 files. The currenly available attributes are `coverage`, `signal`, `signal_sd`, `length` (mean over all reads at each base for last three), and if two FAST5 directories or a statistics file are provided `pvals` (Negative log 10), `qvals` (Negative log 10), and `difference`.

`write_wiggles` Options
+++++++++++++++++++++++

- `--wiggle-types`: Specify which wiggle files should be produced. Options are listed above. (Default: "coverage, pval")
- `--wiggle-basename`: Basename for output files. Two files (plus and minus strand) will be produced for each --wiggle-types supplied (four files will be produced if `--fast5-basedirs2` is supplied for coverage, signal, signal_sd, and length). Filenames will be formatted as "[wiggle-basename].[wiggle-type].[group1/2]?.[plus/minus].wig". Default: Nanopore_data

Example commands
----------------

Write FASTA and wiggles examples::
  
  nanoraw write_most_significant_fasta --fast5-basedirs $g1Dir \
        --fast5-basedirs2 $g2Dir --2d \
        --genome-fasta $genomeFn
  nanoraw write_wiggles --fast5-basedirs $g1Dir \
        --fast5-basedirs2 $g2Dir --2d \
        --wiggle-types coverage signal signal_sd length pvals \
        qvals difference \
        --statistics-filename testing.significance_values.txt
