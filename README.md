## Summary
This package provides tools for the analysis of raw nanopore sequencing data, including correction of basecalls and visualization.

## Requirements

- HDF5 (<http://micro.stanford.edu/wiki/Install_HDF5#Install>)
- graphmap (<https://github.com/isovic/graphmap>)

#### python Requirements:
- python
- numpy
- h5py
- rpy2

#### R Requirements:
- R
- ggplot2 (not actually enfored on install, but plotting commands will fail; install with `install.packages('ggplot2')` from an R prompt)

#### Optionally:
- changepoint (for using R's changepoint package for re-segmentation)

## Usage
Currently, not loaded into pip (to be done soon). To use nanoraw install locally and run with:
```
git clone https://github.com/marcus1487/nanoraw.git
cd nanoraw
pip install .
nanoraw -h
nanoraw [command] [options]
```

## Currently supported commands
#### Main comands:
- correct: Correct annotation of raw signal with genomic aignement of existing basecalls
- write_wiggle: Write wiggle file of genome coverage.

#### Plotting comands:
- plot_signal: Plot signal across selected genomic regions.
- plot_comparison: Plot comparison of two read groups.
- plot_kmer: Plot signal quantiles acorss kmers.

> Get additional help for subcommands with `nanoraw [command] -h`
