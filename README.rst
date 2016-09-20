## Summary

This package provides tools for the analysis of raw nanopore sequencing data, including correction of basecalls and visualization.

## Requirements

- HDF5
- python
- h5py
- numpy

Optionally:
R
- changepoint (for using R's changepoint package for re-segmentation)
- ggplot2 (for plotting scripts)

## Usage

Currently, not loaded into pip (to be done soon). To use nanoraw install locally and run with:
```
git clone https://github.com/marcus1487/nanoraw.git
cd nanoraw
pip install .
nanoraw -h
nanoraw correct -h
nanoraw plot_signal -h
nanoraw plot_kmer -h
```
