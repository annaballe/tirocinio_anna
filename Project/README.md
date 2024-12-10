Project notes
==============

Issue di Truvari 

📝[Truvari issue 212](https://github.com/ACEnglish/truvari/issues/212)


The Benchmark
=============

## Usage 


## 💾 Download
The current release (v1.0) of the benchmark can be found
[here](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/TandemRepeats_v1.0/).

## 📜 README
Included in the above tar-ball is a `README.md` which can be seen
[here](https://github.com/ACEnglish/adotto/tree/main/benchmark/GIABTR_benchmark.6.26/README.md).
This README has details on the files contained in the benchmark as well as instructions on how to compare your
caller using Truvari. Note that Truvari is currently under active development as a part of the GIAB TR project 
so a manual installation of the latest development branch is recommended. The README explains this in detail.


Other Data
==========
Here the confusion matrix

| Dataset | Positive | Negative |
| ------- | ------ | --------------- | ------- |
| Positive | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13987414.svg)](https://doi.org/10.5281/zenodo.13987414) | v1.2.1 | 
| Negative | [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6975244.svg)](https://doi.org/10.5281/zenodo.6975244) | v0.1 | 

This repository
===============
There are a few main sub-parts to this project. Each is contained in a sub-directory. Raw data that's too large to be
kept on github will be made available and documented such that a user can find it and know where to place it within
a clone of this repo in order to run sub-parts of the analysis.

* slides - GIAB team meeting slides
* manuscript - Summary and plotting workflows for the publication
* metadata - Data descriptor files (e.g. download paths of inputs used or sample ancestry information)
* regions - Identification of Tandem-Repeat regions of a reference
* variants - Calling variants from long-read haplotype resolved assemblies
* benchmark - Scripts for consolidating regions and variants to create the benchmark
