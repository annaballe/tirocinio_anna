Report
==============

Step 1: Parsing VCF Files and Converting to BED

**Input:**

- A compressed VCF file with genomic data on short tandem repeats (STRs).

**Objective:**

- Extract key genomic coordinates and information to produce a BED file for downstream analysis.
  
The Benchmark
=============

## Usage 


## 💾 Download
The current release (v1.0) of the benchmark can be found
[here](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/TandemRepeats_v1.0/).

## 📜 README
To update


Other Data
==========
Here the confusion matrix

| Dataset | Positive | Negative |
| ------- | ------ | --------- | 
| Positive | 55 | 67 | 
| Negative | 12 | 47 | 

This repository
===============
Description of the files

* regions - Identification of Tandem-Repeat regions of a reference
* variants - Calling variants from long-read haplotype resolved assemblies
* benchmark - Scripts for consolidating regions and variants to create the benchmark
