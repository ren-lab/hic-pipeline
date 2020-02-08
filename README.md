# hic-pipeline
This is the Hi-C mapping pipeline for Ren lab. The script is still in the development stage, please let us know of any bugs and how to improve!

# Dependencies 
* BWA (v0.7.0, but preferably newer; contains the 'bwa mem' module)
* samtools 
* juicebox.jar (to convert from bam to juicer format)
* straw (to extract contact matrix from juicebox )
* GNU coreutils sort (v8.22 or newer; contains the parallel option) 
* snakemake (workflow management)

# Installation
See INSTALL

# How to run it.
1. Remember to install before running anything (see INSTALL).
2. Create your project folder. 
3. Inside your project folder, create a folder named `fastq`, and put all your .fastq files in that folder.  
4. Copy the snakemake.config.yaml file into your project folder, and edit it to your need. 
5. Run `(path-to-this-directory)/bin/run_hic_vanilla.sh -c snakemake.config.yaml -e your.email@ucsd.edu -s server`. 

6. Run `(path-to-this-directory)/bin/run_hic_options.sh -c snakemake.config.yaml -s server -o [options]` if you only need to generate part of the output, for instance, valid read pairs.

# Modules
## preprocessing (Fastq -> read pairs)
* Map reads from fastq from each read end
* Combine reads from read pair; sort, merge
* PCR duplicate removal
## read pairs to juicer and matrix
* pairs.txt to juicer format
* juicer format to matrix
## calculate directional index and call TADs
* calling TADs from DI requires a **MATLAB license**
## calculate insulation scores
  Use the method from Nagano et al. 2017; Olivares-Chauvet et al 2016
## calculate AB compartments
  Implemented the method from Rao et al. 2009. 

# File formats: 
## Read Pairs: 
```
Columns: 
=======
1)Read Name 
2) R1 mapping Flag
3) R1 chr 
4) R1 position 
5) R1 fragment 
6) R2 mapping Flag
7) R2 chr 
8) R2 position 
9) R2 fragment 
10) R1 mapping quality
11) R2 mapping quality
```

# Contributing Authors
* [Shawn Yanxiao Zhang](https://github.com/shawnzhangyx)
* Initial pieces from [Bin Li](https://github.com/bil022)
