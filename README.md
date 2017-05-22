# hic-pipeline
This is the Hic mapping pipeline for Ren lab. The script is still at the very preliminary stage, please let us know how to improve!

# Dependencies 
## BWA (version that contains `bwa mem` submodule.)
## samtools 
## picard [Dedup bam files]
## juicebox.jar [to convert from bam to juicer format]
## GNU coreutils sort (v8.22 or newer; contains the parallel option) 

# Installation
See INSTALL

# How to run it.
1. Remember to install before running anything (see INSTALL)
2. Copy the run_hic.configure file into your own folder, 
and edit it to your need. 
3. (path-to-the-directory)/bin/run_hic.sh -c run_hic.configure

# Modules
## preprocessing (Fastq -> bam)
* Map reads from fastq from each read end.
* Combine reads from read pair; sort, merge.
* PCR duplicate removal (is it faster to do it here or in a later step?)
## bam to matrix 
* bam to juicer format
* bam to matrix. 
## calculate directional index and call TADs.
## calculate insulation scores. 
## calculate AB compartments 


