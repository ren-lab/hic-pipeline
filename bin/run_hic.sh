#!/bin/bash
### Copy left (c) 2017 Bing Ren Lab
### GNU GPL v3 License. 
###################################


# Reading the arguments to the program 
# if the number of arguments is less than 1 
# print the help document. 
echo $#
if [ $# != 1 ]; then
  echo -e "#Usage:\n\t $0 <HiC conf file>\n"
  echo -e "#HiC conf Example:\n"
  echo -e 'ID_LIST=(Sample1 Sample2)
R1=fastq/${ID}.R1.fastq.bz2
R2=fastq/${ID}.R2.fastq.bz2
ZCAT=bzcat
REF=hg19
NTHREADS=8
EMAIL=shz254@ucsd.edu'
  exit
fi

# load configuration file. 
conf=$1
if ! [ -e "$conf" ]; then echo "$conf not found"; exit; fi
echo "Configuration file: $conf, samtools: $samtools"
. $conf

# root directory
PWD=`pwd`
# bin directory
BIN=`readlink -f $0`
BIN=`dirname $BIN`
# genome index directory
REF_DIR=$BIN/../genome_index/bwa_indices
OUT_DIR=$PWD/$OUT_DIR
R1=$PWD/$R1
R2=$PWD/$R2

echo $PWD
echo $BIN
echo $OUT_DIR
echo $R1 $R2

# make directories 
  [ -e $OUT_DIR/log ] || mkdir $OUT_DIR/log


# test alignment
bwa mem -t $NTHREADS $REF_DIR/${REF}.fa $R1 2> $OUT_DIR/log/${ID}.bwa.R1.log > $OUT_DIR/${ID}.R1.raw.sam 
bwa mem -t $NTHREADS $REF_DIR/${REF}.fa $R2 2> $OUT_DIR/log/${ID}.bwa.R1.log > $OUT_DIR/${ID}.R2.raw.sam


echo success
