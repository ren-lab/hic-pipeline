#! /usr/bin/env bash
# bwa_map2.sh
################
# This function will take a few arguments and map the reads using BWA MEM
################
set -e
set -o pipefail  #will return non-zero status for broken pipe. 
BIN=$(dirname $0)/../bin/

function usage(){
echo -e "Usage: $0" 
}

while getopts ":1:2:e:g:p:s:h:n:" OPT
do
    case $OPT in
    1) R1=$OPTARG;;
    2) R2=$OPTARG;;
    g) REF=$OPTARG;;
    p) NTHREADS=$OPTARG;;
    s) SAMTOOLS=$OPTARG;;
    e) SITE_POS=$OPTARG;;
    n) ID=$OPTARG;;
    h) help ;;
    \?)
         echo "Invalid option: -$OPTARG" >&2
         usage
         exit 1
         ;;
     :)
         echo "Option -$OPTARG requires an argument." >&2
         usage
         exit 1
         ;;
    esac
done


echo $(date) "Entering $(basename $0)"; SECONDS=0

#echo R1: $R1 
#echo R2: $R2
#echo REF: $REF 

if ! [ -e $R1 ]; then echo "$R1 not found"; exit 1; fi
if ! [ -e $R2 ]; then echo "$R2 not found"; exit 1; fi

if [ ${R1: -3} == ".gz" ]; then
  Read="zcat"
  elif [ ${R1: -4} == ".bz2" ]; then
  Read="bzcat"
  elif [ ${R1: -6} == ".fastq" ]; then
  Read="cat"
  else
    echo "File extension of $R1 not recognized.
          Only *.gz *.bz2 or *.fastq are allowed"
    exit 1
fi



FIFO1=$ID.r1.fifo
FIFO2=$ID.r2.fifo
SAM=$ID.raw.sam
NTHREADS=$((NTHREADS/2))

if [ ! -e $FIFO1 ]; then mkfifo $FIFO1; fi
if [ ! -e $FIFO2 ]; then mkfifo $FIFO2; fi

bwa mem -t $NTHREADS $REF.fa <($Read $R1 ) 2> log/$ID.bwa.R1.log |\
  $BIN/chimeric.pl 2> log/$ID.bwa.F1.log > $FIFO1 |\
bwa mem -t $NTHREADS $REF.fa <($Read $R2 ) 2> log/$ID.bwa.R2.log |\
  $BIN/chimeric.pl 2> log/$ID.bwa.F2.log > $FIFO2 |\
  $BIN/mkPE3 -o $SAM -a $ID.valid_pairs.txt -p $SITE_POS -fnt -1 $FIFO1 -2 $FIFO2 -c 1 

rm -f $FIFO1 $FIFO2 
echo -en $(date) "Leaving $(basename $0);\t"
echo $(date -u -d @"$SECONDS" +'%-Hh %-Mm %-Ss') elapsed
