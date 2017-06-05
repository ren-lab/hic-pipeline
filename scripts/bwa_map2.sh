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

while getopts ":1:2:e:g:p:s:h:t:n:" OPT
do
    case $OPT in
    1) R1=$OPTARG;;
    2) R2=$OPTARG;;
    g) REF=$OPTARG;;
    p) NTHREADS=$OPTARG;;
    s) SAMTOOLS=$OPTARG;;
    t) N_LINES=$OPTARG;;
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
#echo N_LINES: $N_LINES

if ! [ -e $R1 ]; then echo "$R1 not found"; exit 1; fi
if ! [ -e $R2 ]; then echo "$R2 not found"; exit 1; fi


FIFO1=r1.fifo
FIFO2=r2.fifo
SAM=$ID.raw.sam
NTHREADS=$((NTHREADS/2))

mkfifo $FIFO1 $FIFO2 
### N_LINES is simply for testing purposes
if [ ! -z ${N_LINES+x} ]; then 
  echo N_LINES: $N_LINES
  bwa mem -t $NTHREADS $REF.fa <(bzcat $R1|head -n $((4*$N_LINES)) ) \
  2> log/bwa.R1.log |
  $BIN/chimeric.pl 2> log/bwa.F1.log > $FIFO1 |\
  bwa mem -t $NTHREADS $REF.fa <(bzcat $R2|head -n $((4*$N_LINES)) ) \
  2> log/bwa.R2.log |\
  $BIN/chimeric.pl 2> log/bwa.F2.log > $FIFO2 |\
  $BIN/mkPE2 -o $SAM -a $ID.valid_pairs.txt -p $SITE_POS -fnt -1 $FIFO1 -2 $FIFO2 -c 1 
else
  bwa mem -t $NTHREADS $REF.fa <(bzcat $R1 ) 2> log/bwa.R1.log |\
  $BIN/chimeric.pl 2> log/bwa.F1.log > $FIFO1 |\
  bwa mem -t $NTHREADS $REF.fa <(bzcat $R2 ) 2> log/bwa.R2.log |\
  $BIN/chimeric.pl 2> log/bwa.F2.log > $FIFO2 |\
  $BIN/mkPE2 -o $SAM -a $ID.valid_pairs.txt -p $SITE_POS -fnt -1 $FIFO1 -2 $FIFO2 -c 1 
fi

rm -f $FIFO1 $FIFO2 
echo -en $(date) "Leaving $(basename $0);\t"
echo $(date -u -d @"$SECONDS" +'%-Hh %-Mm %-Ss') elapsed
