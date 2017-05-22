#! /usr/bin/bash
################
# This function will take a few arguments and map the reads using BWA MEM
################

BIN=$(dirname $0)

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


echo $(date) "Entering $(basename $0)" 
echo R1: $R1 
echo R2: $R2
echo REF: $REF 
echo N_LINES: $N_LINES

if ! [ -e $R1 ]; then echo "$R1 not found"; exit 1; fi
if ! [ -e $R2 ]; then echo "$R2 not found"; exit 1; fi


FIFO1=r1.fifo
FIFO2=r2.fifo
SAM=temp.sam

mkfifo $FIFO1 $FIFO2 $SAM
bwa mem -t $NTHREADS $REF.fa <(bzcat $R1 |head -n $((4*$N_LINES)) ) 2> log/bwa.R1.log | $BIN/chimeric.pl 2> log/bwa.F1.log > $FIFO1 &
bwa mem -t $NTHREADS $REF.fa <(bzcat $R2 |head -n $((4*$N_LINES)) ) 2> log/bwa.R2.log | $BIN/chimeric.pl 2> log/bwa.F2.log > $FIFO2 &
$BIN/mkPE -o $SAM -p $SITE_POS -fnt -1 $FIFO1 -2 $FIFO2 -c 1 &
$SAMTOOLS view -Sb $SAM | $SAMTOOLS sort -@ $NTHREADS - > $ID.raw.bam
rm -f $FIFO1 $FIFO2 $SAM
$SAMTOOLS flagstat $ID.raw.bam > qc/$ID.raw.flagstat 
grep -H ^[1-9] qc/$ID.raw.flagstat

echo $(date) "Leaving $(basename $0)"
