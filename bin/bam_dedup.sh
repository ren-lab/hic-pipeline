#! /usr/bin/bash
################
# This script will dedup the bam file. 

BIN=$(dirname $0)
MARK_DUP="java -jar $BIN/picard.jar MarkDuplicates"

function usage(){
echo -e "Usage: $0" 
}

## processing the command options
while getopts ":n:s:h:" OPT
do
    case $OPT in
    n) ID=$OPTARG;;
    s) samtools=$OPTARG;;
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

## beginning the scripts
echo "$(date) Entering $0"

$MARK_DUP INPUT=$ID.raw.bam  OUTPUT=$ID.dedup.bam ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=. METRICS_FILE=log/$ID.metrics.log &> log/$ID.markdup.log
grep ^LIBRARY -A 1 log/$ID.metrics.log
$samtools index $ID.dedup.bam
$samtools flagstat $ID.dedup.bam > qc/$ID.dedup.flagstat
grep -H ^[1-9] qc/$ID.dedup.flagstat

echo "$(date) Leaving $0"


