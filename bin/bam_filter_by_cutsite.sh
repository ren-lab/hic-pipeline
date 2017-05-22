#! /usr/bin/bash
# bam_filter_by_cutsite.sh
################
# This script will filter reads that are >500 of a enzyme cutting site.

set -o pipefail

while getopts ":n:c:h:p:s:" OPT
do
    case $OPT in
    n) ID=$OPTARG;;
    c) CUT_SITE=$OPTARG;;
    p) NTHREADS=$OPTARG;;
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


echo "$(date) Entering $(basename $0)"
SECONDS=0

bedtools intersect -abam $ID.dedup.bam -b <(zcat ${CUT_SITE}) | $samtools sort -@ $NTHREADS -n - | $samtools view -h - | awk '/^@/{print;next}{if($1==id){print ln"\n"$0}id=$1;ln=$0}' | $samtools view -Sb - | $samtools sort -@ $NTHREADS - > $ID.final.bam
if [ $? -ne 0 ]; then exit 1; fi
$samtools index $ID.final.bam
$samtools flagstat $ID.final.bam > qc/$ID.final.flagstat
grep -H ^[1-9] qc/$ID.final.flagstat

echo -en "$(date) Leaving $(basename $0);\t"
echo $(date -u -d @"$SECONDS" +'%-Hh %-Mm %-Ss') elapsed

