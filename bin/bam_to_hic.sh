#! /usr/bin/env bash
BIN=$(dirname $0)


while getopts ":n:c:g:h:p:s:" OPT
do
    case $OPT in
    n) ID=$OPTARG;;
    c) SITE_POS=$OPTARG;;
    g) REF=$OPTARG;;
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

$samtools view $ID.final.bam | $BIN/sam2hic.pl | awk '{if ($3<$7) print $3$7,$0; else print $7$3,$0}' | $BIN/sort -T. --parallel=$NTHREADS -S 10% | cut -d ' ' -f 2- | gzip > $ID.juicer.txt.gz
java -jar $BIN/Juicebox.jar pre -q 30 -f $SITE_POS $ID.juicer.txt.gz $ID.hic $REF

echo -en "$(date) Leaving $(basename $0);\t"
echo $(date -u -d @"$SECONDS" +'%-Hh %-Mm %-Ss') elapsed

