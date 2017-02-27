#!/bin/bash
  
if [ $# -eq 0 ]; then
  echo -e "#Usage:\n\t $0 <HiC conf file>\n"
  echo -e "#HiC conf Example:\n"
  echo -e 'ID_LIST=(Sample1 Sample2)
R1=fastq/${ID}.R1.fastq.bz2
R2=fastq/${ID}.R2.fastq.bz2
ZCAT=bzcat
REF=mm10
NTHREADS=4
PE_CUTOFF=15000
#RES_ENZYME=AAGCTT
RES_ENZYME=GATC
TAD_BIN=40000
TAD_WIN=50
EMAIL=bil022@ucsd.edu'
  exit
fi

conf=$1
if ! [ -e "$conf" ]; then echo "$conf not found"; exit; fi
echo "Configuration file: $conf, samtools: $samtools"

NCPU=`grep -c ^processor /proc/cpuinfo`
loadavg() { while [ `cat /proc/loadavg | awk '{print int($1)}'` -gt $NCPU ]; do sleep 300; done; }

PWD=`pwd`
. $conf
for ID in ${ID_LIST[@]}; do
  loadavg
  ID=`basename $ID`
  . $conf
  BIN=`readlink -f $0`
  BIN=`dirname $BIN`
  if ! [ -e "$samtools" ]; then samtools=$BIN/samtools; fi
  if [ -z "$TAD_WIN" ]; then TAD_WIN=50; fi
  REF_DIR=$BIN/../ref/
  FILT=$BIN/chimeric.pl
  MAKE_PE=$BIN/mkPE
  MARK_DUP="java -jar $BIN/picard.jar MarkDuplicates"

  R1=$PWD/$R1
  R2=$PWD/$R2
 
  which bwa $samtools bedtools java > /dev/null
  if [ $? -ne 0 ]; then exit; fi

  [ -e $ID ] || mkdir $ID
  if ! [ -e $ID ]; then echo can not create dir $ID; exit; fi

  pushd $ID
  find $FILT $MAKE_PE $R1 $R2 $BIN/ResEnzymeScan $BIN/picard.jar $BIN/sam2mat > /dev/null
  if [ $? -ne 0 ]; then exit; fi
  find $REF_DIR/$REF.fa $REF_DIR/$REF.fa.fai> /dev/null
  if [ $? -ne 0 ]; then exit; fi
  
  [ -e log ] || mkdir log
  if ! [ -e log ]; then echo can not create log dir; exit; fi
 
  FIFO1=$ID.R1.fifo
  FIFO2=$ID.R2.fifo
  SAM=$ID.sam
  HIC=$ID.hic
  SITE_POS=$BIN/juicebox/site_pos/${REF}_$RES_ENZYME.txt
  if ! [ -e "$SITE_POS" ]; then echo $SITE_POS not found; exit; fi
  if ! [ -e $ID.bwa.bam ]; then
    echo -n "Begin $ID BWA: "; date
    if ! [ -e log/$ID.raw.flagstat ]; then
      rm -f $FIFO1 $FIFO2 $SAM $HIC.txt

      mkfifo $FIFO1 $FIFO2 $SAM
      bwa mem -t $NTHREADS $REF_DIR/$REF.fa <($ZCAT $R1) 2> log/$ID.bwa.R1.log | $FILT 2> log/$ID.F1.log > $FIFO1 &
      bwa mem -t $NTHREADS $REF_DIR/$REF.fa <($ZCAT $R2) 2> log/$ID.bwa.R2.log | $FILT 2> log/$ID.F2.log > $FIFO2 &
      $MAKE_PE -o $SAM -p $SITE_POS -fnt -1 $FIFO1 -2 $FIFO2 -c $PE_CUTOFF &
      $samtools view -Sb $SAM | $samtools sort -@ $NTHREADS - > $ID.raw.bam
      rm -f $FIFO1 $FIFO2 $SAM

      #$MAKE_PE -f -1 $FIFO1 -2 $FIFO2 -c $PE_CUTOFF 2> log/$ID.mkPE.log | $samtools view -Sb - | $samtools sort -@ $NTHREADS - $ID.raw
      $samtools flagstat $ID.raw.bam > log/$ID.raw.flagstat
    fi
    grep -H ^[1-9] log/$ID.raw.flagstat
    echo -n "END $ID BWA: "; date; echo

    echo -n "Begin $ID RMDUP: "; date
    $MARK_DUP INPUT=$ID.raw.bam  OUTPUT=$ID.bwa.bam ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT TMP_DIR=. METRICS_FILE=log/$ID.metrics.log >& log/$ID.markdup.log
    grep ^LIBRARY -A 1 log/$ID.metrics.log
    $samtools index $ID.bwa.bam
    $samtools flagstat $ID.bwa.bam > log/$ID.bwa.flagstat
    grep -H ^[1-9] log/$ID.bwa.flagstat
    # rm -f $ID.raw.bam
    echo -n "END $ID RMDUP: "; date; echo
  fi
  
  if ! [ -e $REF_DIR/$REF.$RES_ENZYME.bed.gz ]; then
    cat $REF_DIR/$REF.fa | $BIN/ResEnzymeScan $RES_ENZYME 500 | bedtools merge -i - | gzip > $REF_DIR/$REF.$RES_ENZYME.bed.gz
  fi
  
  if ! [ -e $ID.$RES_ENZYME.bam ]; then
    echo -n "Begin $ID $RES_ENZYME: "; date
    bedtools intersect -abam $ID.bwa.bam -b <(zcat $REF_DIR/$REF.$RES_ENZYME.bed.gz) | $samtools sort -@ $NTHREADS -n - | $samtools view -h - | awk '/^@/{print;next}{if($1==id){print ln"\n"$0}id=$1;ln=$0}' | $samtools view -Sb - | $samtools sort -@ $NTHREADS - > $ID.$RES_ENZYME.bam
    $samtools index $ID.$RES_ENZYME.bam
    $samtools flagstat $ID.$RES_ENZYME.bam > log/$ID.$RES_ENZYME.flagstat
    grep -H ^[1-9] log/$ID.$RES_ENZYME.flagstat
    echo -n "End $ID $RES_ENZYME: "; date; echo
  fi

  if ! [ -s $HIC ]; then
    echo "$samtools view $ID.$RES_ENZYME.bam | $BIN/sam2hic.pl | awk '{if (\$3<\$7) print \$3\$7,\$0; else print \$7\$3,\$0}' | $BIN/sort -T. --parallel=$NTHREADS -S 10% | cut -d ' ' -f 2- | gzip > $HIC.txt.gz"
    echo "java -jar $BIN/Juicebox.jar pre -q 30 -f $SITE_POS $HIC.txt.gz $HIC $REF"
    $samtools view $ID.$RES_ENZYME.bam | $BIN/sam2hic.pl | awk '{if ($3<$7) print $3$7,$0; else print $7$3,$0}' | $BIN/sort -T. --parallel=$NTHREADS -S 10% | cut -d ' ' -f 2- | gzip > $HIC.txt.gz
    java -jar $BIN/Juicebox.jar pre -q 30 -f $SITE_POS $HIC.txt.gz $HIC $REF
    if [ $? -eq 0 ]; then
        rm $HIC.txt.gz
    fi
  fi
 
  [ -e png ] || mkdir png
  if ! [[ $TAD_BIN == +([0-9]) ]]; then TAD_BIN=40000; fi
  gnf=$REF_DIR/gen_features/$REF.$RES_ENZYME.$TAD_BIN.gnf
  if [ -e "$gnf" ]; then
    for chr in `awk '{print $1}' $REF_DIR/$REF.fa.fai`; do
      if ! [ -e png/$ID.$chr.png ]; then
        echo $chr
        # echo "$samtools view -h -F 64 $ID.$RES_ENZYME.bam $chr | $BIN/sam2mat -format sam -pos $chr -o png/$ID.$chr.png -bin $TAD_BIN -fai $REF_DIR/$REF.fa.fai"
        $samtools view -h -F 64 $ID.$RES_ENZYME.bam $chr | $BIN/sam2mat -format sam -pos $chr -o png/$ID.$chr.png -bin $TAD_BIN -fai $REF_DIR/$REF.fa.fai -cut 5 -mat png/$ID.$chr.asc
        $BIN/asc2di png/$ID.$chr.asc $chr $TAD_BIN $TAD_WIN $REF_DIR/$REF.fa.fai >& png/$ID.$chr.DI
        Rscript $BIN/asc2norm.R png/$ID.$chr.asc $gnf $chr > png/$ID.$chr.norm.asc
        $BIN/asc2di png/$ID.$chr.norm.asc $chr $TAD_BIN $TAD_WIN $REF_DIR/$REF.fa.fai >& png/$ID.$chr.norm.DI
      fi
    done
  else
    echo "genome feature $gnf not found"
  fi
  popd
  # time $samtools view -h -F 64 $ID.$RES_ENZYME.bam chr1 | awk '!/^@/{print $3"\t"$4"\t"$8}' | ../bin/sam2mat -format bed -pos chr1:0-249250621 -o $ID.chr1.bed.png -fai $REF_DIR/$REF.fa.fai
done

if ! [ -z "$EMAIL" ]; then
  echo $conf done | mail -s $conf $EMAIL
else
  echo $conf done
fi

#http://stackoverflow.com/questions/1537956/bash-limit-the-number-of-concurrent-jobs
#joblist=($(jobs -p))
#while (( ${#joblist[*]} >= 3 ))
#do
#    sleep 1
#    joblist=($(jobs -p))
#done
#$* &
