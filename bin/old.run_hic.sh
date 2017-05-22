#!/bin/bash
#######################################################################
### Copyleft (c) 2017 Bing Ren Lab
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#######################################################################


# Reading the arguments to the program 
# if the number of arguments is less than 1 
# print the help document. 
if [ $# -ne 1 ]; then
  echo -e "#Usage:\n\t $0 <HiC conf file>\n"
  exit 1
fi

# load configuration file
conf=$1
if ! [ -e "$conf" ]; then echo "$conf not found"; exit; fi
echo "Configuration file: $conf, samtools: $samtools"

PWD=`pwd`
. $conf
for ID in ${ID_LIST[@]}; do
  ID=`basename $ID`
  . $conf
  BIN=`readlink -f $0`
  BIN=`dirname $BIN`
  if ! [ -e "$samtools" ]; then samtools=$BIN/samtools; fi
  if [ -z "$TAD_WIN" ]; then TAD_WIN=50; fi
  REF_DIR=$BIN/../annotation
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
  find $REF_DIR/bwa_indices/$REF.fa $REF_DIR/bwa_indices/$REF.fa.fai> /dev/null
  if [ $? -ne 0 ]; then exit; fi
  
  [ -e log ] || mkdir log
  if ! [ -e log ]; then echo can not create log dir; exit; fi
 
  FIFO1=$ID.R1.fifo
  FIFO2=$ID.R2.fifo
  SAM=$ID.sam
  HIC=$ID.hic
  SITE_POS=$BIN/../annotation/juicebox/site_pos/${REF}_$RES_ENZYME.txt
  if ! [ -e "$SITE_POS" ]; then echo $SITE_POS not found; exit; fi
  if ! [ -e $ID.bwa.bam ]; then
    echo -n "Begin $ID BWA: "; date
    if ! [ -e log/$ID.raw.flagstat ]; then
      rm -f $FIFO1 $FIFO2 $SAM $HIC.txt

      mkfifo $FIFO1 $FIFO2 $SAM
      echo $REF_DIR/bwa_indices/$REF.fa
      bwa mem -t $NTHREADS $REF_DIR/bwa_indices/$REF.fa <($ZCAT $R1) 2> log/$ID.bwa.R1.log | $FILT 2> log/$ID.F1.log > $FIFO1 &
      bwa mem -t $NTHREADS $REF_DIR/bwa_indices/$REF.fa <($ZCAT $R2) 2> log/$ID.bwa.R2.log | $FILT 2> log/$ID.F2.log > $FIFO2 &
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
  
  if ! [ -e $REF_DIR/cut_sites/$REF.$RES_ENZYME.bed.gz ]; then
    cat $REF_DIR/bwa_indices/$REF.fa | $BIN/ResEnzymeScan $RES_ENZYME 500 | bedtools merge -i - | gzip > $REF_DIR/cut_sites/$REF.$RES_ENZYME.bed.gz
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
 
done

if ! [ -z "$EMAIL" ]; then
  echo $conf done | mail -s $conf $EMAIL
else
  echo $conf done
fi

