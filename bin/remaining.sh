
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

