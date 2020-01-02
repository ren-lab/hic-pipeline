for ref in panPan2  panTro6 calJac3 ; do
  for s in  GATC; do
#    if ! [ -e "${ref}_$s.txt" ]; then
      echo $ref $s
      ./generate_site_positions.pl $s /projects/ps-renlab/share/bwa_indices/$ref.fa > ${ref}_$s.txt 
      #./generate_site_positions.pl $s ~/annotations/bwa_index/$ref.fa | sed 's/^chrM/chrMT/' | sed 's/^chr//' > ${ref}_$s.txt
#    fi
  done
done
