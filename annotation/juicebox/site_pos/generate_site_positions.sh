for ref in mm9 mm10 hg18 hg19; do
  for s in HindIII DpnII AAGCTT GATC; do
    if ! [ -e "${ref}_$s.txt" ]; then
      echo $ref $s
      ./generate_site_positions.pl $s /mnt/raid1/IMPORT/Illumina/bowtie_indexes/$ref.fa | sed 's/^chrM/chrMT/' | sed 's/^chr//' > ${ref}_$s.txt
    fi
  done
done
