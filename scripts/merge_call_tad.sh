# merge into one file
cat *.bed |  awk '!/NA|-Inf/' - | tr ' ' '\t' > ../SP108_insulation.bedgraph
cat *.norm.DI |  awk '!/NA|-Inf/' - | tr ' ' '\t' > ../D80_1_norm.DI.bedgraph
# call TAD from DI.
/mnt/silencer2/HiC/bin/DI2TAD.sh png/norm.DI hg19 40000 >& png/norm.log
