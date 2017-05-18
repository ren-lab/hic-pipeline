#!/usr/bin/perl

while (<>) {
    next unless /HC:Z:FRG1=(\d+),FLG2=(\d+),FRG2=(\d+),MPQ2=(\d+)/;
    ($frag1, $flag2, $frag2, $mapq2)=($1, $2, $3, $4);
    ($rid, $flag, $chr, $pos, $mapq, $cigar, $chr2, $pos2, $dist, $seq, $qual, @tail)=split(/\t/, $_);
    $chr2=$chr if $chr2=~/=/;
    $str1=$flag&16; $str2=$flag2&16;
    # <readname> <str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2> <mapq1> <mapq2>
    print "$rid $str1 $chr $pos $frag1 $str2 $chr2 $pos2 $frag2 $mapq $mapq2\n";
}
