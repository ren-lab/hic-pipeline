#!/usr/bin/perl
use strict;
use constant DEBUG => 0;

sub Log {
  print STDERR @_ if DEBUG;
}

sub MDist {
  my ($flag, $cigar)=@_;
  Log "MDIST: ".$flag."\t".$cigar."\n";
  my $off=-1;
  my $len=0;
  die unless $cigar=~/\d+M/; 
  if ($flag&16) { # rev
    $off=0;
    while ($cigar=~/(\d+)([A-Z])/g) {
      my ($n, $X)=($1,$2);
      $len+=$n;
      $off=$len if $X=~/M/;
      Log "off-: -".$off."\t".$X."\n";
    }
    Log "off-: ret ".($len-$off)."\n";
    return $len-$off; 
  }
  while ($cigar=~/(\d+)([A-Z])/g) {
    my ($n, $X)=($1,$2);
    $off=$len if ($off<0 && $X=~/M/);
    Log "off+: ".$off."\t".$X."\n";
    $len+=$n;
  }
  Log "off+: ret ".$off."\n";
  return $off;
}

sub FivePrime {
  my ($a, $b)=@_;
  my @As=split(/\t/, $a);
  my @Bs=split(/\t/, $b);
  my $distA=MDist($As[1],$As[5]);
  my $distB=MDist($Bs[1],$Bs[5]);
  Log "A/B: ".$distA." ".$distB."\n";
  if ($distA<$distB) {
    return $a;
  }
  return $b;
}

my ($prev_id, $prev_ln);
while (<>) { 
  if (/^@/) { print; next; }
  chomp();
  my ($id, $flag, $chr, $pos, $mapq, $cigar, $d1, $d2, $d3, $seq, $qual, @rest) = split /\t/;
  if ($prev_id ne $id) { 
    print "$prev_ln\n" if $prev_ln;
    ($prev_ln, $prev_id)=($_, $id);
    next;
   }
   Log "$prev_ln\n$_\n";
   $prev_ln=FivePrime($prev_ln, $_);
   Log "=> $prev_ln\n\n";
}
print "$prev_ln\n" if $prev_ln;
