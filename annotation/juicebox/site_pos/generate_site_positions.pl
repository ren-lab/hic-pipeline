#!/usr/bin/env perl

($res_enz, $ref)=@ARGV;
die "Usage $0 res_enz ref" unless $ref;

$res_enz="AAGCTT" if $res_enz=~/HindIII/i;
$res_enz="GATC" if $res_enz=~/DpnII/i;
$res_enz="GATC" if $res_enz=~/MboI/i;
$off=length($res_enz)-1;
$res_enz=uc $res_enz;

die $res_enz unless $res_enz=~/^[ACTG]+$/;

die "$ref not found" unless ( -e "$ref" );

open REF, "<$ref" or die "$ref";
while (<REF>) { chomp();
  if (/^>(\S+)/) {
    if ($chr) { search(); $seq=""; }
    $chr=$1; next;
  }
  $seq.=uc $_;
}
if ($chr) { search(); }

sub search {
  print "$chr";
  #print STDERR "$res_enz ".length($seq)."\n";
  while ($seq=~/$res_enz/g) {
    print " ".(pos($seq)-$off); 
  }
  print " ".length($seq)."\n";
}
