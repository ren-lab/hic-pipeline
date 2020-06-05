#!/bin/bash
set -e

di=$1
ref=$2
bin=$3

if ! [ -e "$di" ]; then
  echo "Usage $0 <DI> <ref> <bin>"
  exit
fi

BIN=`readlink -f $0`
BIN=`dirname $BIN`
TAD_BASE="$BIN/domaincall_software"

fai=$ref.fa.fai
if ! [ -e "$fai" ]; then
  echo $fai not found
  exit
fi

TAD_PERL=$TAD_BASE/perl_scripts

export TAD_PATH=$TAD_BASE
export TAD_INPUT=$di
export TAD_OUTPUT=$di.m_out
#export TAD_FINAL-${di%%DI}TAD
#if ! [ -e ${di%%DI}TAD ]; then

#if ! [ -e $TAD_OUTPUT ]; then
matlab -nodisplay -nosplash -nodesktop < $BIN/TAD.m
#fi

perl $TAD_PERL/file_ends_cleaner.pl $TAD_OUTPUT $TAD_INPUT | perl $TAD_PERL/converter_7col.pl | sed 's/chrchr/chr/' > $di.$$.7col

min=2; prob=0.99; binsize=$bin
#for chr in `awk '{print $1}' $fai`; do
for chr in `awk '{print $1}' $di.$$.7col|sort -u`; do
  #echo "perl $TAD_PERL/hmm_probablity_correcter.pl <(awk '/^'$chr'\t/' $di.$$.7col) $min $prob $binsize | perl $TAD_PERL/hmm-state_caller.pl $fai $chr | perl $TAD_PERL/hmm-state_domains.pl"
  perl $TAD_PERL/hmm_probablity_correcter.pl <(awk '/^'$chr'\t/' $di.$$.7col) $min $prob $binsize | perl $TAD_PERL/hmm-state_caller.pl $fai $chr | perl $TAD_PERL/hmm-state_domains.pl
done > $di.TAD.bed
#else
#  echo skip ... ${di%%DI}TAD exists
#fi
