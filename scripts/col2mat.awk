#!/bin/awk -f

BEGIN {
OFS=" ";
#chr="chr19";
#bin_size=1000000;
bin_num=int(chr_size/bin_size)+1;
a[bin_num,bin_num]=0;
#print chr, bin_size, chr_size, bin_num;
print bin_num, bin_num,1,1;
}

{
  if ($3 != "nan"){
  a[int($1/bin_size),int($2/bin_size)] = $3;
  a[int($2/bin_size),int($1/bin_size)] = $3;
  }
  else {
  a[int($1/bin_size),int($2/bin_size)] = 0;
  a[int($2/bin_size),int($1/bin_size)] = 0;
  }

}
END {
  for (i=0;i<bin_num;i++){
    for (j=0;j<bin_num;j++){
      #print i,j;
      idx = i SUBSEP j;
      #print idx;
      printf ((idx in a)?a[idx]:0) OFS
      }
      print ""
      }
}



