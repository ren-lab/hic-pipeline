### convert bam file to valid pairs format. 
###
import sys
import pysam
import gzip 
import time

bamName  = sys.argv[1]
fragName = sys.argv[2]
outName = sys.argv[3]

# get the fragment position from chromosome and position.
def getFragmentPosition(frag_dict,chr,pos):
  frag_list = frag_dict[chr]
  #does not consider the boundary situation yet. 
#  for idx,num in enumerate(frag_list):
#    if num > pos:
    # test which is the closest element.
#      if num-pos > pos - frag_list[idx-1]:
#        return idx 
#      else:
#        return idx+1  
# binary look 
  l = 0
  r = len(frag_list)
  while l!=r-1:
    mid = int( (l+r)/2)
    if pos == frag_list[mid]:
      return mid + 1
    elif pos > frag_list[mid]:
      l=mid
    else: 
      r=mid
  return l + 1
  
def main():
  #read the fragment information. 
  frag_dict = {}
  print("Read the fragment dictionary")
  with open(fragName,'r') as f:
    for line in f: 
      items = line.strip().split(" ")
      frag_dict[items[0]] = [ int(x) for x in items[1:]]
      print(items[0])
  print(len(frag_dict))

  # read sam file. 
  sam = pysam.AlignmentFile(bamName,'rb')
  # open output txt file.
  outF = gzip.open(outName,'wt')  

  tic = time.time()
  for e,line in enumerate(sam.fetch(until_eof=True)):
    if e % 2 ==0:
      readname = line.query_name
      chr1 = line.reference_name
      pos1 = line.pos
      flag1 = line.flag
      mapq1 = line.mapping_quality
      frag1 = getFragmentPosition(frag_dict,chr1,pos1)
    #  print(pos1-frag_dict[chr1][frag1-1],frag_dict[chr1][frag1]-pos1)
    else:
      chr2 = line.reference_name
      pos2 = line.pos
      flag2 = line.flag
      mapq2 = line.mapping_quality
      frag2 = getFragmentPosition(frag_dict,chr2,pos2)
      outF.write("\t".join([readname,str(flag1),chr1,str(pos1),str(frag1),
      str(flag2),chr2,str(pos2),str(frag2),str(mapq1),str(mapq2)])+"\n")
#      outF.write(bytes("\t".join([readname,str(flag1),chr1,str(pos1),str(frag1),
#      str(flag2),chr2,str(pos2),str(frag2),str(mapq1),str(mapq2)])+"\n",encoding="utf-8"))
    #if e % 100000==1: 
      # print(str(e)+" elapsed time:"+ str(time.time()-tic))
      tic=time.time()
#    if e == 50000: break

  outF.close()
if __name__ == "__main__":
  main()
