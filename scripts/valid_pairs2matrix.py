import sys
import gzip

inFileName = sys.argv[1]
outFilePrefix = sys.argv[2]
res = int(sys.argv[3])


def main():
  chr_dict = {}
  inFile = gzip.open(inFileName,'rt')
  for line in inFile:
    items = line.strip().split("\t")
    # trans reads
    if items[2] != items[6]:
      continue
    else:
      chr = items[2][3:]
      pos1 = int(int(items[3])/res)
      pos2 = int(int(items[7])/res)
      if pos1 > pos2:
        pos1,pos2 = pos2,pos1
      if chr in chr_dict:
        if (pos1,pos2) in chr_dict[chr]:
          chr_dict[chr][(pos1,pos2)] += 1
        else:
          chr_dict[chr][(pos1,pos2)] = 1
      else:
        chr_dict[chr] = {}
        chr_dict[chr][(pos1,pos2)] = 1
  # write reads.
  inFile.close()
  for chr in chr_dict:
    outFile = open(outFilePrefix + "_chr" + chr + ".txt",'w')
    for pos1,pos2 in chr_dict[chr]:
      outFile.write("\t".join([str(pos1),str(pos2),str(chr_dict[chr][(pos1,pos2)])])+"\n")
    outFile.close()


if __name__ == "__main__":
  main()
