#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <cassert>
#define cimg_display 0
#include "CImg.h"

using namespace std;
using namespace cimg_library;

int main(int argc, char*argv[]) {
  int bin=cimg_option("-bin", 40000, "bin_size");
  if (bin<=0) bin=40000;
  int cut=cimg_option("-cut", 10, "cut");
  if (cut<=0) cut=10;
  
  const char* format=cimg_option("-format", "sam", "sam or bed");
  bool bed=true;
  if (strcmp(format, "sam")==0) {
    bed=false;
  }
  const char* location=cimg_option("-pos", "chr1", "selection range, ex: chr1, chr2:1000000-8000000");
  const char* output=cimg_option("-o", "output.png", "output png");
  const char* xyz=cimg_option("-xyz", "", "output in 'chr\\tX\\tY\\tN' format");
  const char* mat_asc=cimg_option("-mat", "", "count matrix");
  const char* genome_size=cimg_option("-fai", "", "chrome size file");
  
  if (cimg_option("-h", false, "help")) {
    cerr << "location:\t" << location << endl
      << "output:\t" << output << endl
      << "xyz:\t" << xyz << endl
      << "mat:\t" << mat_asc << endl
      << "format:\t" << format << endl
      << "bin:\t" << bin << endl
      << "cut:\t" << cut << endl
      << "fai:\t" << genome_size << endl;
    return 0;
  }

  string ln;
  char buf[256]; int chr_start, chr_end;
  string chr=location; chr_start=chr_end=0;

  if (sscanf(location, "%127[^:]:%d-%d", buf, &chr_start, &chr_end)==3) {
    chr=buf; 
  } else {
    ifstream ifs(genome_size);
    int len; 
    while (getline(ifs, ln)) {
      if (sscanf(ln.c_str(), "%127s\t%d", buf, &len)==2) {
        if (chr.compare(buf)==0) {
          chr_end=len;
          break;
        }
      }
    }
  }

  if (!chr_end) {
    cerr << chr << " length not found, please privode chrome size files" << endl;
    return -2;
  }

  int chr_size = chr_end - chr_start;
 
  int w=chr_size/bin+1; CImg<int> mat(w,w);
  cerr << chr << ":" << chr_start << "-" << chr_end << " (" << chr_size << " bp) bin_size: " << bin;
  cerr << " matrix: " << w << "^2 ";

  int pairs=0, off1, off2;
  if (bed) {
    while (getline(cin, ln)) {
      if (sscanf(ln.c_str(), "%255s\t%d\t%d", buf, &off1, &off2)!=3) {
        cerr << "sscanf: " << ln << endl;
        continue;
      }
      int x=(off1-chr_start)/bin, y=(off2-chr_start)/bin;
      if (strcmp(buf, chr.c_str())!=0) {
        cerr << buf << " != " << chr << ": " << ln << endl;
        continue;
      }
      if (x<0 || y<0 || x>=w || y>=w) {
        cerr << x << "," << y << " out of range: " << ln << endl;
        continue;
      }
      pairs++;
      mat(x,y)++; mat(y,x)++;
    }
  } else { // sam
// samtools view -h -f 64 $OUTPUT_DIR/$ID.Hind3.bam $pos | awk '{print $3"\t"$4"\t"$8}' | $BIN_DIR/sam2mat -pos $pos -o $png -bin 100000 -mat $mat
// HWI-ST1113:452:HB02EADXX:1:1101:17111:27475	145	chr10	154201	60	100M	=	87595776	87441575	ATGTTGTTCCATCATCAATGTATTTGAAATATTGAAATGACAAAACTCAATGATTTTATGTCTCTTTGCAATGTAAGACATAATCTCTTTATAGGACGAA	DDDEEFFFFFFHGGHHHHIJJIJJJJJJIIIIIIJIIIIIHIIIHJJJJJJJJJIIHIJJJJJJJJJIJJJJIJJJJHJJJJJIJJJHHHHHFFDDFCCB	MD:Z:100	PG:Z:MarkDuplicates	NM:i:0	AS:i:100	XS:i:19
   while (getline(cin, ln)) {
     if (ln[0]=='@')
       continue;
#ifdef DEBUG
     stringstream ss(ln);
     string id, chr1, qual, cigar, eq;
     int flag;
     ss >> id >> flag >> chr1 >> off1 >> qual >> cigar >> eq >> off2;
     if (ss.fail()) {
       cerr << "Can not parse: " << ln << endl;
       return -3;
     }
#else
     char id[256], chr1[256], cigar[256], eq[256];
     int flag, off1, qual, off2;
     // cerr << ln << endl;
     int scan=sscanf(ln.c_str(), "%255s\t%d\t%s\t%d\t%d\t%s\t%s\t%d",
       id, &flag, chr1, &off1, &qual, cigar, eq, &off2);
     if (scan!=8) {
       cerr << "Can not parse: " << ln << endl;
       return -3;
     }
#endif
     if (chr.compare(chr1)!=0) {
       cerr << chr1 << " != " << chr << ": " << ln << endl;
       continue;
     }

     if (flag&64) {
#ifdef DEBUG
       cerr << "flag 64 warn: " << ln << endl;
#endif
       continue;
     }
     if (eq[0]!='=') {
#ifdef DEBUG
       cerr << eq << " eq? " << ln << endl;
#endif
       continue;
     }
    int x=(off1-chr_start)/bin, y=(off2-chr_start)/bin;
    if (x<0 || y<0 || x>=w || y>=w) {
#ifdef DEBUG
      cerr << x << "," << y << " out of range: " << ln << endl;
#endif
      continue;
    }
    pairs++;
    mat(x,y)++; mat(y,x)++;
   } 
  }

  if (strlen(mat_asc))
    mat.save(mat_asc);
  cerr << "read " << pairs << " records\n";
  mat.cut(0,cut).normalize(0,255)*(-1)+255;
  
  CImg<unsigned char> red(mat.width(), mat.height(), 1, 3, 255);
  cimg_forXY(mat, x, y) {
    unsigned char c=mat(x,y);
    red(x, y, 1)=255-c;
    red(x, y, 2)=255-c;
  }
  if (strlen(output))
    red.save(output);

  if (strlen(xyz)) {
    FILE *fp=fopen(xyz, "w");
    assert(fp);
    cimg_forXY(mat, x, y) {
	  if (x<y)
		continue;
	  if (mat(x,y))
        fprintf(fp, "%s\t%d\t%d\t%d\n", chr.c_str(), x, y, mat(x, y));			
    }	
    fclose(fp);
  }

  return 0;
}
