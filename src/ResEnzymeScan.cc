#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <cstring>

using namespace std;
size_t hits=0;

void scan(string& chr, string& seq, char* enzyme, int ext) {
  if (!chr.length())
    return;
  size_t pos=0, hit=0, ext_left=ext, ext_right=ext+strlen(enzyme);
  for (pos=0; (pos=seq.find(enzyme, pos))!=string::npos; pos++, hit++) {
    int left=pos-ext_left, right=pos+ext_right; if (left<0) left=0;
    printf("%s\t%d\t%d\n", chr.c_str(), left, right);
    hits++;
  }
  cerr << chr << "\t" << seq.length() << "\t" << hit << " hits" << endl;
}

int main(int argc, char* argv[]) {
  if (argc!=3) {
    cout << "Usage: " << argv[0] << " <enzyme seq>(AAGCTT|GATC) <ext_len>(500)\n";
    return -1;
  }

  char* enzyme=argv[1];
  int ext=atoi(argv[2]); 

  string ln, chr, seq;
  while (cin >> ln) {
    if (ln[0]=='>') { scan(chr, seq, enzyme, ext); chr=ln.substr(1); seq.clear(); continue; }
    transform(ln.begin(), ln.end(), ln.begin(), ::toupper);
    seq.append(ln);  
  }
  scan(chr, seq, enzyme, ext);
  cerr << hits << " total hits\n";

  return 0;
}
