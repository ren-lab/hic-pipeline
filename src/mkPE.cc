/* mkPE.cc 
 * This is the original version from Bin Li. 
*/
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <algorithm>
#include <vector>
#include <map>
#include <getopt.h>

// https://broadinstitute.github.io/picard/explain-flags.html
#define  PAIRED 1
#define  UNMAPD 4
#define  REVRSE 16
#define  MATERV 32
#define  FIRSTP 64
#define  SECNDP 128
#define  NOTPRM 256
#define  SUPALN 2048

using namespace std;

void usage(char* prog) {
    cout << "Usage: " << prog << " -1 R1 -2 R2 -c cutoff -o sam -p site_pos -snft -j juicbox\n";
    cout << "\ts:\toutput single\n";
    cout << "\tn:\toutput near\n";
    cout << "\tf:\toutput far\n";
    cout << "\tt:\toutput trans\n";
    exit(-1);
}

map<string, vector<int> > sitePos;
map<string, map<int, int> > fragPos;

int getFrag(string& chr, int pos) {
    vector<int>& vec=sitePos[chr];
    vector<int>::iterator low=lower_bound(vec.begin(), vec.end(), pos);
    if (low==vec.end()) { cerr << chr << endl; }
    assert(low!=vec.end());
    int frag=fragPos[chr][*low];
    return frag;
}

int main(int argc, char* argv[]) {
    int c, cutoff=15000;
    char* pR1=NULL, *pR2=NULL, *pSitePos=NULL, *pSam=NULL;
    bool bSingle=false, bNear=false, bFar=false, bTrans=false;
    while ((c = getopt (argc, argv, "1:2:c:p:snfto:h")) != -1) {
        switch (c) {
            case '1': pR1=optarg; break;
            case '2': pR2=optarg; break;
            case 'c': cutoff=atoi(optarg); break;
            case 'p': pSitePos=optarg; break;
            case 's': bSingle=true; break;
            case 'n': bNear=true; break;
            case 'f': bFar=true; break;
            case 't': bTrans=true; break;
            case 'o': pSam=optarg; break;
            default: usage(argv[0]);
        }
    }

    assert(pR1); ifstream R1(pR1); assert(R1.good());
    assert(pR2); ifstream R2(pR2); assert(R2.good());
    assert(cutoff); assert(bFar);

    ofstream ofSam;
    ostream* outp = &cout;
    if (pSam) { ofSam.open(pSam); outp=&ofSam; }
    if (pSitePos) {
        cerr << "Loading site_pos...\n";
        ifstream site_pos(pSitePos); 
        assert(site_pos.good());
        string line, chr; int pos;
        while (getline(site_pos, line)) {
            stringstream ss(line);
            ss>>chr;
            vector<int>& vec=sitePos[chr];
            map<int, int>& m=fragPos[chr];
            vec.push_back(0); 
            int frag=0;
            m[0]=0; 
            while (ss>>pos) { frag++;
                vec.push_back(pos);
                m[pos]=frag;
            }
        }
        cerr << "Done site_pos\n";
    }

    int totalPE=0, mapR1=0, mapR2=0, mapPE=0, cisPE=0, cis15K=0, transPE=0;
    string line1, line2, rid1, rid2, chr1, chr2, tail1, tail2, cigar1, cigar2, rnext1, rnext2;
    int flag1, flag2, pos1, pos2, mapq1, mapq2, pnext1, pnext2, tlen1, tlen2;
    map<int, int> histK;
    cerr << "Begin mkPE" << endl;
    while (getline(R1, line1)) { assert(!R2.eof());
        getline(R2, line2);
        if (line1[0]=='@') { *outp << line1 << endl; continue; }
        totalPE++;

        if ((totalPE%1000000)==0)
            cerr << totalPE << endl;

        stringstream ss1(line1);
        ss1 >> rid1 >> flag1 >> chr1 >> pos1 >> mapq1 >> cigar1 >> rnext1 >> pnext1 >> tlen1; assert(!ss1.fail());
        getline(ss1, tail1);
        stringstream ss2(line2);
        ss2 >> rid2 >> flag2 >> chr2 >> pos2 >> mapq2 >> cigar2 >> rnext2 >> pnext2 >> tlen2; assert(!ss2.fail());
        assert(rid1.compare(rid2)==0);
        getline(ss2, tail2);	

        if (!(flag1&UNMAPD)) mapR1++;
        if (!(flag2&UNMAPD)) mapR2++;

        if ((flag1&UNMAPD)|(flag2&UNMAPD)) {
            if (!(flag1&UNMAPD) && bSingle)
                *outp << line1 << endl;
            if (!(flag2&UNMAPD) && bSingle)
                *outp << line2 << endl;
            continue;
        }
        mapPE++;

        int mate1=(flag2 & REVRSE)?MATERV:0;
        int mate2=(flag1 & REVRSE)?MATERV:0;

        flag1 |= ( PAIRED | mate1 | FIRSTP);
        flag2 |= ( PAIRED | mate2 | SECNDP);

        flag1 &= ~(NOTPRM | SUPALN);
        flag2 &= ~(NOTPRM | SUPALN);

        int frag1=0, frag2=0;
        if (pSitePos) {
            frag1=getFrag(chr1, pos1);
            frag2=getFrag(chr2, pos2);
        }

        if (chr1.compare(chr2)==0) {
            cisPE++;
            int dist=pos1-pos2; if (dist<0) dist=-dist;
            histK[int(dist/1000)]++;
            //$d1[6]="="; $d1[7]=$R2_pos; $d1[8]=$R2_pos-$R1_pos;
            //$d2[6]="="; $d2[7]=$R1_pos; $d2[8]=$R1_pos-$R2_pos;
            bool bPrint=false;
            if (dist>=cutoff && bFar) {
                cis15K++; bPrint=true;
            } else if (bNear) {
                bPrint=true;
            }
            if (bPrint) {
                *outp << rid1 << "\t" << flag1 << "\t" << chr1 << "\t" << pos1 << "\t" << mapq1 << "\t" << cigar1;
                *outp << "\t=\t" << pos2 << "\t" << (pos2-pos1) << tail1;
                if (pSitePos) {
                    *outp << "\tHC:Z:FRG1="<<frag1<<",FLG2="<<flag2<<",FRG2="<<frag2<<",MPQ2="<<mapq2;
                }
                *outp << endl;

                *outp << rid2 << "\t" << flag2 << "\t" << chr2 << "\t" << pos2 << "\t" << mapq2 << "\t" << cigar2;
                *outp << "\t=\t" << pos1 << "\t" << (pos1-pos2) << tail2;
                if (pSitePos) {
                    *outp << "\tHC:Z:FRG1="<<frag2<<",FLG2="<<flag1<<",FRG2="<<frag1<<",MPQ2="<<mapq1;
                }
                *outp << endl;
                /*
                   <readname> <str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2> <mapq1> <mapq2>
                   ofJuiceBox << rid1 <<" "<< (flag1&REVRSE) <<" "<< chr1 <<" "<< pos1 <<" "<< frag1; 
                   ofJuiceBox <<" "<< (flag2&REVRSE) <<" "<< chr2 <<" "<< pos2 <<" "<< frag2;
                   ofJuiceBox << mapq1 <<" "<< mapq2 << endl;
                   */ 
            }
        } else { // [6]:chr [7]:pos [8]:0
            transPE++;
            if (bTrans) {
                *outp << rid1 << "\t" << flag1 << "\t" << chr1 << "\t" << pos1 << "\t" << mapq1 << "\t" << cigar1;
                *outp << "\t" << chr2 << "\t" << pos2 << "\t0" << tail1;
                if (pSitePos) {
                    *outp << "\tHC:Z:FRG1="<<frag1<<",FLG2="<<flag2<<",FRG2="<<frag2<<",MPQ2="<<mapq2;
                }
                *outp << endl;

                *outp << rid2 << "\t" << flag2 << "\t" << chr2 << "\t" << pos2 << "\t" << mapq2 << "\t" << cigar2;
                *outp << "\t" << chr1 << "\t" << pos1 << "\t0" << tail2;
                if (pSitePos) {
                    *outp << "\tHC:Z:FRG1="<<frag2<<",FLG2="<<flag1<<",FRG2="<<frag1<<",MPQ2="<<mapq1;
                }
                *outp << endl;
            }
        }
    }
    cerr << totalPE << endl;
    if (getline(R1, line2)) { cerr << line2 << endl; }

    cerr << "Total\tR1map\tR2map\tPair\tcis\tcis>=cutoff\ttrans\n";
    cerr << totalPE << "\t" << mapR1 <<"\t"<< mapR2 <<"\t"<< mapPE <<"\t"<< cisPE <<"\t"<< cis15K <<"\t"<<transPE<<"\n";
    if (totalPE) {
        fprintf(stderr, "100%\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", 
                mapR1*100/(float)totalPE, mapR2*100/(float)totalPE, mapPE*100/(float)totalPE,
                cisPE*100/(float)totalPE, cis15K*100/(float)totalPE, transPE*100/(float)totalPE);
    }
}
