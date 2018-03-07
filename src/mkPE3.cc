/* mkPE2.cc */
// This version does not remove cutting sites >500 from the restriction enzyme.
// and does not output sam/bam files. 
// Needs more tweaking. 
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
    int dist1=abs(pos-*low);
    int dist2=abs(pos-*(low-1));
    assert(low!=vec.end());
    int frag=fragPos[chr][*low];
    return frag;
}


int main(int argc, char* argv[]) {
    int c, cutoff=15000;
    char* pR1=NULL, *pR2=NULL, *pSitePos=NULL, *pSam=NULL, *pTxt=NULL;
    bool bSingle=false, bNear=false, bFar=false, bTrans=false;
    while ((c = getopt (argc, argv, "1:2:c:p:snfto:a:h")) != -1) {
        switch (c) {
            case '1': pR1=optarg; break;
            case '2': pR2=optarg; break;
            case 'p': pSitePos=optarg; break;
            case 'o': pSam=optarg; break;
            case 'a': pTxt=optarg; break;
            case 'c': cutoff=atoi(optarg); break;
            case 's': bSingle=true; break;
            case 'n': bNear=true; break;
            case 'f': bFar=true; break;
            case 't': bTrans=true; break;
            default: usage(argv[0]);
        }
    }

    assert(pR1); ifstream R1(pR1); assert(R1.good());
    assert(pR2); ifstream R2(pR2); assert(R2.good());
    assert(cutoff); assert(bFar);

    ofstream ofSam;
    ofstream ofTxt;
    ostream* outp = &cout;
    ostream* outq = &cout;
    if (pSam) { ofSam.open(pSam); outp=&ofSam; }
    if (pTxt) { ofTxt.open(pTxt); outq=&ofTxt; }

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

    int totalPE=0, mapR1=0, mapR2=0, mapPE=0, cisPE=0, cisGTcutoff=0, transPE=0, cutOR=0; 
    string line1, line2, rid1, rid2, chr1, chr2, tail1, tail2, cigar1, cigar2, rnext1, rnext2;
    int flag1, flag2, pos1, pos2, mapq1, mapq2, pnext1, pnext2, tlen1, tlen2, dist;
    cerr << "Begin mkPE" << endl;
    //start processing each read pair. 
    while (getline(R1, line1)) { assert(!R2.eof());
        getline(R2, line2);
        if (line1[0]=='@') { //*outp << line1 << endl; 
                            continue; }
        totalPE++;

        if ((totalPE%1000000)==0)
            cerr << totalPE << endl;

        stringstream ss1(line1);
        ss1 >> rid1 >> flag1 >> chr1 >> pos1 >> mapq1 >> cigar1 >> rnext1 >> pnext1 >> tlen1; assert(!ss1.fail());
        getline(ss1, tail1);
        stringstream ss2(line2);
        ss2 >> rid2 >> flag2 >> chr2 >> pos2 >> mapq2 >> cigar2 >> rnext2 >> pnext2 >> tlen2; assert(!ss2.fail());
        //make sure the the read 1 name and read 2 name are the same. 
        assert(rid1.compare(rid2)==0);
        getline(ss2, tail2);	
        // read 1 and read 2 mapped rate. 
        if (!(flag1&UNMAPD)) mapR1++;
        if (!(flag2&UNMAPD)) mapR2++;
        // do not report the read if either of the end does not map.
        // unless bSingle is true. 
        if ((flag1&UNMAPD)|(flag2&UNMAPD)) {
//            if (!(flag1&UNMAPD) && bSingle)
//                *outp << line1 << endl;
//            if (!(flag2&UNMAPD) && bSingle)
//                *outp << line2 << endl;
            continue;
        }
        mapPE++;
        //test if the mate is reversed
        int mate1=(flag2 & REVRSE)?MATERV:0;
        int mate2=(flag1 & REVRSE)?MATERV:0;
        //added paired and mate information
        flag1 |= ( PAIRED | mate1 | FIRSTP);
        flag2 |= ( PAIRED | mate2 | SECNDP);
        // remove suplementary and not-primary bit
        flag1 &= ~(NOTPRM | SUPALN);
        flag2 &= ~(NOTPRM | SUPALN);
        // get fragment information
        int frag1=0, frag2=0;
        if (pSitePos) {
            frag1=getFrag(chr1, pos1);
            frag2=getFrag(chr2, pos2);
        }

        //determine how we want to print the alignments.
        bool bPrint=false;
        if (chr1.compare(chr2)==0) {
            cisPE++;
            dist=pos2-pos1;
            rnext1=rnext2="=";
            if (abs(dist)>=cutoff && bFar) {
                cisGTcutoff++; bPrint=true;
            } else if (bNear) {
                bPrint=true;
            }
        } else if (bTrans) {
            transPE++;
            dist=0;
            rnext1=chr2;
            rnext2=chr1;
            bPrint=true;
            }

            if (bPrint) {
           /*     *outp << rid1 << "\t" << flag1 << "\t" << chr1 << "\t" << pos1 << "\t" << mapq1 << "\t" << cigar1;
                *outp << "\t" << rnext1 << "\t" << pos2 << "\t" << (pos2-pos1) << tail1;
                if (pSitePos) {
                    *outp << "\tHC:Z:FRG1="<<frag1<<",FLG2="<<flag2<<",FRG2="<<frag2<<",MPQ2="<<mapq2;
                }
                *outp << endl;

                *outp << rid2 << "\t" << flag2 << "\t" << chr2 << "\t" << pos2 << "\t" << mapq2 << "\t" << cigar2;
                *outp << "\t" << rnext2 << "\t" << pos1 << "\t" << (pos1-pos2) << tail2;
                if (pSitePos) {
                    *outp << "\tHC:Z:FRG1="<<frag2<<",FLG2="<<flag1<<",FRG2="<<frag1<<",MPQ2="<<mapq1;
                }
                *outp << endl; */
                if (chr1<chr2){ 
                *outq << rid1 << "\t" << (flag1&REVRSE) << "\t" << chr1 << "\t" << pos1 
                << "\t" << frag1 << "\t" << (flag2&REVRSE) << "\t" << chr2 << "\t" 
                << pos2 << "\t" << frag2 << "\t" << mapq1 << "\t" << mapq2 << endl;
                } else {
                *outq << rid1 << "\t" << (flag2&REVRSE) << "\t" << chr2 << "\t" << pos2
                << "\t" << frag2 << "\t" << (flag1&REVRSE) << "\t" << chr1 << "\t"
                << pos1 << "\t" << frag1 << "\t" << mapq2 << "\t" << mapq1 << endl;
                }
                   //<readname> <str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2> <mapq1> <mapq2>
            }
    }
    // print the final qc stats
    cerr << totalPE << endl;
    if (getline(R1, line2)) { cerr << line2 << endl; }

    cerr << "Total\tR1map\tR2map\tPair\tcis\tcis>=cutoff\ttrans\tcutPE>=500\n";
    cerr << totalPE << "\t" << mapR1 <<"\t"<< mapR2 <<"\t"<< mapPE <<"\t"<< cisPE <<"\t"<< cisGTcutoff <<"\t"<<transPE << "\t" << cutOR <<"\n";
    if (totalPE) {
        fprintf(stderr, "100%\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n", 
                mapR1*100/(float)totalPE, mapR2*100/(float)totalPE, mapPE*100/(float)totalPE,
                cisPE*100/(float)totalPE, cisGTcutoff*100/(float)totalPE, transPE*100/(float)totalPE), cutOR*100/(float)totalPE;
    }
}
