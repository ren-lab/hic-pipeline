all: ../bin/mkPE ../bin/ResEnzymeScan ../bin/sam2mat ../bin/mkPE2 

../bin/mkPE: mkPE.cc
	g++ $^ -o $@
../bin/mkPE2: mkPE2.cc
	g++ $^ -o $@
../bin/ResEnzymeScan: ResEnzymeScan.cc
	g++ $^ -o $@
../bin/sam2mat: sam2mat.cc CImg.h
	g++ sam2mat.cc -o $@ -lpthread
