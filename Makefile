### Make file for this pipeline

## Install ##
### build reference genomes. 
### sites of restriction enzymes. 
### sites of restriction enzyme for juicebox

### use the install.configure temporarily

all: build install

## build C executives
build: bin/mkPE bin/ResEnzymeScan bin/sam2mat bin/mkPE2

bin/mkPE: src/mkPE.cc
	g++ $^ -o $@
bin/mkPE2: src/mkPE2.cc
	g++ $^ -o $@
bin/ResEnzymeScan: src/ResEnzymeScan.cc
	g++ $^ -o $@
bin/sam2mat: src/sam2mat.cc src/CImg.h
	g++ src/sam2mat.cc -o $@ -lpthread


SITE_POS := annotation/juicebox/site_pos

install: site_pos

site_pos: $(SITE_POS)/hg19_GATC.txt

$(SITE_POS)/hg19_GATC.txt: 
	cd $(SITE_POS) && \
	wget -r --no-directories --no-parent -A "*.txt"  http://enhancer.sdsc.edu/yanxiazh/software/hic-pip/juicebox_site_pos/ && \
	rm robots.txt

	
