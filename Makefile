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
GENOME_FEATURE := annotation/genome_features/

install: site_pos genome_features domaincall

site_pos: $(SITE_POS)/hg19_GATC.txt

genome_features: $(GENOME_FEATURE)/hg19.GATC.40000.gnf


$(SITE_POS)/hg19_GATC.txt: 
	cd $(SITE_POS) && \
	wget -r --no-directories --no-parent -A "*.txt"  http://enhancer.sdsc.edu/yanxiazh/software/hic-pip/juicebox_site_pos/ && \
	rm robots.txt

$(GENOME_FEATURE)/hg19.GATC.40000.gnf:
	cd $(GENOME_FEATURE) && \
	wget -r -r --no-directories --no-parent -A "*.txt"  http://enhancer.sdsc.edu/yanxiazh/software/hic-pip/genome_features.tar.gz && \
	tar xvf genome_features.tar.gz

