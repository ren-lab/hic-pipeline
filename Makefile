### Make file for this pipeline

## Install ##
### build reference genomes. 
### sites of restriction enzymes. 
### sites of restriction enzyme for juicebox

### use the install.configure temporarily

all: build install

## build C executives
build: bin/ResEnzymeScan bin/sam2mat bin/mkPE2 bin/mkPE3 #bin/mkPE

#bin/mkPE: src/mkPE.cc
#	g++ $^ -o $@
bin/mkPE2: src/mkPE2.cc
	g++ $^ -o $@
bin/mkPE3: src/mkPE3.cc
	g++ $^ -o $@

bin/ResEnzymeScan: src/ResEnzymeScan.cc
	g++ $^ -o $@
bin/sam2mat: src/sam2mat.cc src/CImg.h
	g++ src/sam2mat.cc -o $@ -lpthread


SITE_POS := annotation/juicebox/site_pos
GENOME_FEATURE := annotation/genome_features/
DOMAIN_CALL:= lib/domaincall_software

install: site_pos genome_features domaincall

site_pos: $(SITE_POS)/hg19_GATC.txt

genome_features: $(GENOME_FEATURE)/hg19.GATC.40000.gnf

domaincall: $(DOMAIN_CALL)/HMM_calls.m

$(SITE_POS)/hg19_GATC.txt: 
	cd $(SITE_POS) && \
	wget -r --no-directories --no-parent -A "*.txt"  http://renlab.sdsc.edu/yanxiao/download/hic-pip/juicebox_site_pos/ && \
	rm robots.txt

$(GENOME_FEATURE)/hg19.GATC.40000.gnf:
	mkdir -p $(GENOME_FEATURE) && \
	cd $(GENOME_FEATURE) && \
	wget -r --no-directories --no-parent  http://renlab.sdsc.edu/yanxiao/download/hic-pip/genome_features.tar.gz && \
	tar xvf genome_features.tar.gz

$(DOMAIN_CALL)/HMM_calls.m:
	cd $(DOMAIN_CALL) && \
	wget -r --no-directories --no-parent  http://renlab.sdsc.edu/yanxiao/download/hic-pip/domaincall_software.tar.gz && \
	tar xvf domaincall_software.tar.gz


