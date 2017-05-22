### Make file for this pipeline

## Install ##
### build reference genomes. 
### sites of restriction enzymes. 
### sites of restriction enzyme for juicebox

### use the install.configure temporarily
#all:
#	make -f bin/Makefile CONFIG_FILE=run_hic.configure CONFIG_SYS=install.configure

##

SITE_POS := annotation/juicebox/site_pos

install: site_pos

site_pos: $(SITE_POS)/hg19_GATC.txt

$(SITE_POS)/hg19_GATC.txt: 
	cd $(SITE_POS) && \
	wget -r --no-directories --no-parent -A "*.txt"  http://enhancer.sdsc.edu/yanxiazh/software/hic-pip/juicebox_site_pos/ && \
	rm robots.txt
