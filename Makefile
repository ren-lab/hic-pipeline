### Make file for this pipeline

## Install ##
### build reference genomes. 
### sites of restriction enzymes. 
### sites of restriction enzyme for juicebox

### use the install.configure temporarily
all:
	make -f bin/Makefile CONFIG_FILE=run_hic.configure CONFIG_SYS=install.configure
