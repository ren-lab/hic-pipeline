### configuration of paths.
mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
BIN = $(shell dirname $(mkfile_path) )

ifdef CONFIG_FILE
include $(CONFIG_FILE)
endif

ifdef CONFIG_SYS
include $(CONFIG_SYS)
endif


FASTQ_R1 = $(FASTQ_PREFIX)_$(R1_SUFFIX)
FASTQ_R2 = $(FASTQ_PREFIX)_$(R2_SUFFIX)
REF = $(BWA_INDEX_PATH)/$(GENOME)
RAW_BAM = $(OUT_NAME).raw.bam

SITE_POS = $(BIN)/../annotation/juicebox/site_pos/$(GENOME)_$(RES_ENZYME).txt
### Begining of dependencies. 
all: config_check init bwa_map
	@ echo okay done

config_check: 
ifndef CONFIG_FILE
	$(error CONFIG_FILE is not provided)
endif 
ifndef CONFIG_SYS
	$(error CONFI_SYS is not provided)
endif

init:
	@ if ! [ -d $(OUT_DIR) ]; then mkdir $(OUT_DIR); fi
	@ if ! [ -d $(OUT_DIR)/log ];then  mkdir $(OUT_DIR)/log; fi

bwa_map: $(FASTQ_R1) $(FASTQ_R2)
	cd $(OUT_DIR) &&\
	$(BIN)/bwa_map.sh -1 $(FASTQ_R1) -2 $(FASTQ_R2) -g $(REF) -p $(N_CPUS) -s $(SAMTOOLS_PATH) -e $(SITE_POS) -o $(RAW_BAM) \
	-t 1000  ## this only test the first 1000 records of the fastq file
