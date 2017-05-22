### configuration of paths.
mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
BIN = $(shell dirname $(mkfile_path) )

ifdef CONFIG_FILE
include $(CONFIG_FILE)
endif

ifdef CONFIG_SYS
include $(CONFIG_SYS)
endif

## File names and paths
FASTQ_R1 = $(FASTQ_PREFIX)_$(R1_SUFFIX)
FASTQ_R2 = $(FASTQ_PREFIX)_$(R2_SUFFIX)
REF = $(BWA_INDEX_PATH)/$(GENOME)
RAW_BAM = $(NAME).raw.bam
DEDUP_BAM = $(NAME).dedup.bam
FINAL_BAM = $(NAME).final.bam
SITE_POS = $(BIN)/../annotation/juicebox/site_pos/$(GENOME)_$(RES_ENZYME).txt
ANNOTATION = $(BIN)/../annotation/

#Functions
#TIMEIT = time

### Begining of dependencies. 
all: config_check init bwa_map dedup_bam filter_cutsite_bam bam2hic
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
	@ if ! [ -d $(OUT_DIR)/qc ];then  mkdir $(OUT_DIR)/qc; fi


bwa_map: $(OUT_DIR)/qc/$(NAME).raw.flagstat
dedup_bam: $(OUT_DIR)/qc/$(NAME).dedup.flagstat
filter_cutsite_bam: $(OUT_DIR)/qc/$(NAME).final.flagstat
bam2hic: $(OUT_DIR)/$(NAME).hic

$(OUT_DIR)/qc/$(NAME).raw.flagstat: $(FASTQ_R1) $(FASTQ_R2)
	cd $(OUT_DIR) &&\
	$(TIMEIT) $(BIN)/bwa_map.sh -1 $(FASTQ_R1) -2 $(FASTQ_R2) -g $(REF) -p $(N_CPUS) -s $(SAMTOOLS_PATH) -e $(SITE_POS) -n $(NAME) \
	-t 1000  ## this only test the first 1000 records of the fastq file

$(OUT_DIR)/qc/$(NAME).dedup.flagstat: $(OUT_DIR)/qc/$(NAME).raw.flagstat
	cd $(OUT_DIR) &&\
	$(TIMEIT) $(BIN)/bam_dedup.sh -n $(NAME) -s $(SAMTOOLS_PATH)
	
$(OUT_DIR)/qc/$(NAME).final.flagstat: $(OUT_DIR)/qc/$(NAME).dedup.flagstat
	cd $(OUT_DIR) &&\
	$(TIMEIT) $(BIN)/bam_filter_by_cutsite.sh -n $(NAME) -c $(ANNOTATION)/cut_sites/$(GENOME).$(RES_ENZYME).bed.gz -p $(N_CPUS) -s $(SAMTOOLS_PATH)

$(OUT_DIR)/$(NAME).hic: $(OUT_DIR)/$(FINAL_BAM)
	cd $(OUT_DIR) &&\
	$(TIMEIT) $(BIN)/bam_to_hic.sh -n $(NAME) -g $(GENOME) -s $(SAMTOOLS_PATH) -p $(N_CPUS) -c $(SITE_POS)


