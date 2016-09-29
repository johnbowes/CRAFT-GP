#!/bin/bash

#
# LOAD iCSF SOFTWARE MODULES
#
module load apps/gcc/ruby/2.2.3
module load apps/binapps/anaconda/2.3.0
module load apps/gcc/R/3.3.0
module load apps/binapps/ensembl-api/84
module load apps/gcc/tabix/0.2.6

#
# DEFINE USER INPUT VARIABLES
#
NAME='test'
INDEX_SNP_FILE='input/test/test_index_snps.txt'		
SUMMARY_STATS_FILE='input/test/test_summary_stats.txt'
DISTANCE_TYPE='-m'
DISTANCE='0.1'
CASES='1962'
CONTROLS='8923'
EPIGENOMES='CD8_Naive_Primary_Cells'
EPIGENOMES_TYPE='list'

#
# CREATE OUTPUT DIRECTORY FOR RESULTS
#
OUTPUT_DIR="${NAME}_results"
mkdir -p ${OUTPUT_DIR}/{regions,credible_snps,plots,annotation,bed}

#
# CALCULATE GENOMIC REGIONS BASED ON DEFINED DISTANCE
#
OUTPUT="${OUTPUT_DIR}/regions/"
ruby scripts/define_regions_main.rb -i $INDEX_SNP_FILE $DISTANCE_TYPE $DISTANCE -o $OUTPUT

# remove unnecessary dir
rm -r ${OUTPUT}supplementary/

#
# CALCULATE CREDIBLE SNP SETS
#

# subset GWAS data
REGIONS="${OUTPUT_DIR}/regions/region_boundaries_01cm.txt"
GWAS_SUBSET="${OUTPUT_DIR}/credible_snps/summary_stats_subset.txt"

python scripts/filter_summary_stats.py \
	--regions $REGIONS \
	--stats $SUMMARY_STATS_FILE \
	--out $GWAS_SUBSET

# SNP sets
OUTPUT="${OUTPUT_DIR}/credible_snps/"
BED_OUTPUT="${OUTPUT_DIR}/bed/"

Rscript --vanilla scripts/credible_snps_main.R \
	-r $REGIONS \
	-a $CASES \
	-u $CONTROLS \
	-s $GWAS_SUBSET \
	-o $OUTPUT \
	-b $BED_OUTPUT

# clean up
rm $GWAS_SUBSET

#
# ANNOTATION
#
CREDIBLE_SNPS="${OUTPUT_DIR}/credible_snps/credible_snps_0.99.txt"
OUTPUT="${OUTPUT_DIR}/annotation/annotation"

python scripts/annotation.py \
	--input $CREDIBLE_SNPS \
	--output $OUTPUT \
	--epi_names $EPIGENOMES \
	--epi_type $EPIGENOMES_TYPE

#
# VISUALISATION
#
REGIONS="${OUTPUT_DIR}/credible_snps/summary_table_0.99.txt"
SNPS="${OUTPUT_DIR}/credible_snps/credible_snps_0.99.txt"
OUTPUT="${OUTPUT_DIR}/plots/"
EPIGENOMES="${OUTPUT_DIR}/annotation/annotation.epigenomes"

Rscript --vanilla scripts/visualisation_main.R \
	-r ${REGIONS} \
	-s ${SNPS} \
	-e $EPIGENOMES \
	-o ${OUTPUT}

#
# COMPRESS RESULTS
#
tar -czvf ${NAME}_results.tar.gz ${NAME}_results
