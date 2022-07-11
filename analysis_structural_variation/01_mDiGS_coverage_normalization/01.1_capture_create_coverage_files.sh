#/usr/bin/env bash

##############################################################################################
##############################################################################################
###
###
###    01.1_CAPTURE_CREATE_COVERAGE_FILES.SH
###
###    This script runs samtools depth to evaluate coverage over the targeted capture regions
###    for each individual in the HS cohort.
###
###
###

INFO_DIR=../../info
BAM_DIR=../../data/bam
COVERAGES_DIR=../../output/raw_capture_coverages
SAMTOOLS_DIR=/home/davemh/src/samtools-1.9
BAM_EXTENSION=_rgid_nodup.bam

BED_FILE=$INFO_DIR/HS_capture_fosmid_regions_no_overlaps.bed

for BAM_FILE in `ls $BAM_DIR/*.bam`
do
	SAMPLE=${BAM_FILE#$BAM_DIR/}
	SAMPLE=${SAMPLE%$BAM_EXTENSION}

	OUTFILE=$COVERAGES_DIR/${SAMPLE}_NoDuplicates_Coverages.txt

	$SAMTOOLS_DIR/samtools view -u -f 3 -F 1024 -L $BED_FILE $BAM_FILE | $SAMTOOLS_DIR/samtools depth -b $BED_FILE /dev/stdin > $OUTFILE	
done
