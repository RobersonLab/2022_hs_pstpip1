#/usr/bin/env bash

###########################################################################################
###########################################################################################
###
###    S.4_CAPTURE_READ_ANALYSIS.SH
###
###    This script runs samtools view -c with various parameters on bam files from all HS capture
###    individuals to extract total read counts and on-target read counts in order to evaluate
###    capture on-target performance.
### 
###
###

INFO_DIR=../../info
BAM_DIR=../../data/bam
SAMTOOLS_DIR=/home/davemh/src/samtools-1.9
COVERAGES_DIR=../../output/capture_coverages

BAM_EXTENSION=_rgid_nodup.bam

CAPTURE_BED=$INFO_DIR/HS_capture_fosmid_regions_no_overlaps.bed

ANALYSIS_OUTFILE=$COVERAGES_DIR/HS_capture_on_target_and_pcr_duplicates_ALL_SAMPLES.txt

#SAMPLE_FILE=$PROJECT_PATH/$INFO_DIR/sample_test.txt

echo -e 'SAMPLE\tTOTAL_READS\tREADS_ON_TARGET\tDUPLICATE_READS\tTOTAL_READS_NO_DUPLICATES\tREADS_ON_TARGET_NO_DUPLICATES' > $ANALYSIS_OUTFILE

for BAM_FILE in `ls $BAM_DIR/*.bam`
do
	SAMPLE=${BAM_FILE#$BAM_DIR/}
	SAMPLE=${SAMPLE%$BAM_EXTENSION}
	
	READS__MAPPED_PROPER_PAIR=`$SAMTOOLS_DIR/samtools view -c -f 3 $BAM_FILE`
	READS__MAPPED_PROPER_PAIR_DUPLICATE=`$SAMTOOLS_DIR/samtools view -c -f 1027 $BAM_FILE`
	READS__MAPPED_PROPER_PAIR_ON_TARGET=`$SAMTOOLS_DIR/samtools view -c -f 3 -L $CAPTURE_BED $BAM_FILE`
	READS__MAPPED_PROPER_PAIR_NO_DUPS=`$SAMTOOLS_DIR/samtools view -c -f 3 -F 1024 $BAM_FILE`
	READS__MAPPED_PROPER_PAIR_ON_TARGET_NO_DUPS=`$SAMTOOLS_DIR/samtools view -c -f 3 -F 1024 -L $CAPTURE_BED $BAM_FILE`

	echo -e $SAMPLE'\t'$READS__MAPPED_PROPER_PAIR'\t'$READS__MAPPED_PROPER_PAIR_ON_TARGET'\t'$READS__MAPPED_PROPER_PAIR_DUPLICATE'\t'$READS__MAPPED_PROPER_PAIR_NO_DUPS'\t'$READS__MAPPED_PROPER_PAIR_ON_TARGET_NO_DUPS >> $ANALYSIS_OUTFILE
done

#samtools view -c -f 3 -F1027 -L ./info/HS_capture_fosmid_regions_no_overlaps.bed ./data/bam/OSU_HS_2018-08-J_rgid_nodup.bam
#samtools view -uf 3 ./data/bam/OSU_HS_2018-08-J_rgid_nodup.bam | coverageBed -sorted -g ./info/chromList -a ./info/HS_capture_fosmid_regions_no_overlaps.bed -b stdin  > test_cov.txt
#samtools view -u -f 3 -F1027 ./data/bam/OSU_HS_2018-08-J_rgid_nodup.bam | samtools view -c -F 1027 -L./info/HS_capture_fosmid_regions_no_overlaps.bed ./data/bam/OSU_HS_2018-08-J_rgid_nodup.bam
