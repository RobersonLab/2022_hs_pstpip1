#/usr/bin/env bash

##################################################################################
##################################################################################
###
###       01.1_VARIANT_SITES_EXTRACTION.SH
###
###       gnomAD has different files containing variant information for exomes and genomes. This script
###       extracts the information from those files for genomic positions that fall in the regions of our
###       targeted capture.
###
###       It also extracts on-target variant data from our HS cohort targeted capture.
###
###

BCF_PATH=/home/davemh/src/bcftools-1.9

INFO_DIR=../../info
GNOMAD_GENOME_SITES_DIR=../../data/gnomad_sites_files/genome_sites
GNOMAD_EXOME_SITES_DIR=../../data/gnomad_sites_files/exome_sites
CAPTURE_SITES_DIR=../../data/capture_file
OUTPUT_SITES_DIR=../../output/sites_files

# Inputs
BED_FILE=$INFO_DIR/HS_capture_fosmid_regions_no_overlaps.bed
EXOMES_SITES_FILE=$GNOMAD_EXOME_SITES_DIR/gnomad.exomes.r2.1.1.sites.vcf.bgz
CAPTURE_SITES_FILE=$CAPTURE_SITES_DIR/DJMH_016_029_joint_call.decomp.norm.snpeff.vcf.gz

# Outputs
GENOMES_OUT_FILE=$OUTPUT_SITES_DIR/gnomad_genomes_sites_CAPTURE_REGION.txt
EXOMES_OUT_FILE=$OUTPUT_SITES_DIR/gnomad_exomes_sites_CAPTURE_REGION.txt
CAPTURE_OUT_FILE=$OUTPUT_SITES_DIR/DJMH_016_029_joint_call_sites.txt

##
##  Filter gnomad GENOME sites files to only keep genomic positions within fosmid capture region.
##

FIRST_FILE=True

for GENOMES_SITES_FILE in `ls $GNOMAD_GENOME_SITES_DIR`;
do

	TEMP_FILE=$GNOMAD_GENOME_SITES_DIR/$GENOMES_SITES_FILE

	if [ "$FIRST_FILE" = "True" ]

	then

		$BCF_PATH/bcftools query -R $BED_FILE -H -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/AC\t%INFO/AN\t%INFO/AC_nfe\t%INFO/AN_nfe\t%INFO/AC_nfe_nwe\t%INFO/AN_nfe_nwe\t%INFO/AC_afr\t%INFO/AN_afr\t%INFO/AC_amr\t%INFO/AN_amr\t%INFO/AC_eas\t%INFO/AN_eas[\t%FORMAT/GT][\t%FORMAT/DP][\t%FORMAT/GQ]\n' -i 'FILTER="PASS"' $TEMP_FILE > $GENOMES_OUT_FILE
		FIRST_FILE=False

	else

		$BCF_PATH/bcftools query -R $BED_FILE  -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/AC\t%INFO/AN\t%INFO/AC_nfe\t%INFO/AN_nfe\t%INFO/AC_nfe_nwe\t%INFO/AN_nfe_nwe\t%INFO/AC_afr\t%INFO/AN_afr\t%INFO/AC_amr\t%INFO/AN_amr\t%INFO/AC_eas\t%INFO/AN_eas[\t%FORMAT/GT][\t%FORMAT/DP][\t%FORMAT/GQ]\n' -i 'FILTER="PASS"' $TEMP_FILE >> $GENOMES_OUT_FILE

	fi

done

##
##  Filter gnomad EXOME sites files to only keep genomic positions within fosmid capture region.
##

$BCF_PATH/bcftools query -R $BED_FILE -H -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/AC\t%INFO/AN\t%INFO/AC_nfe\t%INFO/AN_nfe\t%INFO/AC_nfe_nwe\t%INFO/AN_nfe_nwe\t%INFO/AC_afr\t%INFO/AN_afr\t%INFO/AC_amr\t%INFO/AN_amr\t%INFO/AC_eas\t%INFO/AN_eas\t%INFO/AC_sas\t%INFO/AN_sas[\t%FORMAT/GT][\t%FORMAT/DP][\t%FORMAT/GQ]\n' -i 'FILTER="PASS"' $EXOMES_SITES_FILE > $EXOMES_OUT_FILE

##
##  Filter HS TARGETED CAPTURE sites files to only keep genomic positions within fosmid capture region; this removes off target sites
##

$BCF_PATH/bcftools query -R $BED_FILE -H -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%GT][\t%DP][\t%GQ]\n' -i 'FILTER="PASS"' $CAPTURE_SITES_FILE > $CAPTURE_OUT_FILE
