#/usr/bin/env bash

##################################################################################
##################################################################################
###
###       03.1_GNOMAD_STRUCTURAL_VARIANTS_EXTRACTION.SH
###
###		  This script extracts all of the structural variants from the gnomAD structural
###       variants sites file that fall within the HS-targeted-capture region.
###
###       NOTE: this script only collects variants that overlap with the targeted capture
###       regions based on the 5' position provided in the database. So for example, a variant
###       whose provided positin is 500bp upstream of a capture region will not be included, even
###       if the variant is a 3kb deletion that does overlap with the targeted capture region.
###
###

BCF_PATH=/home/davemh/src/bcftools-1.9

INFO_DIR=../../info
SV_SITES_DIR=../../data/gnomad_structural_variants_sites_file
OUTPUT_SITES_DIR=../../output/gnomad_structural_variants

BED_FILE=$INFO_DIR/HS_capture_fosmid_regions_no_overlaps.bed

STRUCTURAL_VARIANTS_FILE=$SV_SITES_DIR/gnomad_v2.1_sv.sites.vcf.gz
STRUCTURAL_OUT_FILE=$OUTPUT_SITES_DIR/gnomad_structural_variants_in_capture_region.txt

$BCF_PATH/bcftools query -R $BED_FILE -H -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/SVLEN\t%INFO/AF\t%INFO/AC\t%INFO/AN\t%INFO/EUR_AF\t%INFO/AFR_AF\t%INFO/AMR_AF\t%INFO/EAS_AF[\t%FORMAT/GT][\t%FORMAT/DP][\t%FORMAT/GQ]\n' -i 'FILTER="PASS"' $STRUCTURAL_VARIANTS_FILE > $STRUCTURAL_OUT_FILE
