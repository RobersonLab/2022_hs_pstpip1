#!/usr/bin/env python

##########
# Import #
##########
import numpy
import sys

from burdens_code import *

###################################################################
###################################################################
##
##  02.2_ENSEMBL_FORMAT_CONVERSION.PY
##
##
##
##

def filter_out_nonHS_patients( variant_list , cohortSet ):
	passed_variants = []
	filtered_variants = []

	for x in variant_list:
		affected_individuals = set(x.showAltIndividuals())

		if len( affected_individuals.intersection( cohortSet ) ) > 0:
			passed_variants.append(x)
		else:
			filtered_variants.append(x)


	PASS_SET = VARIANT_SET( passed_variants )
	PASS_RETURN = PASS_SET.returnSet()

	FILTER_SET = VARIANT_SET( filtered_variants )
	FILTER_RETURN = FILTER_SET.returnSet()

	return ( PASS_RETURN , FILTER_RETURN )

def variant_format_conversion_for_VEP_submission( variant_file ):
	temp_variant_vars = extractVariants( variant_file )
	temp_VARIANTS = [VARIANT( x , temp_variant_vars[0] ) for x in temp_variant_vars[1:]]

	ensembl_out_list = []

	for var in temp_VARIANTS:
		CHROM = var.info['CHROM']
		START = int(var.info['POS'])
		REF = var.info['REF']
		ALT = var.info['ALT']
		STRAND = '+'

		temp_entry = []

		if len(REF) == len(ALT):
			if len(REF) == 1:
				temp_entry = [ CHROM , START , START , '%s/%s' %( REF , ALT ) , STRAND ]
			elif len(REF) > 1:

				#for multibase substitutions, VEP wants the positions of the first and last changed bases.

				temp_entry = [ CHROM , START , (START + len(ALT) -1) , '%s/%s' %( REF , ALT ) , STRAND ]

		elif len(REF) < len(ALT):
			#for insertions, VEP wants the variant POSITIONS to be reported as START = ( position_of_last_base_before_insertion + 1) , END/STOP = ( position_of_last_base_before_insertion )

			temp_entry = [ CHROM , (START + 1) , START , '%s/%s' %( '-' , ALT[1:] ) , STRAND ]
		elif len(REF) > len(ALT):

			#for deletions, VEP wants the positions of the deleted bases; the POSITION of our sites/variant info is reported as the last reference base that is not deleted.

			temp_entry = [ CHROM , (START + 1) , ( START + len(REF) -1 ) , '%s/%s' %( REF[1:] , '-' ) , STRAND ]

		ensembl_out_list.append(temp_entry)

	return( ensembl_out_list )

#
# Read in file list of analysis cohort.
#

info_dir = '../../info'
output_samples_info_dir = '../../output/samples_info_output'
variant_dir = '../../output/VEP_variant_files'
sites_dir = '../../output/sites_files'

hs_capture_file = '%s/%s' % ( sites_dir , 'HS_CAPTURE_allQC_variants_capturedHSpatients_only.txt' )
gnomad_file = '%s/%s' % ( sites_dir , 'gnomad_COMPILED_QC_passed_sites.txt' )

hs_capture_enseml_format = variant_format_conversion_for_VEP_submission( hs_capture_file )
gnomad_enseml_format = variant_format_conversion_for_VEP_submission( gnomad_file )

outPut_NODATE( hs_capture_enseml_format , '%s/%s' %( variant_dir , 'hs_capture_allsites_ensembl_submission.txt') )
outPut_NODATE( gnomad_enseml_format , '%s/%s' %( variant_dir , 'gnomad_allsites_ensembl_submission.txt') )
