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
##  03.2_PCA_ANALYSIS_gnomad_AAF_append.PY
##
##  This script appends gnomAD exome alt allele frequencies to the
##  HS targeted capture sites dataset.  For performing PCA ancestry
##  analysis, we need both the genotypes of the HS targeted capture
##  individuals, as well as the gnomAD subpopulation alt allele frequencies
##  so that we can simulate individual genotypes from gnomAD information.
##
##

##########################
# Read in variant files
#

pca_dir = '../../output/pca_results'

gnomad_EXOME_variants_file = '%s/%s' % ( pca_dir , 'PCA_gnomad_EXOME_PASSED_variants.txt')
gnomad_EXOME_vars = extractVariants(gnomad_EXOME_variants_file)

capture_sites_file = '%s/%s' % ( pca_dir , 'PCA_capture_PASSED_variants_capturedHSpatients_only.txt')
capture_sites_vars = extractVariants(capture_sites_file)

########################################################################################################
#  Convert variant files into list of VARIANT classes.
#  Create dictionary to be able to hash directly to VARIANT given a key made up of variant information;
#              key = CHROM_POS_REF_ALT (ex,  14_1231412_A_T ) ; value = VARIANT
#
#

EXOME_VARIANTS = [VARIANT(x,gnomad_EXOME_vars[0]) for x in gnomad_EXOME_vars[1:]]
CAPTURE_VARIANTS = [VARIANT(x , capture_sites_vars[0]) for x in capture_sites_vars[1:]]

exome_dictionary = {}
capture_dictionary = {}

for x in EXOME_VARIANTS:
	tmp_ID = '%s_%s_%s_%s' % (x.info['CHROM'] , x.info['POS'] , x.info['REF'] , x.info['ALT'])
	exome_dictionary[tmp_ID] = x

for x in CAPTURE_VARIANTS:
	tmp_ID = '%s_%s_%s_%s' % (x.info['CHROM'] , x.info['POS'] , x.info['REF'] , x.info['ALT'])
	capture_dictionary[tmp_ID] = x

####################################################################
#
#  Add gnomad alt allele frequencies to HS_CAPTURE variant data set.
#     If Capture variant is not present in gnomAD report all subpopulation
#     frequencies as -1.
#

OUT_VARIANTS = []

for cap_variant in CAPTURE_VARIANTS:
	tmp_ID = '%s_%s_%s_%s' % (cap_variant.info['CHROM'] , cap_variant.info['POS'] , cap_variant.info['REF'] , cap_variant.info['ALT'])

	if tmp_ID in exome_dictionary:
		exome_variant = exome_dictionary[tmp_ID]

		AC_nfe = int(exome_variant.info['AC_nfe'])
		AN_nfe = int(exome_variant.info['AN_nfe'])
		AC_afr = int(exome_variant.info['AC_afr'])
		AN_afr = int(exome_variant.info['AN_afr'])
		AC_amr = int(exome_variant.info['AC_amr'])
		AN_amr = int(exome_variant.info['AN_amr'])
		AC_eas = int(exome_variant.info['AC_eas'])
		AN_eas = int(exome_variant.info['AN_eas'])
		AC_sas = int(exome_variant.info['AC_sas'])
		AN_sas = int(exome_variant.info['AN_sas'])

		AC_all_subsets = AC_nfe + AC_afr + AC_amr + AC_eas + AC_sas
		AN_all_subsets = AN_nfe + AN_afr + AN_amr + AN_eas + AN_sas

		if float(AN_nfe) == 0.0:
			aaf_gnomad_nfe = '-1'
		else:
			aaf_gnomad_nfe = '%.15f' % ( float(AC_nfe) / float(AN_nfe) )      ## Not sure if we need to adjust this to acount for proper significant figures.
																			  ## Right now I just have it set to report 15 decimal places so as not to accidentally
		if float(AN_afr) == 0.0:                                              ## throw out a value.  But theoretically if 1 alt was found in 400,000 alleles,
			aaf_gnomad_afr = '-1'                                             ## the frequency would be 5.0 x 10^-6.  (So for that if we did %.4f, I believe it would)
		else:                                                                 ## end up just reporting an AAF of 0.0000, right?)  We could adjust the decimal places reported
			aaf_gnomad_afr = '%.15f' % ( float(AC_afr) / float(AN_afr) )      ## based on the total alleles count, and that would be a way to adjust significant figures
																			  ## for each variant and each subpopulation within each variant.
		if float(AN_sas) == 0.0:
			aaf_gnomad_sas = '-1'
		else:
			aaf_gnomad_sas = '%.15f' % ( float(AC_sas) / float(AN_sas) )

		if float(AN_eas) == 0.0:
			aaf_gnomad_eas = '-1'
		else:
			aaf_gnomad_eas = '%.15f' % ( float(AC_eas) / float(AN_eas) )

		if float(AN_amr) == 0.0:
			aaf_gnomad_amr = '-1'
		else:
			aaf_gnomad_amr = '%.15f' % ( float(AC_amr) / float(AN_amr) )

		if float(AN_all_subsets) == 0.0:
			aaf_gnomad_all_subsets = '-1'
		else:
			aaf_gnomad_all_subsets = '%.15f' % ( float(AC_all_subsets) / float(AN_all_subsets) )

		########################################################
		#
		#  Append calculated alt_allele_frequencies to VARIANT information
		#

		cap_variant.appendInfo('AC_all_subsets', AC_all_subsets )
		cap_variant.appendInfo('AN_all_subsets', AN_all_subsets )

		cap_variant.appendInfo('aaf_gnomad_nfe', aaf_gnomad_nfe )
		cap_variant.appendInfo('aaf_gnomad_afr', aaf_gnomad_afr )
		cap_variant.appendInfo('aaf_gnomad_amr', aaf_gnomad_amr )
		cap_variant.appendInfo('aaf_gnomad_eas', aaf_gnomad_eas )
		cap_variant.appendInfo('aaf_gnomad_sas', aaf_gnomad_sas )
		cap_variant.appendInfo('aaf_gnomad_all_subsets', aaf_gnomad_all_subsets )

		OUT_VARIANTS.append( cap_variant )

	else:

		None

out_set = VARIANT_SET(OUT_VARIANTS)
out_return = out_set.returnSet()

outPut_NODATE( out_return , '%s/%s' % ( pca_dir , 'PCA_capture_PASSED_variants_capturedHSpatients_only__gnomad_aaf_appended.txt' ) )
