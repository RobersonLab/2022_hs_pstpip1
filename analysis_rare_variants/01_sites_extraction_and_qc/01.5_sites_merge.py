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
##  01.5_SITES_MERGE.PY
##
##  This script merges the gnomAD EXOME and GENOME variants information.
##  Additionally, this script calculates the alternate allele frequencies
##  for the merged/compiled variants, and it appends the gnomAD AAF
##  to the HS_CAPTURE variant set.
##

##########################
# Read in variant files
#

sites_dir = '../../output/sites_files'

gnomad_GENOME_SITES_file = '%s/%s' % ( sites_dir , 'gnomad_GENOME_QC_passed_sites.txt')
gnomad_GENOME_vars = extractVariants(gnomad_GENOME_SITES_file)

gnomad_EXOME_SITES_file = '%s/%s' % ( sites_dir , 'gnomad_EXOME_QC_passed_sites.txt')
gnomad_EXOME_vars = extractVariants(gnomad_EXOME_SITES_file)

capture_SITES_file = '%s/%s' % ( sites_dir , 'CAPTURE_QC_passed_sites.txt')
capture_vars = extractVariants(capture_SITES_file)

########################################################################################################
#  Convert variant files into list of VARIANT classes.
#  Create dictionary to be able to hash directly to VARIANT given a key made up of variant information;
#              key = CHROM_POS_REF_ALT (ex,  14_1231412_A_T ) ; value = VARIANT
#
#

GENOME_VARIANTS = [VARIANT(x,gnomad_GENOME_vars[0]) for x in gnomad_GENOME_vars[1:]]
EXOME_VARIANTS = [VARIANT(x,gnomad_EXOME_vars[0]) for x in gnomad_EXOME_vars[1:]]
CAPTURE_VARIANTS = [VARIANT(x , capture_vars[0]) for x in capture_vars[1:]]

genome_dictionary = {}
exome_dictionary = {}
capture_dictionary = {}

for x in GENOME_VARIANTS:
	tmp_ID = '%s_%s_%s_%s' % (x.info['CHROM'] , x.info['POS'] , x.info['REF'] , x.info['ALT'])
	genome_dictionary[tmp_ID] = x

for x in EXOME_VARIANTS:
	tmp_ID = '%s_%s_%s_%s' % (x.info['CHROM'] , x.info['POS'] , x.info['REF'] , x.info['ALT'])
	exome_dictionary[tmp_ID] = x

for x in CAPTURE_VARIANTS:
	tmp_ID = '%s_%s_%s_%s' % (x.info['CHROM'] , x.info['POS'] , x.info['REF'] , x.info['ALT'])
	capture_dictionary[tmp_ID] = x

########################################
#  Merge gnomAD EXOME and GENOME sites.
#

MERGED_VARIANTS = []

IN_BOTH_COUNTER = 0

for ex_variant in EXOME_VARIANTS:
	tmp_ID = '%s_%s_%s_%s' % (ex_variant.info['CHROM'] , ex_variant.info['POS'] , ex_variant.info['REF'] , ex_variant.info['ALT'])

	# If the a variant is present in both the gnomAD EXOME and gnomAD GENOME data sets:
	#     Sum the gnomad EXOME AC counts and the gnomad GENOME AC counts.
	#     Sum the gnomad EXOME AN counts and the gnomad GENOME AN counts.

	if tmp_ID in genome_dictionary:
		IN_BOTH_COUNTER +=1
		gen_variant = genome_dictionary[tmp_ID]

		for field in ex_variant.header_field_list:
			if 'AC_' in field:
				ex_variant.info[field] = int(ex_variant.info[field]) + int(gen_variant.info[field])
			elif 'AN_' in field:
				ex_variant.info[field] = int(ex_variant.info[field]) + int(gen_variant.info[field])
			else:
				None
	else:
		None

	##################################################
	##  Calculate merged gnomAD alt allele frequencies
	##

	AC_nfe = int(ex_variant.info['AC_nfe'])
	AN_nfe = int(ex_variant.info['AN_nfe'])
	AC_afr = int(ex_variant.info['AC_afr'])
	AN_afr = int(ex_variant.info['AN_afr'])
	AC_amr = int(ex_variant.info['AC_amr'])
	AN_amr = int(ex_variant.info['AN_amr'])
	AC_eas = int(ex_variant.info['AC_eas'])
	AN_eas = int(ex_variant.info['AN_eas'])
	AC_sas = int(ex_variant.info['AC_sas'])
	AN_sas = int(ex_variant.info['AN_sas'])

	AC_all_subsets = AC_nfe + AC_afr + AC_amr + AC_eas + AC_sas
	AN_all_subsets = AN_nfe + AN_afr + AN_amr + AN_eas + AN_sas
	#AC_all_subsets = AC_nfe + AC_afr + AC_amr + AC_sas                     ## I believe that I removed gnomad_EAS from the 'all_subsets'
	#AN_all_subsets = AN_nfe + AN_afr + AN_amr + AN_sas                     ## because we didn't have any affected patients classified as EAS by PCA
	                                                                       ## might want to double check that... but we don't end up using this
	                                                                       ## 'all_subsets' value in the burdens analysis anyways.


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

	ex_variant.appendInfo('AC_all_subsets', AC_all_subsets )
	ex_variant.appendInfo('AN_all_subsets', AN_all_subsets )

	ex_variant.appendInfo('aaf_gnomad_nfe', aaf_gnomad_nfe )
	ex_variant.appendInfo('aaf_gnomad_afr', aaf_gnomad_afr )
	ex_variant.appendInfo('aaf_gnomad_amr', aaf_gnomad_amr )
	ex_variant.appendInfo('aaf_gnomad_eas', aaf_gnomad_eas )
	ex_variant.appendInfo('aaf_gnomad_sas', aaf_gnomad_sas )
	ex_variant.appendInfo('aaf_gnomad_all_subsets', aaf_gnomad_all_subsets )

	MERGED_VARIANTS.append( ex_variant )

###################################################################################################################
#
#  Unique EXOME variants and variants in both the EXOME and GENOME data sets have been added to MERGED_VARIANT set.
#  Now we need to add the variants that are unique to the GENOME data set into the MERGED_VARIANT set.
#
#  Process is essentially identical to above.  Notes from the section above also apply to this section.
#

for gen_variant in GENOME_VARIANTS:
	tmp_ID = '%s_%s_%s_%s' % (gen_variant.info['CHROM'] , gen_variant.info['POS'] , gen_variant.info['REF'] , gen_variant.info['ALT'])

	if tmp_ID in exome_dictionary:
		None
	else:
		AC_nfe = int(gen_variant.info['AC_nfe'])
		AN_nfe = int(gen_variant.info['AN_nfe'])
		AC_afr = int(gen_variant.info['AC_afr'])
		AN_afr = int(gen_variant.info['AN_afr'])
		AC_amr = int(gen_variant.info['AC_amr'])
		AN_amr = int(gen_variant.info['AN_amr'])
		AC_eas = int(gen_variant.info['AC_eas'])
		AN_eas = int(gen_variant.info['AN_eas'])
		AC_sas = int(gen_variant.info['AC_sas'])
		AN_sas = int(gen_variant.info['AN_sas'])

		AC_all_subsets = AC_nfe + AC_afr + AC_amr + AC_eas + AC_sas
		AN_all_subsets = AN_nfe + AN_afr + AN_amr + AN_eas + AN_sas
		#AC_all_subsets = AC_nfe + AC_afr + AC_amr + AC_sas
		#AN_all_subsets = AN_nfe + AN_afr + AN_amr + AN_sas

		if float(AN_nfe) == 0.0:
			aaf_gnomad_nfe = '-1'
		else:
			aaf_gnomad_nfe = '%.15f' % ( float(AC_nfe) / float(AN_nfe) )

		if float(AN_afr) == 0.0:
			aaf_gnomad_afr = '-1'
		else:
			aaf_gnomad_afr = '%.15f' % ( float(AC_afr) / float(AN_afr) )

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

		gen_variant.appendInfo('AC_all_subsets', AC_all_subsets )
		gen_variant.appendInfo('AN_all_subsets', AN_all_subsets )

		gen_variant.appendInfo('aaf_gnomad_nfe', aaf_gnomad_nfe )
		gen_variant.appendInfo('aaf_gnomad_afr', aaf_gnomad_afr )
		gen_variant.appendInfo('aaf_gnomad_amr', aaf_gnomad_amr )
		gen_variant.appendInfo('aaf_gnomad_eas', aaf_gnomad_eas )
		gen_variant.appendInfo('aaf_gnomad_sas', aaf_gnomad_sas )
		gen_variant.appendInfo('aaf_gnomad_all_subsets', aaf_gnomad_all_subsets )

		MERGED_VARIANTS.append( gen_variant )

####################################################################
#
#  Add gnomad alt allele frequencies to HS_CAPTURE variant data set.
#     If Capture variant is not present in gnomAD report all subpopulation
#     frequencies as -1.
#

gnomad_merged_dictionary = {}

for x in MERGED_VARIANTS:
	tmp_ID = '%s_%s_%s_%s' % (x.info['CHROM'] , x.info['POS'] , x.info['REF'] , x.info['ALT'])
	gnomad_merged_dictionary[tmp_ID] = x

for cap_variant in CAPTURE_VARIANTS:
	tmp_ID = '%s_%s_%s_%s' % (cap_variant.info['CHROM'] , cap_variant.info['POS'] , cap_variant.info['REF'] , cap_variant.info['ALT'])

	if tmp_ID in gnomad_merged_dictionary:
		merge_variant = gnomad_merged_dictionary[tmp_ID]

		cap_variant.appendInfo('aaf_gnomad_nfe', merge_variant.info['aaf_gnomad_nfe'] )
		cap_variant.appendInfo('aaf_gnomad_afr', merge_variant.info['aaf_gnomad_afr'] )
		cap_variant.appendInfo('aaf_gnomad_amr', merge_variant.info['aaf_gnomad_amr'] )
		cap_variant.appendInfo('aaf_gnomad_eas', merge_variant.info['aaf_gnomad_eas'] )
		cap_variant.appendInfo('aaf_gnomad_sas', merge_variant.info['aaf_gnomad_sas'] )
		cap_variant.appendInfo('aaf_gnomad_all_subsets', merge_variant.info['aaf_gnomad_all_subsets'] )
	else:
		cap_variant.appendInfo('aaf_gnomad_nfe', -1 )
		cap_variant.appendInfo('aaf_gnomad_afr', -1 )
		cap_variant.appendInfo('aaf_gnomad_amr', -1 )
		cap_variant.appendInfo('aaf_gnomad_eas', -1 )
		cap_variant.appendInfo('aaf_gnomad_sas', -1 )
		cap_variant.appendInfo('aaf_gnomad_all_subsets', -1 )

print 'Number of overlapping exome/genome variants:\t%i' % (IN_BOTH_COUNTER)

merged_set = VARIANT_SET(MERGED_VARIANTS)
merged_return = merged_set.returnSet()

cap_set = VARIANT_SET(CAPTURE_VARIANTS)
cap_return = cap_set.returnSet()

outPut_NODATE(merged_return, '%s/%s' % ( sites_dir , 'gnomad_COMPILED_QC_passed_sites.txt' ) )
outPut_NODATE(cap_return, '%s/%s' % ( sites_dir , 'CAPTURE_QC_passed_sites_gnoAAF_appended.txt') )

'''
same_variant_counter = 0
for x in genome_dictionary.keys():
	if x in exome_dictionary.keys():
		same_variant_counter += 1
	else:
		None

print 'number of variants in genome dictionary:\t%i' % (len(genome_dictionary.keys()))
print 'number of variants in exome dictionary:\t%i' % (len(exome_dictionary.keys()))
print 'number of variants in both sets:\5%i' % (same_variant_counter)
'''
