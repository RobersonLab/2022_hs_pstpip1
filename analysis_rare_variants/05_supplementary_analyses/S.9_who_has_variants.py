#!/usr/bin/env python

##########
# Import #
##########
import numpy
import sys

from burdens_code import *

variants_dir = '../../output/VEP_variant_files'
sample_info_dir = '../../output/samples_info_output/'
results_dir = '../../results'

variants_file = '%s/%s' % (variants_dir , 'HS_CAPTURE_missense_variants_capturedHSpatients_only.txt')

missense_vars = extractVariants(variants_file)
mis_VARIANTS = [ VARIANT(x,missense_vars[0]) for x in missense_vars]

cohort_file = '%s/%s' % (sample_info_dir , 'all_cohorts_HS_AFFECTED_captured_samples_withPCAinfo.txt')
cohort_vars = extractVariants(cohort_file)
cohort_dict = {}
for x in cohort_vars[1:]:
	cohort_dict[x[0]] = x[-1]

pstpip1_variant_list = []
header = ['CHROM' , 'POS' , 'REF' , 'ALT' , 'GENE' , 'AA_CHANGE' , 'number_affected' , 'number_alts_observed', 'affected_individuals' , 'aaf_gnomad_nfe', 'aaf_gnomad_afr' , 'aaf_gnomad_amr' , 'aaf_gnomad_sas' , 'aaf_gnomad_eas']
pstpip1_variant_list.append(header)

for variant in mis_VARIANTS:
	if 'PSTPIP1' in variant.info['GENE']:
		raw_affected_individuals = variant.showAltIndividuals()
		filtered_hs_affected_individuals = []
		total_alts_count = 0
		total_individuals_count = 0
		for individual in raw_affected_individuals:
			if individual in cohort_dict.keys():
				individual_entry = []
				genotype = ''
				individual_entry.append(individual)

				if variant.patientAltCount(individual) == 1:
					genotype = '0/1'
				elif variant.patientAltCount(individual) == 2:
					genotype = '1/1'
				else:
					None

				individual_entry.append(genotype)
				individual_entry.append(cohort_dict[individual])

				total_individuals_count += 1
				total_alts_count += variant.patientAltCount(individual)

				filtered_hs_affected_individuals.append( individual_entry )
			else:
				None

		chrom = variant.info['CHROM']
		pos = variant.info['POS']
		ref = variant.info['REF']
		alt = variant.info['ALT']
		gene = variant.info['GENE'].split('/')[-1]
		aa_change = variant.info['AA_CHANGE'].split('/')[0]
		aaf_NFE = variant.info['aaf_gnomad_nfe']
		aaf_AFR = variant.info['aaf_gnomad_afr']
		aaf_AMR = variant.info['aaf_gnomad_amr']
		aaf_SAS = variant.info['aaf_gnomad_sas']
		aaf_EAS = variant.info['aaf_gnomad_eas']

		temp_entry = []
		temp_entry.append(chrom)
		temp_entry.append(pos)
		temp_entry.append(ref)
		temp_entry.append(alt)
		temp_entry.append(gene)
		temp_entry.append(aa_change)
		temp_entry.append(total_individuals_count)
		temp_entry.append(total_alts_count)
		temp_entry.append(filtered_hs_affected_individuals)
		temp_entry.append(aaf_NFE)
		temp_entry.append(aaf_AFR)
		temp_entry.append(aaf_AMR)
		temp_entry.append(aaf_SAS)
		temp_entry.append(aaf_EAS)

		pstpip1_variant_list.append( temp_entry )
	else:
		None

for x in pstpip1_variant_list:
	print x

outPut_NODATE( pstpip1_variant_list , '%s/%s' % (results_dir , 'who_has_PSTPIP1_variants.txt'))
