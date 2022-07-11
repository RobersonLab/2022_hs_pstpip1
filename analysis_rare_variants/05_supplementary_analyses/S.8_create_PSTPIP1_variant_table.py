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
##  S.8_CREATE_PSTPIP1_VARIANT_TABLE.PY
##
##  Creates summary table of PSTPIP1 variants.
##
##  Summarizes information such as: Chrom, position, ref, alt, domain, gnomadAAF, number of alts observed for variant,
##  number of individuals affected by variant, etc.
##
##

variants_dir = '../../output/VEP_variant_files'
sites_dir = '../../output/sites_files'
results_dir = '../../results'
samples_info_dir = '../../output/samples_info_output'

#set up Capture files
CAPTURE_missense_file = '%s/%s' % ( variants_dir , 'HS_CAPTURE_missense_variants_capturedHSpatients_only.txt' )
capture_missense_list = extractVariants( CAPTURE_missense_file )

HS_cohort_file = '%s/%s' % ( samples_info_dir , 'all_cohorts_HS_AFFECTED_captured_samples.txt' )
HS_cohort_vars = extractVariants( HS_cohort_file )

missense_VARIANTS = [ VARIANT( x , capture_missense_list[0] ) for x in capture_missense_list[1:] ]

HS_samples_list = [ x[0] for x in HS_cohort_vars[1:] ]

pstpip1_variants = []

for x in missense_VARIANTS:
	if 'PSTPIP1' in x.info['GENE']:
		pstpip1_variants.append( x )
	else:
		None

out_list = []
header = ['CHROMOSOME' , 'POSITION' , 'REF' , 'ALT', 'GENE' , 'DOMAIN' , 'AMINO_ACID_CHANGE' , 'DNA_CHANGE' , 'ALLELES_OBSERVED' , 'INDIVIDUALS_AFFECTED', 'NFE' , 'AFR', 'AMR' , 'SAS' , 'EAS' ]
out_list.append( header )

variant_summary_list = []

for x in pstpip1_variants:
	temp_chrom = x.info['CHROM']
	temp_position = x.info['POS']
	temp_ref = x.info['REF']
	temp_alt = x.info['ALT']
	temp_NFE = x.info['aaf_gnomad_nfe']
	temp_AFR = x.info['aaf_gnomad_afr']
	temp_AMR = x.info['aaf_gnomad_amr']
	temp_SAS = x.info['aaf_gnomad_sas']
	temp_EAS = x.info['aaf_gnomad_eas']
	temp_gene = x.info['GENE'].split('/')[0]
	temp_AAchange = x.info['AA_CHANGE'].split('/')[0]
	temp_dna_change = x.info['HGVSc']

	temp_total_alts = 0
	for patient in HS_samples_list:
		temp_total_alts += x.patientAltCount( patient )

	variant_affecteds = set(x.showAltIndividuals())
	samples_set = set(HS_samples_list)
	temp_total_affected = len(variant_affecteds.intersection(samples_set))

	if len(x.info['GENE'].split('/')) > 1:
		temp_domain = x.info['GENE'].split('/')[-1].split('_')[-2]
	else:
		temp_domain = '--'

	temp_entry = []
	temp_entry.append( temp_chrom )
	temp_entry.append( temp_position )
	temp_entry.append(temp_ref)
	temp_entry.append(temp_alt)
	temp_entry.append(temp_gene)
	temp_entry.append(temp_domain)
	temp_entry.append(temp_AAchange)
	temp_entry.append(temp_dna_change)
	temp_entry.append(temp_total_alts)
	temp_entry.append(temp_total_affected)
	temp_entry.append(temp_NFE)
	temp_entry.append(temp_AFR)
	temp_entry.append(temp_AMR)
	temp_entry.append(temp_SAS)
	temp_entry.append(temp_EAS)

	variant_summary_list.append(temp_entry)

#sort list by genomic position
variant_summary_list.sort(key=lambda z: int(z[1]))

out_list += variant_summary_list

outPut_NODATE( out_list , '%s/%s' % ( results_dir , 'manuscript_PSTPIP1_variant_table.txt') )
