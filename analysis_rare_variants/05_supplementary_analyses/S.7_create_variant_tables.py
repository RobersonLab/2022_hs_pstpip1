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
##  02.3_CREATE_VARIANT_TABLES.PY
##
##  Creates summary tables of genetic variants for the manuscript.
##
##  1. Creates tables with quantitative summaries of gnomAD and HS_capture variants. Total number of variants,
##     in capture, total number of missense/synonymous/noncoding. And also the missense/synonymous/noncoding
##     variant number breakdown for each gene of interest
##
##  2. Creates tables with the actualy HS_capture variant information.  One table each for missense, synonymous, and non-coding.
##     Tables have chrom, position, ref, alt, aa_change, and gnomad_aaf information.
##
##

def variant_organizer( variants_file ):

	table_vars = extractVariants( variants_file )
	table_VARIANTS = [ VARIANT( x , table_vars[0] ) for x in table_vars[1:] ]

	var_out_list = []
	header = ['CHROM' , 'POS' , 'GENE' , 'REF' , 'ALT' , 'AA_CHANGE' , 'HGVSc' , 'aaf_gnomad_afr' , 'aaf_gnomad_amr' , 'aaf_gnomad_eas' , 'aaf_gnomad_nfe' , 'aaf_gnomad_sas' ]
	var_out_list.append( header )

	unsorted_variants_list = []

	for curr_variant in table_VARIANTS:
		temp_entry = []

		for field in header:
			if field == 'GENE' or field == 'AA_CHANGE':
				temp_entry.append( curr_variant.info[ field ].split('/')[0] )
			else:
				temp_entry.append( curr_variant.info[ field ] )

		unsorted_variants_list.append( temp_entry )

	unsorted_variants_list.sort( key = lambda z: ( int(z[0]) , int(z[1]) ) )

	var_out_list = var_out_list + unsorted_variants_list

	return( var_out_list )

def variant_summary_parser( all_sites_file , missense_variants_file , synonymous_variants_file , noncoding_variants_file ):

	all_sites_vars = extractVariants( all_sites_file )
	missense_vars = extractVariants( missense_variants_file )
	synonymous_vars = extractVariants( synonymous_variants_file )
	noncoding_vars = extractVariants( noncoding_variants_file )

	missense_VARIANTS = [ VARIANT( x , missense_vars[0] ) for x in missense_vars[1:] ]
	synonymous_VARIANTS = [ VARIANT( x , synonymous_vars[0] ) for x in synonymous_vars[1:] ]
	noncoding_VARIANTS = [ VARIANT( x , noncoding_vars[0] ) for x in noncoding_vars[1:] ]

	all_VARIANTS_dictionary = {}
	all_VARIANTS_dictionary['missense'] = missense_VARIANTS
	all_VARIANTS_dictionary['synonymous'] = synonymous_VARIANTS
	all_VARIANTS_dictionary['noncoding'] = noncoding_VARIANTS

	all_genes_totals_dictionary = {}
	all_genes_totals_dictionary['total'] = 0
	all_genes_totals_dictionary['missense'] = 0
	all_genes_totals_dictionary['synonymous'] = 0
	all_genes_totals_dictionary['noncoding'] = 0

	variant_summary_dictionary = {}

	for var_type in all_VARIANTS_dictionary.keys():
		for curr_variant in all_VARIANTS_dictionary[var_type]:
			gene = curr_variant.info['GENE'].split('/')[0]

			if gene in variant_summary_dictionary:
				variant_summary_dictionary[ gene ][ var_type ] += 1
				variant_summary_dictionary[ gene ]['total'] += 1

				all_genes_totals_dictionary[ var_type ] += 1
				all_genes_totals_dictionary['total'] += 1
			else:
				variant_summary_dictionary[ gene ] = {}
				variant_summary_dictionary[ gene ]['missense'] = 0
				variant_summary_dictionary[ gene ]['synonymous'] = 0
				variant_summary_dictionary[ gene ]['noncoding'] = 0
				variant_summary_dictionary[ gene ]['total'] = 0

				variant_summary_dictionary[ gene ][ var_type ] += 1
				variant_summary_dictionary[ gene ]['total'] += 1

				all_genes_totals_dictionary[ var_type ] += 1
				all_genes_totals_dictionary['total'] += 1

	header = [ 'GENE' , 'TOTAL', 'NON_CODING' , 'SYNONYMOUS' , 'MISSENSE' ]
	out_list = []
	out_list.append( header )
	all_sites_entry = [ 'ALL_SITES' , '%i' % (len(all_sites_vars[1:])) , 'NA', 'NA', 'NA' ]
	out_list.append( all_sites_entry )

	gene_summary_list = []

	for gene in variant_summary_dictionary.keys():
		temp_total = variant_summary_dictionary[ gene ]['total']
		temp_missense = variant_summary_dictionary[ gene ]['missense']
		temp_synonymous = variant_summary_dictionary[ gene ]['synonymous']
		temp_noncoding = variant_summary_dictionary[ gene ]['noncoding']

		temp_entry = [ gene , temp_total , temp_noncoding , temp_synonymous , temp_missense ]

		gene_summary_list.append(temp_entry)

	all_GENES_entry = [ 'ALL_GENES' , all_genes_totals_dictionary['total'] , all_genes_totals_dictionary['noncoding'] , all_genes_totals_dictionary['synonymous'] , all_genes_totals_dictionary['missense'] ]
	out_list.append( all_GENES_entry )

	gene_summary_list.sort( key = lambda z: z[0])

	out_list = out_list + gene_summary_list

	return( out_list )

##########
##########
##########
variants_dir = '../../output/VEP_variant_files'
sites_dir = '../../output/sites_files'
results_dir = '../../results'

#set up Capture files
CAPTURE_all_QC_passed_sites_file = '%s/%s' % ( sites_dir , 'HS_CAPTURE_allQC_variants_capturedHSpatients_only.txt' )
CAPTURE_missense_file = '%s/%s' % ( variants_dir , 'HS_CAPTURE_missense_variants_capturedHSpatients_only.txt' )
CAPTURE_synonymous_file = '%s/%s' % ( variants_dir , 'HS_CAPTURE_synonymous_variants_capturedHSpatients_only.txt' )
CAPTURE_noncoding_file = '%s/%s' % ( variants_dir , 'HS_CAPTURE_noncoding_variants_capturedHSpatients_only.txt' )

#set up Gnomad files
GNOMAD_all_QC_passed_sites_file = '%s/%s' % ( sites_dir , 'gnomad_COMPILED_QC_passed_sites.txt' )
GNOMAD_missense_file = '%s/%s' % ( variants_dir , 'gnomad_MISSENSE_variants.txt' )
GNOMAD_synonymous_file = '%s/%s' % ( variants_dir , 'gnomad_SYNONYMOUS_variants.txt' )
GNOMAD_noncoding_file = '%s/%s' % ( variants_dir , 'gnomad_NONCODING_variants.txt' )

##########################################################################
###
###  Generate summary statistics for GNOMAD and CAPTURE variants... ie, X # noncoding, X # synonymous , X# missense

CAPTURE_variant_summary = variant_summary_parser( CAPTURE_all_QC_passed_sites_file , CAPTURE_missense_file , CAPTURE_synonymous_file , CAPTURE_noncoding_file )
GNOMAD_variant_summary = variant_summary_parser( GNOMAD_all_QC_passed_sites_file , GNOMAD_missense_file , GNOMAD_synonymous_file , GNOMAD_noncoding_file )

outPut_NODATE( CAPTURE_variant_summary , '%s/%s' % ( results_dir , 'HS_CAPTURE_variant_summary_table.txt' ) )
outPut_NODATE( GNOMAD_variant_summary , '%s/%s' % ( results_dir , 'GNOMAD_variant_summary_table.txt' ) )

##########################################################################
###
###  Generate tables for missense, synonymous, and noncoding variants
capture_missense_list = variant_organizer( CAPTURE_missense_file )
capture_synonymous_list = variant_organizer( CAPTURE_synonymous_file )
capture_noncoding_list = variant_organizer( CAPTURE_noncoding_file )

outPut_NODATE( capture_missense_list , '%s/%s' % ( results_dir , 'capture_MISSENSE_variants_table.txt' ) )
outPut_NODATE( capture_synonymous_list , '%s/%s' % ( results_dir , 'capture_SYNONYMOUS_variants_table.txt' ) )
outPut_NODATE( capture_noncoding_list , '%s/%s' % ( results_dir , 'capture_NONCODING_variants_table.txt' ) )
