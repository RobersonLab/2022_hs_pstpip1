#!/usr/bin/env python

##########
# Import #
##########
import numpy
import re
import sys

from burdens_code import *

###################################################################
###################################################################
##
##  02.3_VEP_ANNOTATION_APPEND.PY
##
##  This script takes in the filtered variant lists for the HS capture and/or the gnomAD variants,
##  the VEP annotation file generated on the Ensembl VEP web browser, the transcipt file of interest, and
##  a domains information file. It then adds GENE, amino_acid_change, and HGVS_dna_change information to
##  all of the variants.
##
##  Then in a second step, the code separates variants for both HS and gnomAD variants into lists of
##  missense, synonymous, non-coding, and not-in-a-gene-region variants.
##
##

def append_VEP_annotations( variant_list , annotations_list , gene_list , transcript_list , domains_dictionary ):
	TEST_RETURN = 'EVERYTHING GOOD'

	# Create a header index dictionary for the VEP annotation file so that we can
	# access fields by requesting line[ header_dictionary['field_of_interest'] ]
	# instead of having to know the index number for the field.

	annotation_tracker = 0
	vep_annotation_dict = {}
	for x in annotations_list[0]:
		col = x.strip('#')
		vep_annotation_dict[ col ] = annotation_tracker
		annotation_tracker += 1

	annotations_appended_list = []

	# Create dictionary of variants from VEP annotation file; this will only include variants
	# annotated in the transcript version that we are assessing for the burdens analysis.
	#
	# The dictionary is set up so that keys are 'variant ids' which are all constructed as: 'CHROM_POSITION_REF_ALT'
	# Each key is for a sub-dictionary which contains the information for the 'GENE' the variant is located
	# in, the 'AA_CHANGE', the 'AA_POSITION', and the 'DNA_CHANGE' (hgvs dna change nomenclature.)
	#
	#
	annotations_aa_change_dict = {}

	for var in annotations_list[1:]:
		variation = var[ vep_annotation_dict['Uploaded_variation'] ]
		REF = variation.split('_')[-1].split('/')[0]
		ALT = variation.split('_')[-1].split('/')[1]

		location = var[ vep_annotation_dict['Location'] ]

		CHROM = location.split(':')[0]
		POSITION = location.split(':')[1].split('-')[0]

		var_id = '%s_%s_%s_%s' % ( CHROM , POSITION , REF , ALT)

		gene = var[ vep_annotation_dict['SYMBOL'] ]

		transcript = var[ vep_annotation_dict['Feature'] ]

		temp_DNA_var = var[ vep_annotation_dict['HGVSc'] ]

		if temp_DNA_var == '-':
			HGVS_dna = 'None'
		else:
			if ':' in temp_DNA_var:
				HGVS_dna = temp_DNA_var.split(':')[1]
			else:
				print 'unaccounted for HGVSc annotation'

		temp_aa_change = var[ vep_annotation_dict['HGVSp'] ]

		temp_protein_position = var[ vep_annotation_dict['Protein_position'] ].split('-')

		# start and last aa_position of an indel are separated by '-' in the VEP annotation... so deletion 100-105, deletion of
		#
		# amino acids 100 through 105; default to -9 for non-indel. NOTE, I am keeping track
		# of indel ranges here because later we will use the range to determine if any
		# position in the indel is located in a domain of interest.

		indel_range = -9

		if len(temp_protein_position) == 1:
				protein_position = temp_protein_position[0]
		elif len(temp_protein_position) == 2:
			if temp_protein_position[0] == '' and temp_protein_position[1] == '':
				protein_position = 'None'
			elif '?' in temp_protein_position:   #VEP reports aa_positions for deletions spanning intron-exon junctions as ?-340, or 340-?, and VEP aa_change is reported as '-' since it is not clear how deletion of a splice site will affect protein

				if temp_protein_position[0] == '?' and temp_protein_position[1] == '?':
					print 'oops, somehow both protein positions are \"?\".  Check annotations.'

				if temp_protein_position[0] == '?' and temp_protein_position[1] != '?':
					protein_position = temp_protein_position[1]

				elif temp_protein_position[1] == '?' and temp_protein_position[0] != '?':
					protein_position = temp_protein_position[0]

			elif '?' not in temp_protein_position:
				protein_position = temp_protein_position[0]
				indel_range = '%s-%s' % ( temp_protein_position[0] , temp_protein_position[1] )

		else:
			print 'unaccounted for protein position annotation\t',temp_protein_position

		if temp_aa_change == '-':
			if protein_position == 'None':
				aa_change = 'None'
			else:
				aa_change = 'Splice_Affecting'
		else:
			if len(temp_aa_change.split(':')) == 2:

				# '%3D' is the annotation code returned by VEP for a synonymous variant
				if '%3D' in temp_aa_change.split(':')[1]:
					aa_change = '%s%s' % ( temp_aa_change.split(':')[1].split('%')[0] , '=')
				else:
					aa_change = temp_aa_change.split(':')[1]

			else:
				print 'unaccounted for protein annotation\t',temp_aa_change

		if var_id in annotations_aa_change_dict:
			None

		else:
			if transcript in transcript_list and gene in gene_list:
				annotations_aa_change_dict[ var_id ] = {}
				annotations_aa_change_dict[ var_id ]['GENE'] = gene
				annotations_aa_change_dict[ var_id ]['AA_CHANGE'] = aa_change
				annotations_aa_change_dict[ var_id ]['DNA_CHANGE'] = HGVS_dna

				if indel_range == -9:
					annotations_aa_change_dict[ var_id ]['AA_POSITION'] = protein_position
				else:
					annotations_aa_change_dict[ var_id ]['AA_POSITION'] = indel_range

			else:
				None

	##
	## Now that we have the annotations dictionary.  We can loop through our HS-capture variant list or
	## the gnomAD variant list.  For each variant we determine the variant ID ('CHROM_POSITION_REF_ALT').
	## We can then use that variant ID to extract the VEP annotation information from our annotations dictionary
	## that was generated above.  This annotation info is then appended to the variant information line.
	##
	##

	for variant in variant_list:
		CHROM = variant.info['CHROM']
		POSITION = variant.info['POS']
		temp_ref = variant.info['REF']
		temp_alt = variant.info['ALT']

		# convert variant features to VEP annotation nomenclature;
		# VEP uses dashes for indel annotation.
		if len(temp_ref) == len(temp_alt):
			REF = temp_ref
			ALT = temp_alt
		elif len(temp_ref) > len(temp_alt):
			REF = temp_ref[1:]
			ALT = '-'
		elif len(temp_ref) < len(temp_alt):
			REF = '-'
			ALT = temp_alt[1:]

		temp_variant_ID = '%s_%s_%s_%s' % ( CHROM , POSITION , REF , ALT)

		if temp_variant_ID in annotations_aa_change_dict:
			variant.appendInfo( 'GENE' , annotations_aa_change_dict[ temp_variant_ID ]['GENE'] )
			variant.appendInfo( 'AA_CHANGE' , annotations_aa_change_dict[ temp_variant_ID ]['AA_CHANGE'])
			variant.appendInfo( 'HGVSc' , annotations_aa_change_dict[ temp_variant_ID ]['DNA_CHANGE'])

			current_gene = annotations_aa_change_dict[ temp_variant_ID ]['GENE']
			current_aa_change = annotations_aa_change_dict[ temp_variant_ID ]['AA_CHANGE']

			if annotations_aa_change_dict[ temp_variant_ID ]['AA_POSITION'] == 'None':
				current_aa_position = annotations_aa_change_dict[ temp_variant_ID ]['AA_POSITION']

			elif len( annotations_aa_change_dict[ temp_variant_ID ]['AA_POSITION'].split('-')) == 2:
				current_aa_position = annotations_aa_change_dict[ temp_variant_ID ]['AA_POSITION']

			else:
				current_aa_position = int( annotations_aa_change_dict[ temp_variant_ID ]['AA_POSITION'] )

			if current_aa_position == 'None':
				None

			else:

				# Determine if AA_POSITION or if any position within the AA_POSITION range for an indel
				# is located in one of our domains of interest within a gene of interest.

				for domain in domains_dictionary.keys():

					current_domain_start = int( domains_dictionary[domain]['START'] )
					current_domain_stop = int( domains_dictionary[domain]['STOP'] )

					if len( annotations_aa_change_dict[ temp_variant_ID ]['AA_POSITION'].split('-')) == 1:

						if current_gene in domain and current_aa_position >= current_domain_start and current_aa_position <= current_domain_stop:

							variant.appendInfo( 'GENE' , domain )
							variant.appendInfo( 'AA_CHANGE' , annotations_aa_change_dict[ temp_variant_ID ]['AA_CHANGE'] )

						else:
							None

					elif len( annotations_aa_change_dict[ temp_variant_ID ]['AA_POSITION'].split('-')) == 2:
						#print annotations_aa_change_dict[ temp_variant_ID ]['AA_POSITION']
						aa_first_position = int(current_aa_position.split('-')[0])
						aa_last_position = int(current_aa_position.split('-')[1])

						indel_range_list = range(aa_first_position , (aa_last_position + 1 ))

						if current_gene in domain and any( pos >= current_domain_start and pos <= current_domain_stop for pos in indel_range_list ):

							variant.appendInfo( 'GENE' , domain )
							variant.appendInfo( 'AA_CHANGE' , annotations_aa_change_dict[ temp_variant_ID ]['AA_CHANGE'] )

						else:
							None

		else:
			variant.appendInfo( 'GENE' , 'None' )
			variant.appendInfo( 'AA_CHANGE' , 'None' )
			variant.appendInfo( 'HGVSc' , 'None' )

		annotations_appended_list.append( variant )

	return(annotations_appended_list)

def categorize_variant_effect_type( variant_list , gene_list ):

	# This function just takes in the VEP annotated variant lists, and it splits the list up into
	# separate lists of variants categorized by missense, synonymous, non-coding, non-gene-region.

	temp_NONCODING_list = []
	temp_MISSENSE_list = []
	temp_SYNONYMOUS_list = []
	temp_NO_GENE_list = []    # variants that are on target in that they are within the capture regions, but are not within a gene of interest for the capture panel

	for variant in variant_list:
		aa_change = variant.info['AA_CHANGE'].split('/')[0]
		gene = variant.info['GENE'].split('/')[0]

		if gene in gene_list and aa_change == 'None':
			temp_NONCODING_list.append( variant )
		elif gene in gene_list and aa_change == 'Splice_Affecting':
			temp_MISSENSE_list.append( variant )
		elif gene in gene_list and aa_change[:2] == 'p.':
			if '=' in aa_change:
				temp_SYNONYMOUS_list.append( variant )
			else:
				temp_MISSENSE_list.append( variant )
		else:
			temp_NO_GENE_list.append( variant )
	return( temp_NONCODING_list , temp_MISSENSE_list , temp_SYNONYMOUS_list , temp_NO_GENE_list)

###############################################################
###############################################################
###############################################################

##
##
##  Implement functions above to annotated variant lists.
##
##

#
# Read in file list of analysis cohort.
#

info_dir = '../../info'
output_samples_info_dir = '../../output/samples_info_output'

sites_dir = '../../output/sites_files'
transcript_dir = '../../output/transcript_file'
vep_variant_dir = '../../output/VEP_variant_files'

HS_sites_file = '%s/%s' % ( sites_dir , 'HS_CAPTURE_allQC_variants_capturedHSpatients_only.txt' )
HS_sites_vars = extractVariants( HS_sites_file )
HS_sites_VARIANTS = [VARIANT( x , HS_sites_vars[0] ) for x in HS_sites_vars[1:]]

gnomad_sites_file = '%s/%s' % ( sites_dir , 'gnomad_COMPILED_QC_passed_sites.txt' )
gnomad_sites_vars = extractVariants( gnomad_sites_file )
gnomad_sites_VARIANTS = [VARIANT( x , gnomad_sites_vars[0] ) for x in gnomad_sites_vars[1:]]

HScapture_vep_annotation_file = '%s/%s' % ( vep_variant_dir , 'hs_capture_vep_annotations.txt' )
HScapture_vep_annotation_vars = extractVariants( HScapture_vep_annotation_file )

gnomad_vep_annotation_file = '%s/%s' % ( vep_variant_dir , 'gnomad_vep_annotations.txt' )
gnomad_vep_annotation_vars = extractVariants( gnomad_vep_annotation_file )

transcript_file = '%s/%s' % ( transcript_dir , 'HS_capture_primary_transcript_exon_sequences_WITH_DOMAINS.txt' )
transcript_vars = extractVariants( transcript_file )

domains_file = '%s/%s' % ( info_dir , 'capture_gene_domain_regions_to_analyze.txt' )
domains_vars = extractVariants( domains_file )

#
#  Extract transcript information, gene information and domain information for
#  capture genes of interest.  This information will need to be fed into the
#  annotation function.
#

transcript_tracker = 0
transcript_header_dictionary = {}
for x in transcript_vars[0]:
	transcript_header_dictionary[x] = transcript_tracker
	transcript_tracker += 1

gene_list = []
transcript_list = []

temp_tracker = 0
domain_header_dict = {}
for x in domains_vars[0]:
	domain_header_dict[x] = temp_tracker
	temp_tracker +=1

domains_dictionary = {}
for domain in domains_vars[1:]:
	domain_name = domain[domain_header_dict['DOMAIN']]
	domain_start = domain[domain_header_dict['START']]
	domain_stop = domain[domain_header_dict['STOP']]

	if domain_name in domains_dictionary:
		print 'oops, domain is listed more than once in domain regions file'
	else:
		domains_dictionary[domain_name] = {}
		domains_dictionary[domain_name]['START'] = domain_start
		domains_dictionary[domain_name]['STOP'] = domain_stop

for x in transcript_vars[1:]:
	if x[transcript_header_dictionary['Gene']] in gene_list:
		None
	else:
		gene_list.append( x[transcript_header_dictionary['Gene']] )
		#print gene_list.append( x[transcript_header_dictionary['Gene']] )
	if x[transcript_header_dictionary['Transcript_ID']] in transcript_list:
		None
	else:
		transcript_list.append( x[transcript_header_dictionary['Transcript_ID']] )

#############################################
###  HS VARIANT ANNOTATION AND CATEGORIZATION

hs_annotated_list = append_VEP_annotations( HS_sites_VARIANTS , HScapture_vep_annotation_vars , gene_list , transcript_list , domains_dictionary)

hs_noncoding , hs_missense , hs_synonymous , hs_noGene = categorize_variant_effect_type( hs_annotated_list , gene_list )

if len(hs_annotated_list) > 0:
	hs_annotated_variant_set = VARIANT_SET(hs_annotated_list)
	hs_annotated_return = hs_annotated_variant_set.returnSet()
	outPut_NODATE( hs_annotated_return , '%s/%s' % ( vep_variant_dir , 'HS_CAPTURE_allQC_variants_capturedHSpatients_only_VEP_ANNOTATED.txt' ) )
else:
	print 'No variants to report'

if len(hs_noncoding) > 0:
	hs_NONCODING_variant_set = VARIANT_SET(hs_noncoding)
	hs_NONCODING_return = hs_NONCODING_variant_set.returnSet()
	outPut_NODATE( hs_NONCODING_return , '%s/%s' % ( vep_variant_dir , 'HS_CAPTURE_noncoding_variants_capturedHSpatients_only.txt' ) )
else:
	print 'No variants to report'

if len(hs_missense) > 0:
	hs_MISSENSE_variant_set = VARIANT_SET( hs_missense )
	hs_MISSENSE_return = hs_MISSENSE_variant_set.returnSet()
	outPut_NODATE( hs_MISSENSE_return , '%s/%s' % ( vep_variant_dir , 'HS_CAPTURE_missense_variants_capturedHSpatients_only.txt' ) )
else:
	print 'No variants to report'

if len(hs_synonymous) > 0:
	hs_SYNONYMOUS_variant_set = VARIANT_SET(hs_synonymous)
	hs_SYNONYMOUS_return = hs_SYNONYMOUS_variant_set.returnSet()
	outPut_NODATE( hs_SYNONYMOUS_return , '%s/%s' % ( vep_variant_dir , 'HS_CAPTURE_synonymous_variants_capturedHSpatients_only.txt' ) )
else:
	print 'No variants to report'

if len(hs_noGene) > 0:
	hs_NOGENE_variant_set = VARIANT_SET(hs_noGene)
	hs_NOGENE_return = hs_NOGENE_variant_set.returnSet()
	outPut_NODATE( hs_NOGENE_return , '%s/%s' % ( vep_variant_dir , 'HS_CAPTURE_noGene_variants_capturedHSpatients_only.txt' ) )
else:
	print 'No variants to report'


#############################################
###  GNOMAD VARIANT ANNOTATION AND CATEGORIZATION

gnomad_annotated_list = append_VEP_annotations( gnomad_sites_VARIANTS , gnomad_vep_annotation_vars , gene_list , transcript_list , domains_dictionary)

gnomad_noncoding , gnomad_missense , gnomad_synonymous , gnomad_noGene = categorize_variant_effect_type( gnomad_annotated_list , gene_list )

if len(gnomad_annotated_list) > 0:
	gnomad_annotated_variant_set = VARIANT_SET(gnomad_annotated_list)
	gnomad_annotated_return = gnomad_annotated_variant_set.returnSet()
	outPut_NODATE( gnomad_annotated_return , '%s/%s' % ( vep_variant_dir , 'gnomad_COMPILED_QC_passed_sites_VEP_ANNOTATED.txt' ) )
else:
	print 'No variants to report'

if len(gnomad_noncoding) > 0:
	gnomad_NONCODING_variant_set = VARIANT_SET(gnomad_noncoding)
	gnomad_NONCODING_return = gnomad_NONCODING_variant_set.returnSet()
	outPut_NODATE( gnomad_NONCODING_return , '%s/%s' % ( vep_variant_dir , 'gnomad_NONCODING_variants.txt' ) )
else:
	print 'No variants to report'

if len(gnomad_missense) > 0:
	gnomad_MISSENSE_variant_set = VARIANT_SET( gnomad_missense )
	gnomad_MISSENSE_return = gnomad_MISSENSE_variant_set.returnSet()
	outPut_NODATE( gnomad_MISSENSE_return , '%s/%s' % ( vep_variant_dir , 'gnomad_MISSENSE_variants.txt' ) )
else:
	print 'No variants to report'

if len(gnomad_synonymous) > 0:
	gnomad_SYNONYMOUS_variant_set = VARIANT_SET(gnomad_synonymous)
	gnomad_SYNONYMOUS_return = gnomad_SYNONYMOUS_variant_set.returnSet()
	outPut_NODATE( gnomad_SYNONYMOUS_return , '%s/%s' % ( vep_variant_dir , 'gnomad_SYNONYMOUS_variants.txt' ) )
else:
	print 'No variants to report'

if len(gnomad_noGene) > 0:
	gnomad_NOGENE_variant_set = VARIANT_SET(gnomad_noGene)
	gnomad_NOGENE_return = gnomad_NOGENE_variant_set.returnSet()
	outPut_NODATE( gnomad_NOGENE_return , '%s/%s' % ( vep_variant_dir , 'gnomad_noGene_variants.txt' ) )
else:
	print 'No variants to report'
