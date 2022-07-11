#!/usr/bin/env python

##########
# Import #
##########
import numpy
import sys

from burdens_code import *

#################################################################################################
#################################################################################################
###
###    PRE.3_TRANSCRIPT_TABLE_GENERATOR.PY
###
###    This script compiles the information from the biomart transcript information file into a table
###    for the paper which report the genes/protein domains analyzed in the variant burdens analysis,
###    their transcript_ID's, their start and stop amino acid positions, and protein lengths.
###
###

def generate_cds_sequence( coding_positions_dictionary ):
	ordered_positions_list = []
	for key in coding_positions_dictionary.keys():
		chrom = key.split('_')[0]
		position = int(key.split('_')[-1])
		ordered_positions_list.append( [ key , position , chrom ] )
	ordered_positions_list.sort( key=lambda z: z[1])

	sequence = ''
	for x in ordered_positions_list:
		sequence += coding_positions_dictionary[x[0]]

	#print sequence
	return( sequence )

def translate_sequence( nucleotide_sequence , STRAND ):

	codon_string = ''
	aa_sequence = ''
	aa_counter = 0

	if STRAND == 1:
		sequence_to_translate = nucleotide_sequence
	elif STRAND == -1:
		sequence_to_translate = reversed( nucleotide_sequence )

	for base in sequence_to_translate:

		codon_string += base

		if len(codon_string)>0 and len(codon_string)%3 == 0:

			aa_sequence+= codon_dictionary[codon_string]
			aa_counter += 1
			codon_string = ''
		else:
			None

	if aa_sequence[-1] == '*':
		aa_counter -=1
	else:
		None

	if len( codon_string ) > 0:
		print 'oops, something is wrong. there are leftover nucleotides'
	else:
		None

	#print 'aa_length: %i' % (aa_counter)
	#print aa_sequence
	return( aa_sequence , aa_counter )

def determine_positions_of_domains( sequences_dictionaries , PROTEIN, DOMAIN ):

	PROTEIN_sequence = sequences_dictionaries[PROTEIN]['amino_acid_sequence']
	PROTEIN_length = len(PROTEIN_sequence)

	DOMAIN_sequence = sequences_dictionaries[DOMAIN]['amino_acid_sequence']

	temp_protein_sequence = PROTEIN_sequence

	first_base_position = 0
	while DOMAIN_sequence in temp_protein_sequence:
		first_base_position += 1
		temp_protein_sequence = temp_protein_sequence[1:]

	temp_protein_sequence = PROTEIN_sequence
	last_base_position = (PROTEIN_length + 1)
	while DOMAIN_sequence in temp_protein_sequence:
		last_base_position -= 1
		temp_protein_sequence = temp_protein_sequence[:-1]

	return( first_base_position , last_base_position)

#############################################
#############################################

transcript_dir = '../../output/transcript_file'
variant_dir = '../../output/variant_files/'
results_dir = '../../results'

transcript_file = '%s/%s' % ( transcript_dir , 'HS_capture_primary_transcript_exon_sequences_WITH_DOMAINS.txt' )
transcript_vars = extractVariants(transcript_file)

index_tracker = 0
transcript_index_dictionary = {}

for x in transcript_vars[0]:
	transcript_index_dictionary[x]= index_tracker
	index_tracker += 1

gene_list = []
gene_info_dictionary = {}

for x in transcript_vars[1:]:
	biomart_GENE = x[transcript_index_dictionary['Gene']]
	biomart_CHROM = x[transcript_index_dictionary['Chrom']]
	biomart_TRANS_ID = x[transcript_index_dictionary['Transcript_ID']]

	if biomart_GENE in gene_list:
		None
	else:
		gene_list.append( biomart_GENE)

		gene_info_dictionary[biomart_GENE] = {}
		gene_info_dictionary[biomart_GENE]['chrom'] = biomart_CHROM
		gene_info_dictionary[biomart_GENE]['transcript_id'] = biomart_TRANS_ID

REGULATORY_BUFFER_WINDOW = 0

###
### Establish the coding base length for each gene in capture.
###

gene_length_dictionary = {}

for gene_to_analyze in gene_list:
	# Establish transcript dictionaries for current gene being analyzed.

	exon_positions_dictionary,coding_positions_dictionary,gene_position_dictionary, codon_positions_dictionary, strandedness = GenerateTranscriptGenomicPositionDictionaries(transcript_vars,gene_to_analyze,REGULATORY_BUFFER_WINDOW)

	gene_length_dictionary[gene_to_analyze] = len(coding_positions_dictionary)

	temp_nuc_sequence = generate_cds_sequence( coding_positions_dictionary )
	temp_aa_sequence , aa_length = translate_sequence( temp_nuc_sequence , strandedness )

	if temp_aa_sequence[-1] == '*':
		temp_aa_sequence = temp_aa_sequence[:-1]
	else:
		None

	gene_info_dictionary[gene_to_analyze]['nucleotide_sequence'] = temp_nuc_sequence
	gene_info_dictionary[gene_to_analyze]['amino_acid_sequence'] = temp_aa_sequence
	gene_info_dictionary[gene_to_analyze]['protein_length'] = aa_length

temp_out_list = []

for gene_to_analyze in gene_list:

	first_aa , last_aa = determine_positions_of_domains( gene_info_dictionary , gene_to_analyze.split('_')[0], gene_to_analyze )

	temp_trans_id = gene_info_dictionary[gene_to_analyze]['transcript_id']
	temp_protein_length = gene_info_dictionary[gene_to_analyze]['protein_length']

	temp_entry = []
	temp_entry.append(gene_to_analyze)
	temp_entry.append(temp_trans_id)
	temp_entry.append(first_aa)
	temp_entry.append(last_aa)
	temp_entry.append(temp_protein_length)

	temp_out_list.append( temp_entry )

	#print '%s\t%s\t%i\t%i\t%i' % (gene_to_analyze , temp_trans_id , first_aa , last_aa , temp_protein_length)

temp_out_list.sort(key=lambda z: z[0])

header = [['GENE','TRANSCRIPT_ID','START_AMINO_ACID','END_AMINO_ACID','PROTEIN_LENGTH']]

final_out = header + temp_out_list

outPut_NODATE( final_out , '%s/%s' % ( results_dir , 'manuscript_transcript_table.txt') )
