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
###    PRE.2_DOMAIN_ANNOTATOR.PY
###
###    This script takes in a file with a list of protein domains of interest in the form of
###
###    domain_name \t start_amino_acid_position \t end_amino_acid_position \t domain_length
###
###    The script then takes the biomart transcript information and creates new sequence and genomic position
###    information for the domains of interest. This will be used at a future step to analyze whether variants
###    fall in these domains, and if so what the coding impact is.
###
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

def translate_sequence( nucleotide_sequence ):

	codon_string = ''
	aa_sequence = ''
	aa_counter = 0

	for base in nucleotide_sequence:

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

def identify_biomart_exons_for_domain( biomart_data , domain , domains_dictionary ):

	biomart_lines_list = []

	coding_start_position = domains_dictionary[domain]['CODING_START']
	coding_end_position = domains_dictionary[domain]['CODING_STOP']
	base_gene = domains_dictionary[domain]['BASE_GENE']

	biomart_index_dictionary = {}
	tracker = 0
	for x in biomart_data[0]:
		biomart_index_dictionary[x] = tracker
		tracker += 1

	for line in biomart_data[1:]:

		biomart_CHROM = line[ biomart_index_dictionary['Chrom'] ]
		biomart_TRANS_ID = line[ biomart_index_dictionary['Transcript_ID'] ]
		biomart_START = int(line[ biomart_index_dictionary['Exon_Start'] ])
		biomart_STOP = int(line[ biomart_index_dictionary['Exon_End'] ])

		if line[ biomart_index_dictionary['Coding_Start'] ] == 'None':          ## Note that the transcript file might have '' for an entry, but 'extractVariants()' recodes '' to 'None'
			biomart_CODING_START = line[ biomart_index_dictionary['Coding_Start']]
		else:
			biomart_CODING_START = int(line[ biomart_index_dictionary['Coding_Start'] ])

		if line[ biomart_index_dictionary['Coding_End'] ] == 'None':
			biomart_CODING_END = line[ biomart_index_dictionary['Coding_End']]
		else:
			biomart_CODING_END = int(line[ biomart_index_dictionary['Coding_End'] ])

		biomart_STRAND = int(line[ biomart_index_dictionary['Strand'] ])
		biomart_GENE = line[ biomart_index_dictionary['Gene'] ]
		biomart_SEQUENCE = line[ biomart_index_dictionary['Sequence'] ]

		if biomart_GENE == base_gene:

			if biomart_STRAND == 1:

				## Note that biomart lists all regions in +1 strand reference frame.  Whereas for this code I
				## determined 'domain coding start' as being > 'domain coding end' for -1 strand genes. That makes
				## the > and < evaluations a little counterintuitive for this sections.

				new_line = []
				new_line.append(biomart_CHROM)
				new_line.append(biomart_TRANS_ID)

				# if domain region starts in exon
				if biomart_START <= coding_start_position and coding_start_position <= biomart_STOP and coding_end_position >= biomart_STOP:

					new_line.append(coding_start_position)
					new_line.append(biomart_STOP)
					new_line.append(coding_start_position)
					new_line.append(biomart_CODING_END)
					new_line.append(biomart_STRAND)
					new_line.append(domain)

					start_bases_to_trim = ( coding_start_position - biomart_START )
					new_seq = biomart_SEQUENCE[ start_bases_to_trim: ]

					new_line.append(new_seq)

					biomart_lines_list.append( new_line )

				# if domain region ends in exon
				elif biomart_START <= coding_end_position and coding_end_position <= biomart_STOP and coding_start_position <= biomart_START:

					new_line.append(biomart_START)
					new_line.append(coding_end_position)
					new_line.append(biomart_CODING_START)
					new_line.append(coding_end_position)
					new_line.append(biomart_STRAND)
					new_line.append(domain)

					end_bases_to_trim = ( biomart_STOP - coding_end_position )
					new_seq = biomart_SEQUENCE[ :-(end_bases_to_trim) ]

					new_line.append(new_seq)

					biomart_lines_list.append( new_line )

				# if domain boundaries span exon region
				elif coding_start_position <=  biomart_START and biomart_STOP <= coding_end_position:

					new_line.append(biomart_START)
					new_line.append(biomart_STOP)
					new_line.append(biomart_CODING_START)
					new_line.append(biomart_CODING_END)
					new_line.append(biomart_STRAND)
					new_line.append(domain)
					new_line.append(biomart_SEQUENCE)

					biomart_lines_list.append( new_line )

				# if exon boundaries span domain region
				elif biomart_START <= coding_start_position and coding_end_position <= biomart_STOP:

					new_line.append(coding_start_position)
					new_line.append(coding_end_position)
					new_line.append(coding_start_position)
					new_line.append(coding_end_position)
					new_line.append(biomart_STRAND)
					new_line.append(domain)

					start_bases_to_trim = ( coding_start_position - biomart_START )
					end_bases_to_trim = ( biomart_STOP - coding_end_position )
					new_seq = biomart_SEQUENCE[ start_bases_to_trim:-(end_bases_to_trim) ]

					new_line.append(new_seq)

					biomart_lines_list.append( new_line )

				else:
					None

			elif biomart_STRAND == -1:

				## Note that biomart lists all regions in +1 strand reference frame.  Whereas for this code I
				## determined 'domain coding start' as being > 'domain coding end' for -1 strand genes. That makes
				## the > and < evaluations a little counterintuitive for this sections.

				new_line = []
				new_line.append(biomart_CHROM)
				new_line.append(biomart_TRANS_ID)

				# if exon starts in domain region
				if biomart_START <= coding_end_position and coding_end_position <= biomart_STOP and coding_start_position >= biomart_STOP:

					new_line.append(coding_end_position)
					new_line.append(biomart_STOP)
					new_line.append(coding_end_position)
					new_line.append(biomart_CODING_END)
					new_line.append(biomart_STRAND)
					new_line.append(domain)

					end_bases_to_trim = ( coding_end_position - biomart_START )

					new_seq = biomart_SEQUENCE[ :-(end_bases_to_trim) ]
					new_line.append(new_seq)

					biomart_lines_list.append( new_line )

				# if exon ends in domain region
				elif biomart_START <= coding_start_position and coding_start_position <= biomart_STOP and coding_end_position <= biomart_START:

					new_line.append(biomart_START)
					new_line.append(coding_start_position)
					new_line.append(biomart_CODING_START)
					new_line.append(coding_start_position)
					new_line.append(biomart_STRAND)
					new_line.append(domain)

					start_bases_to_trim = ( biomart_STOP - coding_start_position )

					new_seq = biomart_SEQUENCE[ start_bases_to_trim: ]
					new_line.append(new_seq)

					biomart_lines_list.append( new_line )

				# if domain boundaries span exon region
				elif coding_end_position <=  biomart_START and biomart_STOP <= coding_start_position:

					new_line.append(biomart_START)
					new_line.append(biomart_STOP)
					new_line.append(biomart_CODING_START)
					new_line.append(biomart_CODING_END)
					new_line.append(biomart_STRAND)
					new_line.append(domain)
					new_line.append(biomart_SEQUENCE)

					biomart_lines_list.append( new_line )

				# if exon boundaries span domain region
				elif biomart_START <= coding_end_position and coding_start_position <= biomart_STOP:

					new_line.append(coding_end_position)
					new_line.append(coding_start_position)
					new_line.append(coding_end_position)
					new_line.append(coding_start_position)
					new_line.append(biomart_STRAND)
					new_line.append(domain)

					end_bases_to_trim = ( coding_end_position - biomart_START )
					start_bases_to_trim = ( biomart_STOP - coding_start_position )

					new_seq = biomart_SEQUENCE[ start_bases_to_trim:-(end_bases_to_trim) ]
					new_line.append(new_seq)

					biomart_lines_list.append( new_line )

				else:
					None

	return( biomart_lines_list )

transcript_dir = '../../output/transcript_file'
variant_dir = '../../output/variant_files/'
results_dir = '../../results'
info_dir = '../../info'

transcript_file = '%s/%s' % ( transcript_dir , 'HS_capture_primary_transcript_exon_sequences.txt' )
transcript_vars = extractVariants(transcript_file)

domain_file = '%s/%s' % ( info_dir , 'capture_gene_domain_regions_to_analyze.txt' )
domain_vars = extractVariants( domain_file )

index_tracker = 0
domain_index_dictionary = {}

for x in domain_vars[0]:
	domain_index_dictionary[x]= index_tracker
	index_tracker += 1

domain_list = []

domains_dictionary = {}
for x in domain_vars[1:]:
	DOMAIN = x[ domain_index_dictionary['DOMAIN'] ]
	BASE_GENE = DOMAIN.split('_')[0]
	START = int(x[ domain_index_dictionary['START'] ])
	STOP = int(x[ domain_index_dictionary['STOP'] ])
	LENGTH = int(x[ domain_index_dictionary['LENGTH'] ])

	domain_list.append(DOMAIN)

	domains_dictionary[DOMAIN] = {}
	domains_dictionary[DOMAIN]['BASE_GENE'] = BASE_GENE
	domains_dictionary[DOMAIN]['START'] = START
	domains_dictionary[DOMAIN]['STOP'] = STOP
	domains_dictionary[DOMAIN]['LENGTH'] = LENGTH

REGULATORY_BUFFER_WINDOW = 0

###
### Establish the coding base length for each gene in capture.
###
gene_length_dictionary = {}
sequences_dictionary = {}

positions_list = []

domain_exons_list = []

for domain in domain_list:

	gene_to_analyze = domains_dictionary[domain]['BASE_GENE']

	# Establish transcript dictionaries for current gene being analyzed.

	exon_positions_dictionary,coding_positions_dictionary,gene_position_dictionary, codon_positions_dictionary, strandedness = GenerateTranscriptGenomicPositionDictionaries(transcript_vars , gene_to_analyze , REGULATORY_BUFFER_WINDOW)

	first_AA_in_domain = domains_dictionary[domain]['START']
	last_AA_in_domain = domains_dictionary[domain]['STOP']

	##
	## Codon dictionary is: key = chr_position, value = [ [chr_pos1, chr_pos2, chr_pos3] , aa_# ]
	##

	for key in codon_positions_dictionary.keys():

		if int(codon_positions_dictionary[key][1]) == first_AA_in_domain:
			temp_pos_list = []
			for chr_pos in codon_positions_dictionary[key][0]:
				temp_pos_list.append(int(chr_pos.split('_')[-1]))

			if strandedness == -1:
				start_nuc = max(temp_pos_list)
			elif strandedness == 1:
				start_nuc = min(temp_pos_list)

		elif int(codon_positions_dictionary[key][1]) == last_AA_in_domain:
			temp_pos_list = []
			for chr_pos in codon_positions_dictionary[key][0]:
				temp_pos_list.append(int(chr_pos.split('_')[-1]))

			if strandedness == -1:
				stop_nuc = min(temp_pos_list)
			elif strandedness == 1:
				stop_nuc = max(temp_pos_list)
		else:
			None

	domains_dictionary[domain]['CODING_START'] = int( start_nuc )
	domains_dictionary[domain]['CODING_STOP'] = int( stop_nuc )

	temp_exons_list = identify_biomart_exons_for_domain( transcript_vars , domain , domains_dictionary )

	for x in temp_exons_list:
		domain_exons_list.append( x )

for x in domain_exons_list:
	transcript_vars.append( x )

outPut_NODATE( transcript_vars , '%s/%s' % ( transcript_dir , 'HS_capture_primary_transcript_exon_sequences_WITH_DOMAINS.txt' ) )
