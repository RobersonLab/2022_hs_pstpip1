#!/usr/bin/env python

##########
# Import #
##########
import numpy
import re
import sys

from burdens_code import *

##################################################################################
##################################################################################
###
###   03.3_exon_overlapping_HS_CNV.py
###
###   This script takes the called CNV from our HS-targeted capture and identifies
###   any in the list that affect regions that overlap with exonic regions in our
###   target-capture genes of interest list.
###
###

gnomad_sv_vars_dir = '../../output/gnomad_structural_variants'
info_dir = '../../info'
cnv_results_dir = '../../results/mDiGS_CNV'
results_dir = '../../results'

transcript_file = '%s/%s' % ( info_dir , 'HS_capture_primary_transcript_exon_sequences.txt' )
transcript_vars = extractVariants( transcript_file )

cnv_file = '%s/%s' % ( cnv_results_dir , 'manual_correction_ALL_BATCHES_COMPILED_cnv_distributions.txt' )
cnv_vars = extractVariants( cnv_file )

transcript_index_dictionary = {}
index_tracker = 0
for x in transcript_vars[0]:
	transcript_index_dictionary[x] = index_tracker
	index_tracker += 1

gene_index = transcript_index_dictionary['Gene']
exon_start_index = transcript_index_dictionary['Exon_Start']
exon_end_index = transcript_index_dictionary['Exon_End']

##
## Generate an exons dictionary where keys = gene, and values are list of lists
## containing the start and end positions of all of the exons for the key gene.
##
## ie:
##
## key = 'PSTPIP1', value = [ [exon_1_start , exon_1_end] , [exon_2_start , exon_2_end] , [Etc..] ]
##
##

transcript_exons_DICTIONARY = {}
for x in transcript_vars[1:]:

	current_gene = x[ gene_index ]
	current_START = int(x[ exon_start_index ])
	current_END = int(x[ exon_end_index ])

	if current_gene in transcript_exons_DICTIONARY.keys():
		transcript_exons_DICTIONARY[current_gene].append( [current_START, current_END] )

	else:
		transcript_exons_DICTIONARY[current_gene] = []
		transcript_exons_DICTIONARY[current_gene].append( [current_START, current_END] )

##
## create index dictionary to identify index position for fields of interest
##

cnv_index_dictionary = {}
index_tracker = 0
for x in cnv_vars[0]:
	cnv_index_dictionary[x] = index_tracker
	index_tracker += 1

cnv_gene_index = cnv_index_dictionary['GENE']
cnv_START_index = cnv_index_dictionary['start.pos']
cnv_END_index = cnv_index_dictionary['end.pos']

##
## Identify CNV that overlap exons
##

exon_overlapping_CNV = []
exon_overlapping_CNV.append(cnv_vars[0])

for cnv in cnv_vars[1:]:

	#
	# Identify range of genomic positions that the CNV spans
	#

	cnv_gene = cnv[ cnv_gene_index ]
	cnv_start = int(cnv[ cnv_START_index ])
	cnv_end = int(cnv[ cnv_END_index ])
	cnv_range = range( cnv_start, (cnv_end + 1), 1 )

	# loop through exons for the gene for the capture region of interest

	for exon in transcript_exons_DICTIONARY[ cnv_gene ]:
		# Identify genomic positions that the current exon spans
		temp_exon_range = range( exon[0], (exon[1] + 1), 1)

		# Determine if there is an overlap between CNV-range and exon-range

		if len( set(cnv_range).intersection(set(temp_exon_range))) > 0:
			exon_overlapping_CNV.append( cnv )
			break
		else:
			None

for x in exon_overlapping_CNV:
	print x

outPut_NODATE( exon_overlapping_CNV , '%s/%s' % ( cnv_results_dir , 'manual_correction_EXON_OVERLAPPING__ALL_BATCHES_COMPILED_cnv_distributions.txt') )
