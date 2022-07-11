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
###    S.6_SINGLETON_ANALYSIS.PY
###
###    This script evaluates the missense singletons in the gnomad data set.  It evaluates
###    the number of variants in the data set that were singletons compared to variants that were
###    observed more frequently in populations.  It also calculates a 'singleton frequency' for
###    each gene in the targetd capture.  This singleton frequency is not 100% accurate, but a best
###    approximation.  Different genomic positions in the gnomad data set have different total
###    observed allele counts.  And even if two sites have the same total allele count, it cannot be
###    determined if the allele counts are coming from an identical pool of individuals between two
###    sites.  Ideally, every position would have identical observed allele numbers from the same
###    pool of individuals.  We chose to calculate singleton frequencies by totaling up the number
###    of singletons observed, and then dividing that by the average total allele count at all of
###    singleton positions.
###

transcript_dir = '../../output/transcript_file'
variant_dir = '../../output/VEP_variant_files'

transcript_file = '%s/%s' % ( transcript_dir , 'HS_capture_primary_transcript_exon_sequences.txt' )
transcript_vars = extractVariants(transcript_file)

REGULATORY_BUFFER_WINDOW = 0
capture_gene_list = ['MAML3','NOTCH4','PSTPIP1','APH1A','APH1B','PSEN1','PSEN2','NCSTN','PSENEN','PSTPIP2','NOTCH1','NOTCH2','NOTCH3','DLL1','DLL3','DLL4','JAG1','JAG2','RBPJ','MAML1','MAML2']

###
### Establish the coding base length for each gene in capture.
###
gene_length_dictionary = {}

for gene_to_analyze in capture_gene_list:

	# Establish transcript dictionaries for current gene being analyzed.

	exon_positions_dictionary,coding_positions_dictionary,gene_position_dictionary, codon_positions_dictionary, strandedness = GenerateTranscriptGenomicPositionDictionaries(transcript_vars,gene_to_analyze,REGULATORY_BUFFER_WINDOW)

	gene_length_dictionary[gene_to_analyze] = len(coding_positions_dictionary)

###
### Set up singleton_dictionary.  This dictionary will keep track of the number of variants in gnomad that are singletons.
###
###

var_file = '%s/%s' % ( variant_dir , 'gnomad_MISSENSE_variants.txt' )
vars = extractVariants(var_file)

gem_vars = [ VARIANT( x , vars[0] ) for x in vars[1:] ]

sub_groups_list = [ 'nfe' , 'afr' , 'amr' , 'sas' , 'eas' ]

singleton_dictionary = {}
singleton_dictionary['all_subsets'] = {}
for gene in capture_gene_list:
	singleton_dictionary['all_subsets'][gene] = {}
	singleton_dictionary['all_subsets'][gene]['singleton_counts'] = 0
	singleton_dictionary['all_subsets'][gene]['singleton_depths'] = []

for each in sub_groups_list:
	singleton_dictionary['aaf_gnomad_%s' % (each) ] = {}

	for gene in capture_gene_list:

		singleton_dictionary['aaf_gnomad_%s' % (each) ][gene] = {}
		singleton_dictionary['aaf_gnomad_%s' % (each) ][gene]['singleton_counts'] = 0
		singleton_dictionary['aaf_gnomad_%s' % (each) ][gene]['singleton_depths'] = []

for temp_var in gem_vars:

	gene = temp_var.info['GENE'].split('/')[0]

	all_AC = 0
	all_AN = 0

	for each in sub_groups_list:
		tmp_AC = int(temp_var.info['AC_%s' % (each)])
		tmp_AN = int(temp_var.info['AN_%s' % (each)])

		all_AC += tmp_AC
		all_AN += tmp_AN

	if all_AC == 1:
		for each in sub_groups_list:

			tmp_AC = int(temp_var.info['AC_%s' % (each)])
			tmp_AN = int(temp_var.info['AN_%s' % (each)])

			if tmp_AC == 1:
				singleton_dictionary['aaf_gnomad_%s' % (each)][gene]['singleton_counts'] += tmp_AC
				singleton_dictionary['aaf_gnomad_%s' % (each)][gene]['singleton_depths'].append(tmp_AN)
			else:
				None

		singleton_dictionary['all_subsets'][gene]['singleton_counts'] += all_AC
		singleton_dictionary['all_subsets'][gene]['singleton_depths'].append( all_AN	)

	else:
		None

###
###
### Calculate overall and population specific singleton frequencies.  Frequencies are still estimates as singleton frequencies
### are going to be a function of the number of alleles sampled.  As different genomic positions have different total allele
### denominators, we calculated the singleton frequencies as the number of singletons observed for the average total allele
### number for all singletons in observed in the gene.
###

out_list = []
header = [ 'POPULATION', 'GENE', 'GENE_LENGTH', 'DEPTH_MEAN', 'DEPTH_STD', 'DEPTH_VARIANCE', 'SINGLETON_COUNT' , 'SINGLETON_PER_ALLELE' ]
out_list.append(header)

sub_groups_list.append( 'all_subsets' )

for gene in capture_gene_list:
	for each in sub_groups_list:
		if each == 'all_subsets':
			number_of_singletons = singleton_dictionary[ 'all_subsets' ][ gene ][ 'singleton_counts' ]
		else:
			number_of_singletons = singleton_dictionary[ 'aaf_gnomad_%s' % ( each )][ gene ][ 'singleton_counts' ]

		print 'singletons:\t%i' % ( number_of_singletons )
		print 'depths:'         #NOTE: 'depth' here was a bad wordage choice.  'depth' really means total number of alleles... ie, what 'depth' of the population do you have to sequence to find a singleton.

		if number_of_singletons > 0:
			depths_array = numpy.zeros(number_of_singletons)

			for x in xrange(0,number_of_singletons):
				if each == 'all_subsets':
					depths_array[x] = singleton_dictionary[ 'all_subsets' ][ gene ][ 'singleton_depths' ][x]
				else:
					depths_array[x] = singleton_dictionary[ 'aaf_gnomad_%s' % ( each )][ gene ][ 'singleton_depths' ][x]

			depths_mean = depths_array.mean()
			depth_std = depths_array.std()
			depth_variance = depths_array.var()
			print 'depth mean:\t', depths_mean
			print 'depth std:\t', depth_std
			print 'std/depth pct:\t', (depth_std/depths_mean)*100
			print 'POP_SPECIFIC_SINGLETON_RATE:\t', number_of_singletons/depths_mean
			print

			tmp_entry = []

			tmp_entry.append(each)
			tmp_entry.append(gene)
			tmp_entry.append(gene_length_dictionary[gene])
			tmp_entry.append('%.2f' % (depths_mean))
			tmp_entry.append('%.2f' % (depth_std))
			tmp_entry.append('%.2f' % (depth_variance))
			tmp_entry.append(number_of_singletons)
			tmp_entry.append('%.8f' % (number_of_singletons/depths_mean) )
			out_list.append(tmp_entry)

		else:
			tmp_entry = []
			tmp_entry.append(each)
			tmp_entry.append(gene)
			tmp_entry.append(gene_length_dictionary[gene])
			tmp_entry.append('0')
			tmp_entry.append('0')
			tmp_entry.append('0')
			tmp_entry.append('0')
			tmp_entry.append('0')
			out_list.append(tmp_entry)
outPut_NODATE( out_list , '%s/%s' % ( variant_dir , 'singleton_variant_summary.txt' ) )

###
###  gnomad variant singleton histogram.  Loop through variants and determine which bin a variant falls into
###  with regards to how many alt alleles were observed for that variant in all of the populations.
###

var_histo_dict = {}

var_histo_dict['1'] = 0
var_histo_dict['2-10'] = 0
var_histo_dict['>10'] = 0
var_histo_dict['total'] = 0

for x in gem_vars:
	AC_count = 0

	for group in sub_groups_list:
		if group == 'all_subsets':
			None
		else:
			AC_count += int( x.info[ 'AC_%s' % (group) ] )

	if AC_count == 1:
		var_histo_dict['1'] += 1
		var_histo_dict['total'] += 1
	elif AC_count >= 2 and AC_count <= 10:
		var_histo_dict['2-10'] += 1
		var_histo_dict['total'] += 1
	elif AC_count > 10:
		var_histo_dict['>10'] += 1
		var_histo_dict['total'] += 1

histo_out_list = []
histo_header = [ 'AC_RANGE' , 'COUNTS' ]

histo_out_list.append( histo_header )

for range in var_histo_dict.keys():
	histo_out_list.append( [ range , var_histo_dict[ range ] ] )

outPut_NODATE( histo_out_list , '%s/%s' % ( variant_dir , 'gnomad_singletons_histogram.txt' ) )
