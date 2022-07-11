#!/usr/bin/env python

##########
# Import #
##########
import numpy
import sys

from burdens_code import *

#############################################################################################
#############################################################################################
###
###    S.2_CAPTURE_COVERAGES_COMPILE.PY
###
###    This script compiles all of the individual coverage files generated by S.1 into one large
###    cohort coverages file.
###
###
###

info_dir = '../../info'
coverages_dir = '../../output/capture_coverages'
samples_info_dir = '../../output/samples_info_output'

bed_FILE = '%s/%s' % ( info_dir , 'HS_capture_fosmid_regions_no_overlaps.bed' )
bed_VARS = extractVariants( bed_FILE )

samples_file = '%s/%s' % ( samples_info_dir , 'all_cohorts_HS_AFFECTED_captured_samples.txt' )
samples_VARS = extractVariants( samples_file )

####
####Create 'big dictionary', a dictionary with an entry for every position in the targeted capture;
####

BIG_DICTIONARY = {}
ALL_BD_entries = []
fosmid_LIST = []

for x in bed_VARS:
	chrom = x[0]
	region_start_position = int(x[1])
	region_end_position = int(x[2])

	for y in xrange( int(region_start_position) , int(region_end_position) + 1 ):

		temp_entry = '%s_%s' % (chrom,str(y))

		if temp_entry in BIG_DICTIONARY:
			None
		else:
			BIG_DICTIONARY[temp_entry]={}
			ALL_BD_entries.append( temp_entry )

####
###
### Establish list of sample IDs and sample file names
###
samples_coverage_file_dictionary = {}

for x in samples_VARS[1:]:
	sample_id = x[0]
	temp_coverage_file = '%s/%s_NoDuplicates_Coverages.txt' % ( coverages_dir , sample_id )
	samples_coverage_file_dictionary[ sample_id ] = temp_coverage_file

for x in samples_coverage_file_dictionary.keys():
	print '%s\t%s' % ( x , samples_coverage_file_dictionary[ x ] )

print 'number of samples %i' % ( len( samples_coverage_file_dictionary.keys() ) )

####
####Fill 'big dictionary', each entry will be another dictionary where the keys are the targeted capture sample names
####and the entries in the dictionary are the coverage depth for the corresponding genomic position
####

for x in samples_coverage_file_dictionary.keys():
	sample_name = x

	temp_coverage_file = samples_coverage_file_dictionary[ sample_name ]
	temp_VARS = extractVariants(temp_coverage_file)

	for y in temp_VARS:
		chrom = y[0]
		position = y[1]
		depth = int(y[2])
		chrm_pos = '%s_%s' % (chrom,position)
		try:
			BIG_DICTIONARY[chrm_pos][sample_name] = depth
		except:
			print chrm_pos
			print 'oops dictionary is not filling properly'
			sys.exit(1)

#print 'big dictionary is filled'
####
#### Output 'big dictionary' data into a list format;  each entry contains ['chrom','position','sample_a_depth','sample_b_depth','sample_c_depth',...]
####

final_output_list = []
final_output_list.append(['CHROM','POSITION'])

for x in samples_coverage_file_dictionary.keys():
	sample_name = x
	final_output_list[0].append( sample_name )

unordered_out_list = []

for x in ALL_BD_entries:
	temp_entry = x.split('_')

	for y in final_output_list[0][2:]:
		sample_name = y

		if sample_name in BIG_DICTIONARY[ x ]:
			temp_entry.append( BIG_DICTIONARY[ x ][ sample_name ])
		else:
			temp_entry.append(0)

	unordered_out_list.append(temp_entry)

unordered_out_list.sort( key=lambda z: ( int(z[0]) , int(z[1]) ) )

for x in unordered_out_list:
	final_output_list.append( x )

outPut_NODATE( final_output_list , '%s/%s' % ( coverages_dir , 'HS_capture_all_AFFECTED_samples_region_coverages.txt' ) )
