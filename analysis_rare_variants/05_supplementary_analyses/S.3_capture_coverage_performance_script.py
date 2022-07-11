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
###    S.3_CAPTURE_COVERAGE_PERFORMANCE_SCRIPT.PY
###
###    This script parses the HS cohort compiled coverages file.  It determines mean coverages for
###    each individual over each fosmid used in the targeted capture.
###
###
###
###
###

info_dir = '../../info'
coverages_dir = '../../output/capture_coverages'
samples_info_dir = '../../output/samples_info_output'
results_dir = '../../results'

capture_batch_file = '%s/%s' % ( info_dir , 'targeted_capture_samples_batch_info.txt' )
capture_batch_vars = extractVariants(capture_batch_file)

capture_batch_dictionary = {}
index_tracker = 0
for x in capture_batch_vars[1:]:
	capture_batch_dictionary[x[0]] = {}
	capture_batch_dictionary[x[0]]['capture_batch'] = x[1]
	capture_batch_dictionary[x[0]]['sequencing_batch'] = x[2]


cohort_file = '%s/%s' % ( samples_info_dir , 'all_cohorts_HS_AFFECTED_captured_samples.txt' )
cohort_vars = extractVariants( cohort_file )
cohort_sample_list = [ x[0] for x in cohort_vars[1:] ]

capture_coverage_file = '%s/%s' % ( coverages_dir , 'HS_capture_all_AFFECTED_samples_region_coverages.txt' )
coverage_vars = extractVariants(capture_coverage_file)

##
## Set up quick header dictionary to get index position of each sample name in the header line.
## Later will be able to call LINE[ sample_index_dict["sample_id"] ] to hash access coverage for
## the individual of interest at the genomic position... genomic position == what line you are on in file.
##
sample_index_dict = {}
tracker = 0
for x in coverage_vars[0]:
	sample_index_dict[x] = tracker
	tracker += 1

bed_file = '%s/%s' % ( info_dir , 'HS_fosmid_set.bed' )
bed_vars = extractVariants( bed_file )

##
## Set up fosmid dictionary. Each key is a fosmid name, and the corresponding value is sub-dictionary.
## The keys for the sub-dictionary within each fosmid key are the sample_ids, and the value will be the
## a list of coverages for the individual for the fosmid that the sub-dictionary is under.  This will
## Then allow us to calculate an average coverage per individual per fosmid.
##

fosmid_information_dictionary = {}
fosmids_dict = {}
out_dictionary = {}

for x in bed_vars:
	chrom = x[0]
	fosmid_name = x[-1]
	start = int(x[1])
	end = int(x[2])

	fosmid_information_dictionary[ fosmid_name ] = {}
	fosmid_information_dictionary[ fosmid_name ][ 'chrom' ] = chrom
	fosmid_information_dictionary[ fosmid_name ][ 'start' ] = start
	fosmid_information_dictionary[ fosmid_name ][ 'end' ] = end

	fosmids_dict[ fosmid_name ] = {}
	out_dictionary[ fosmid_name ] = {}

	for individual in cohort_sample_list:

		fosmids_dict[fosmid_name][individual] = []

##
## Range dictionary:  The fosmids used to capture a discrete genomic region were designed
## to overlap by ~5000 bp at their ends to ensure that seamless coverage of the genomic
## region were obtained.  As a result, there will be particular genomic positions that
## fall within 2 fosmids instead of just one fosmid.  We need to identify those postions in order
## to attribute the coverages at those positions to both fosmids when we are calculating
## the coverage for each fosmid.
##
## The range_dictionary just has chrm_positions for keys, and the value is a list of fosmids
## that the position falls within.
##

print 'making range_dictionary'
print
range_dictionary = {}
for x in fosmid_information_dictionary.keys():
	fosmid_name = x
	chrom = fosmid_information_dictionary[ x ][ 'chrom' ]
	start = fosmid_information_dictionary[ x ][ 'start' ]
	end = fosmid_information_dictionary[ x ][ 'end' ]

	for each in range( start , end + 1 ):
		position = '%s_%i' % ( chrom , each )

		if position in range_dictionary:
			range_dictionary[position].append(fosmid_name)
		else:
			range_dictionary[position]= []
			range_dictionary[position].append(fosmid_name)

print 'filling fosmid_dictionary'
print

for x in coverage_vars[1:]:
	# identify the chrom_position
	chrom = x[ sample_index_dict['CHROM'] ]
	pos = x[ sample_index_dict['POSITION'] ]
	chrm_pos = '%s_%s' % (chrom , pos)

	# identify which fosmids the chrom_position falls under
	if chrm_pos in range_dictionary:

		# for each fosmid that the chrom_position is in, cycle through and add the coverage
			# for each individual at that position to the fosmid_dict
		for fosmid in range_dictionary[ chrm_pos ]:
			for sample in cohort_sample_list:
				fosmids_dict[fosmid][sample].append( int( x[ sample_index_dict[ sample ]] ) )
	else:
		print 'ooops... this chrm_position is not in the range dictionary for some reason'
		None

##
## sample_coverage_dictionary is a dictionary of dictionaries. Highest level keys is
## sample_ids, then there are fosmid_name keys for each sample_id dictionary. Values
## for the fosmid_name keys are the mean coverage for the sample_id over the entire fosmid.
##

print 'starting fosmid evaluation output'
print
sample_coverage_dictionary = {}
for x in cohort_sample_list:
	sample_coverage_dictionary[x] = {}

for fosmid in fosmids_dict.keys():
	####
	#### CHECK
	#### There should be as many values in the MEAN COVERAGE calculation as there are positions
	#### in the fosmid.  Ensure that the length of the fosmid region is equal to the number of
	#### values that are going into the MEAN COVERAGE calculation.

	fos_start = fosmid_information_dictionary[ fosmid ][ 'start' ]
	fos_end = fosmid_information_dictionary[ fosmid ][ 'end' ]
	fos_length = ( int(fos_end) - int(fos_start) + 1 )

	first_coverage_entry = True
	for each in fosmids_dict[ fosmid ].keys():
		if first_coverage_entry:
			coverages_entries = len( fosmids_dict[ fosmid ][each])
			first_coverage_entry = False
		else:
			if len( fosmids_dict[ fosmid ][each]) == coverages_entries:
				None
			else:
				sys.exit( 'oops, different samples have different numbers of coverage entries for fosmid:\t%s' % (fosmid) )

	if coverages_entries == fos_length:
		None
		#print 'EVERYTHING LOOKS GOOD !!'
	else:
		sys.exit('Something is wrong. Fosmid:\t%s,\tExpected_length:\t%i,\tCoverage_entries:\t%i' % ( fosmid , fos_length , coverages_entries ) )

	for sample in fosmids_dict[fosmid].keys():

		temp_array = numpy.zeros(len(fosmids_dict[fosmid][sample]))
		counter = 0
		for position in fosmids_dict[fosmid][sample]:
			temp_array[counter] = position
			counter += 1

		sample_coverage_dictionary[sample][fosmid] = temp_array.mean()
		out_dictionary[fosmid][sample] = temp_array.mean()

##
## Parse and output mean coverage data.
##

out_list = []
out_header = ['INDIVIDUAL']

for x in fosmids_dict.keys():
	out_header.append(x)

temp_coverage_list = []

for x in cohort_sample_list:
	temp_entry = []
	temp_entry.append(x)
	for fosmid in out_header[1:]:
		temp_entry.append(sample_coverage_dictionary[x][fosmid])

	temp_entry.append(capture_batch_dictionary[x]['capture_batch'])
	temp_entry.append(capture_batch_dictionary[x]['sequencing_batch'])

	temp_coverage_list.append(temp_entry)

out_header.append('capture_batch')
out_header.append('sequencing_batch')

out_list.append(out_header)
for x in temp_coverage_list:
	out_list.append(x)

outPut_NODATE( out_list , '%s/%s' % ( coverages_dir , 'HS_capture_mean_fosmid_coverages.txt' ) )
