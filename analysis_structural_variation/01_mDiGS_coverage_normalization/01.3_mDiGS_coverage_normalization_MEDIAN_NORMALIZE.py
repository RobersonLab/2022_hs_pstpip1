#!/usr/bin/env python

##########
# Import #
##########
import argparse
import numpy
import sys

from burdens_code import *

def batch_normalize_coverages( full_coverage_list , batch_index_list , bed_vars , WINDOW_SIZE , batch_group ):
	##
	## Convert the coverage list of lists values into floats.  They are read in as strings, but we will
	## need floats in order to perform our normalization calculations.
	##

	covg_float = []

	for line in full_coverage_list[1:]:

		temp_entry = []
		for value in line:
			temp_entry.append( float(value) )

		covg_float.append( temp_entry)

	full_coverage_array = numpy.array(covg_float)

	##
	## Create sliced array containing just the coverage info for samples in batch being analyzed
	##

	covg_array = full_coverage_array[:,batch_index_list]

	fosmid_normalized_coverage_list = []

	##
	## Generate sliding window data
	##

	bed_index_dictionary = {}
	bed_index_dictionary['CHROM'] = 0
	bed_index_dictionary['START'] = 1
	bed_index_dictionary['END'] = 2
	bed_index_dictionary['GENE'] = 3

	for region in bed_vars:
		bed_chrom = float( region[bed_index_dictionary['CHROM']] )
		bed_start = float( region[bed_index_dictionary['START']] )
		bed_end = float( region[bed_index_dictionary['END']] )
		bed_gene = region[bed_index_dictionary['GENE']]

		window_list = range( int(bed_start), (int(bed_end) + 1) , WINDOW_SIZE )

		## Generate sliding windows over which to average sequencing coverage; if the last sliding window (remainder) is less than half the size of the chosen
		## sliding window, then the remainder window is merged into the last sliding window to create a slightly larger window.
		##
		## ie, if sliding window is 100bp, and the capture region results in 3000 100bp sliding windows plus one 10bp sliding window, then the code will analyze the
		## data using 2999 100bp sliding windows and 1 110bp sliding window instead of 3000 100 bp sliding windows and 1 10bp sliding window.

		if (int(bed_end) - int(bed_start)) % WINDOW_SIZE < ( WINDOW_SIZE/2.0 ) :
			window_list = window_list[:-1]
		else:
			None

		for window in window_list:

			reported_window = ( float(window) + (WINDOW_SIZE/2) )
			# This isn't necessary but was a choice that I made to report the genomic position for each sliding window
			# as the middle position of the sliding window, instead of the first position in the sliding window.
			temp_entry = []
			temp_entry.append( bed_chrom )
			temp_entry.append( reported_window )

			#
			# Extract the data for just the sliding window we are currently working with.
			#

			if window == window_list[-1]:
				window_array = covg_array[(covg_array[:,0] == bed_chrom) & (covg_array[:,1] >= window) & (covg_array[:,1] <= bed_end)][:,2:]
				print '%s:\tlast window length:\t%i' % ( bed_gene , len(window_array[:,0]))
			else:
				window_array = covg_array[(covg_array[:,0] == bed_chrom) & (covg_array[:,1] >= window) & (covg_array[:,1] <= (window+WINDOW_SIZE-1))][:,2:]

			#
			# Convert the window data frame into one value for each sample: the average coverage for the sample over the region.
			# In other words, [column_sums]/[number_of_columns]
			#
			mean_window_coverage = ( window_array.sum(axis=0) ) / len(window_array[:,0])

			#temp_colsums = list(window_array.sum(axis=0))

			temp_entry += list( mean_window_coverage )

			#
			#  apppend sliding window to new list, and ultimately convert that list of lists into a new numpy array
			#

			fosmid_normalized_coverage_list.append( temp_entry )

	norm_fos_array = numpy.array( fosmid_normalized_coverage_list )

	##
	## Normalize across samples within fosmid sliding windows
	##

	for position in norm_fos_array:

		# use [2:] here because first column is 'chrom', second column is 'position', and subsequent columns are sample coverages
		if position[2:].mean() == 0.0:
			position[2:] = position[2:]

		else:

			median_list = list(position[2:])
			median_list.sort()

			median = median_list[ (len(median_list)/2) ]

			if median == 0:
				position[2:] = position[2:]

			else:
				position[2:] = ( position[2:] / median )

	##
	## Normalize coverages within each individual across entire capture region
	##
	columns_number = len( norm_fos_array[0,:] )

	for col in range( 2 , (columns_number) ):

		norm_fos_array[:,col] = 2*( norm_fos_array[:,col] / norm_fos_array[:,col].mean() )

		###
		### NOTE:  This sample normalization step also multiplies everything by 2 to convert the normalized
		### coverages back to normalized estimated copy-number for each sliding window.
		###

	for region in bed_vars:
		bed_chrom = float( region[bed_index_dictionary['CHROM']] )
		bed_start = float( region[bed_index_dictionary['START']] )
		bed_end = float( region[bed_index_dictionary['END']] )
		bed_gene = region[bed_index_dictionary['GENE']]

		temp_sliced_array = norm_fos_array[ (norm_fos_array[:,0] == bed_chrom) & (norm_fos_array[:,1] >= bed_start) & (norm_fos_array[:,1] <= bed_end) ]

		header = []
		for index in batch_index_list:
			header.append( full_coverage_list[0][index])

		temp_out_list = []
		temp_out_list.append( header )

		for each_position in temp_sliced_array:
			temp_pos_list = list( each_position )
			temp_pos_list[0] = int( temp_pos_list[0] )
			temp_pos_list[1] = int( temp_pos_list[1] )

			temp_out_list.append( temp_pos_list )


		temp_file_name = '%s_%s_batch_normalized_coverages.txt' % ( batch_group , bed_gene )
		temp_out_file = '%s/%s' % ( norm_coverage_out_dir , temp_file_name )

		outPut_NODATE( temp_out_list , temp_out_file)

#############################################################################################
#############################################################################################
###
###    01.3_mDiGS_COVERAGE_NORMALIZATION.PY
###
###    This script normalizes the individual raw coverage files generated by 01.1_capture_create_coverage_files.sh and 01.2_capture_coverages_compile.py
###    using the mDiGS strategy.  The idea behind the mDiGS CNV analysis is that variations
###    in coverages between samples is mainly due to variation introduced as a result of processing samples in different
###    capture batches and different sequencing bacthes.  Because mDiGS multiplexes samples in the same capture batch, which
###    is then all analyzed on the same sequencing run, the bulk of this variation is controlled for.  So sample coverages
###    can be normalized within the cohort, and then the average coverage of the cohort can be used as the reference
###    coverage to establish approximately what 2X diploid sequencing coverage looks like.
###
###    Our samples were captured in 3 batches and sequenced in two batches.  All captures were performed using the same
###    original prep of capture probes.
###
###    For normalization, samples are separated into capture batches, and separate normalizations are performed
###    for each capture batch.
###
###    Normalization procedure:
###         1. Separate coverage data into the discrete genomic regions that were captured
###				- Generate sliding window average coverages for each sample
###
###         2. For each genomic sliding window, normalize each sample's coverages to the mean coverage for that window over the entire cohort
###
###			2. Normalize the coverage for each fosmid-normalized sliding window within an individual to the mean
###            sliding window coverage for all sliding windows within that individual (over all fosmid regions)
###
###    Normalizing both within each fosmid and within each individual should control for unequal loading of fosmids during capture probes
###    generation as well as unequal loading of sample libraries during the capture hybridization.
###
###
###    *** NOTE ***
###    In the last normalization step that normalizes within an individual across all regions of the targete capture, I multiply by 2 to
###    convert the normalized relative coverage into normalized estimated copy-number.
###
###
###

norm_coverage_out_dir = '../../output/MEDIAN_NORMALIZED_coverages'
info_dir = '../../info'
sample_info_dir = '../../output/samples_info_output'
raw_coverages_dir = '../../output/raw_capture_coverages'

coverages_file = '%s/%s' % ( raw_coverages_dir , 'HS_capture_all_AFFECTED_samples_region_coverages.txt' )
bed_file = '%s/%s' % ( info_dir , 'HS_capture_fosmid_regions_no_overlaps_WITH_GENE.bed' )
batch_file = '%s/%s' % ( info_dir , 'targeted_capture_samples_batch_info.txt' )
samples_info_file = '%s/%s' % ( sample_info_dir , 'all_cohorts_HS_AFFECTED_captured_samples_withPCAinfo.txt' )

coverage_vars = extractVariants( coverages_file )
bed_vars = extractVariants( bed_file )
batch_vars = extractVariants( batch_file )
samples_info_vars = extractVariants( samples_info_file )

##
## Create samples list; not all samples in the batch info file are HS affected or included in the burdens analysis
##

samples_list = []
for x in samples_info_vars:
	if x[0] in samples_list:
		None
	elif x[0] != 'SAMPLE_ID':
		samples_list.append( x[0])

##
## generate index dictionary for sample batch_info file
##

index_tracker = 0
seq_batch_dictionary = {}
for x in batch_vars[0]:
	seq_batch_dictionary[x] = index_tracker
	index_tracker += 1

##
## generate sample batch dictionary; key = batch#, value = list of samples in that batch
##

sample_batch_dictionary = {}
for x in batch_vars:
	sample_name = x[seq_batch_dictionary['sample']]
	curr_sample_sequencing_batch = x[seq_batch_dictionary['capture_batch']]
	if sample_name in samples_list:
		if curr_sample_sequencing_batch in sample_batch_dictionary:
			sample_batch_dictionary[curr_sample_sequencing_batch].append( sample_name )
		else:
			sample_batch_dictionary[curr_sample_sequencing_batch] = []
			sample_batch_dictionary[curr_sample_sequencing_batch].append( sample_name )

	else:
		None

##
## generate index dictionary for coverage file header; value = column_name (sample_id) , value = index position of the sample in the header row
##

coverage_index_dictionary = {}
index_tracker = 0
for x in coverage_vars[0]:
	coverage_index_dictionary[x] = index_tracker
	index_tracker += 1

batch_1_index_list = []
batch_1_index_list.append(0)
batch_1_index_list.append(1)
for sample in sample_batch_dictionary['1']:
	sample_index = coverage_index_dictionary[sample]
	batch_1_index_list.append( sample_index )
batch_1_index_list.sort()

print 'Samples in sequencing batch #1:\t%i' % (len(batch_1_index_list))
print batch_1_index_list
for x in batch_1_index_list:
	print coverage_vars[0][x]
print

batch_2_index_list = []
batch_2_index_list.append(0)
batch_2_index_list.append(1)
for sample in sample_batch_dictionary['2']:
	sample_index = coverage_index_dictionary[sample]
	batch_2_index_list.append( sample_index )
batch_2_index_list.sort()

print 'Samples in sequencing batch #2:\t%i' % (len(batch_2_index_list))
print batch_2_index_list
for x in batch_2_index_list:
	print coverage_vars[0][x]
print

batch_3_index_list = []
batch_3_index_list.append(0)
batch_3_index_list.append(1)
for sample in sample_batch_dictionary['3']:
	sample_index = coverage_index_dictionary[sample]
	batch_3_index_list.append( sample_index )
batch_3_index_list.sort()

print 'Samples in sequencing batch #3:\t%i' % (len(batch_3_index_list))
print batch_3_index_list
for x in batch_3_index_list:
	print coverage_vars[0][x]
print

gene_list = []
for x in bed_vars:
	if x[-1] in gene_list:
		None
	else:
		gene_list.append(x[-1])

WINDOW_SIZE = 50

batch_normalize_coverages( coverage_vars , batch_1_index_list , bed_vars , WINDOW_SIZE , 'BATCH_1')

batch_normalize_coverages( coverage_vars , batch_2_index_list , bed_vars , WINDOW_SIZE , 'BATCH_2')

batch_normalize_coverages( coverage_vars , batch_3_index_list , bed_vars , WINDOW_SIZE , 'BATCH_3')
