#!/usr/bin/env python

##########
# Import #
##########
import numpy
import sys

from burdens_code import *

##############################################################################
##############################################################################
###
###
###    S.5_READ_ANALYSIS_COMPILE.PY
###
###    This script takes in the S.4 generated file with total read counts, on-target
###    read counts, pcr duplicate percentage, etc. and filters the data to only keep
###    individuals that were included in the HS analysis. It also adds information about
###    which capture batch and which sequencing batch the sample was processed in.
###
###

info_dir = '../../info'
coverages_dir = '../../output/capture_coverages'
samples_info_dir = '../../output/samples_info_output'

capture_reads_file = '%s/%s' % ( coverages_dir , 'HS_capture_on_target_and_pcr_duplicates_ALL_SAMPLES.txt' )
capture_reads_vars = extractVariants( capture_reads_file )

sample_info_file = '%s/%s' % ( samples_info_dir , 'all_cohorts_HS_AFFECTED_captured_samples.txt' )
sample_info_vars = extractVariants( sample_info_file )

capture_batch_info_file = '%s/%s' % ( info_dir , 'targeted_capture_samples_batch_info.txt' )
capture_batch_info_vars = extractVariants( capture_batch_info_file )

##
## Establish list of samples that are HS affected and were in the capture.
##

captured_HS_samples_list = []
for x in sample_info_vars[1:]:
	sample_id = x[0]
	captured_HS_samples_list.append( sample_id )

##
## Establish sample_batch_dictionary, where key is sample_id, and values
## are the capture and sequencing batch information for the sample.
##

sample_batch_dictionary = {}
for x in capture_batch_info_vars[1:]:
	sample_id = x[0]
	capture_batch = x[1]
	sequence_batch = x[2]

	if sample_id in captured_HS_samples_list:
		sample_batch_dictionary[ sample_id ] = {}
		sample_batch_dictionary[ sample_id ]['CAPTURE_BATCH'] = capture_batch
		sample_batch_dictionary[ sample_id ]['SEQUENCE_BATCH'] = sequence_batch
	else:
		None

##
## Filter on_target file and add in information about capture/sequence batch.
##
out_list = []

new_header = capture_reads_vars[0]
new_header.append( 'CAPTURE_BATCH' )
new_header.append( 'SEQUENCE_BATCH' )

out_list.append( new_header )

for x in capture_reads_vars[1:]:
	sample_id = x[0]

	if sample_id in captured_HS_samples_list:

		temp_entry = x
		temp_entry.append( sample_batch_dictionary[ sample_id ]['CAPTURE_BATCH'] )
		temp_entry.append( sample_batch_dictionary[ sample_id ]['SEQUENCE_BATCH'] )
		out_list.append( temp_entry )
	else:
		None

outPut_NODATE( out_list , '%s/%s' % ( coverages_dir , 'HS_capture_on_target_and_pcr_duplicates_HS_CAPTURED_ONLY.txt' ) )
