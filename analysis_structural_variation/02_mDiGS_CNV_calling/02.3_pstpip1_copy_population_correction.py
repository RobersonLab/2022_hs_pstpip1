#!/usr/bin/env python

##########
# Import #
##########
import argparse
import numpy
import sys

from burdens_code import *


#############################################################################################
#############################################################################################
###
###    02.3_pstpip1_copy_population_correction.py
###
###    Because the PSTPIP1 deletion CNV polymorphism was so prevalent in our cohort (>50%),
###    the mean coverage normalization process is dramatically skewed as the average coverage is closer to 1-copy than 2-copies. This
###	   skewing is so extreme that median-normalization of the CNV calls does not correct for the skewing, and in Batches 1 and 3 the only
###	   copy number populations called in the cohort were 0,2, and 4.  This CNV is called in the gnomAD structural database as a deletion,
###	   and not as an amplification.  Furthermore, it is highly unlikely that our cohort would have large populations of homozygous deletions
###	   and homozygous amplifications without observing one individual who was heterozygous for at deletion or insertion.  Rather what makes
###	   more sense is that as a result of the high prevalence of the CNV, the copy-number data has been normalized to the 1-copy population,
###	   making 1-copy individuals be reported as 2-copy individuals... and likewise, actual 2-copy individuals (who have 2X coverage as the
###	   1-copy individuals) are reporting as 4-copy because the 1-copies are reporting as 2-copies.  Thus I believe that the copy number calls
###	   need to be corrected from 0,2,4 to 0,1,2 respectively.
###
###	   I would have loved to find a way to correct for high prevalent polymorphic CNV calls in the normalization or filtering processes, but
###    without an external 2-copy reference, I don't know if it is possible to identify the true 2-copy population other than the logic were
###    we have taken here. This python script is essentially equivalent to a manual correction of these values. It isn't the most elegant solution.
###    But I strongly believe that the 0,2,4 copy population calls is a mis-interpretation of the data and needs to be corrected.
###
###
###

info_dir = '../../info'
cnv_dir = '../../output/mDiGS_CNV'

batch_1_cnv_file = '%s/%s' % ( cnv_dir , 'BATCH_1_ALL_GENES_COMPILED_filtered_SMOOTHED_CNV_copy_number_calls_TWO_STAGE_CNV_FILTER.txt' )
batch_3_cnv_file = '%s/%s' % ( cnv_dir , 'BATCH_3_ALL_GENES_COMPILED_filtered_SMOOTHED_CNV_copy_number_calls_TWO_STAGE_CNV_FILTER.txt' )

batch_1_vars = extractVariants( batch_1_cnv_file )
batch_3_vars = extractVariants( batch_3_cnv_file )

def convert_skewed_CN_populations( copy_number_vars ):

	converted_out_list = []
	converted_out_list.append( copy_number_vars[0] )

	info_field = [ 'chrom' , 'arm' , 'start.pos' , 'end.pos' , 'n.probes' , 'GENE' ]
	SAMPLES_FOUND = False

	first_sample_index = 'not_found'

	header_dictionary = {}
	index_tracker = 0
	for x in copy_number_vars[0]:
		header_dictionary[ x ] = index_tracker

		if SAMPLES_FOUND:
			None

		else:
			if x in info_field:
				None
			else:
				first_sample_index = index_tracker
				SAMPLES_FOUND = True

		index_tracker += 1

	for cnv in copy_number_vars[1:]:

		if cnv[ header_dictionary['GENE'] ] == 'PSTPIP1':

			temp_entry = cnv[0:first_sample_index]

			for sample in cnv[ first_sample_index: ]:
				if int(sample) == 0:
					None
				elif int(sample) == 2:
					sample = '1'
				elif int(sample) == 4:
					sample = '2'
				else:
					print 'oops, a CNV copy number call is something other than [0,2,4] '
					break

				temp_entry.append( sample)

			converted_out_list.append( temp_entry )

		else:
			converted_out_list.append(cnv)

	return(converted_out_list)

converted_batch_1 = convert_skewed_CN_populations( batch_1_vars )
converted_batch_3 = convert_skewed_CN_populations( batch_3_vars )

outPut_NODATE( converted_batch_1 , '%s/%s' % ( cnv_dir , 'manual_cnv_correction_BATCH_1_ALL_GENES_COMPILED_filtered_SMOOTHED_CNV_copy_number_calls_TWO_STAGE_CNV_FILTER.txt'))
outPut_NODATE( converted_batch_3 , '%s/%s' % ( cnv_dir , 'manual_cnv_correction_BATCH_3_ALL_GENES_COMPILED_filtered_SMOOTHED_CNV_copy_number_calls_TWO_STAGE_CNV_FILTER.txt'))
