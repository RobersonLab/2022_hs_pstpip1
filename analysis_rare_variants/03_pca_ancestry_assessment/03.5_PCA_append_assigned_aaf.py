#!/usr/bin/env python

##########
# Import #
##########
import numpy
import sys

from burdens_code import *

###################################################################
###################################################################
##
##  03.5_PCA_APPEND_ASSIGNED_AAF.PY
##
##  This script just adds the gnomAD ancestry assigned to each HS cohort individual, by
##  the kNN PCA analysis, to the HS-affected-captured-sample file list that will be used as
##  a reference for which samples to include in the burdens analysis.
##

pca_dir = '../../output/pca_results'
sample_info_dir = '../../output/samples_info_output'

###
### generate PCA table
###

pca_results_file = '%s/%s' % ( pca_dir , 'HS_capture_cohort_PCA_assignment_knn25.txt' )
pca_results_vars = extractVariants( pca_results_file )

samples_info_file = '%s/%s' % ( sample_info_dir , 'all_cohorts_HS_AFFECTED_captured_samples.txt')
samples_info_vars = extractVariants( samples_info_file )

pca_dictionary = {}
for x in pca_results_vars[1:]:
	sample = x[0]
	pca_assignment = x[-1]
	pca_dictionary[ sample ] = pca_assignment

samples_info_vars[0].append( 'PCA_ASSIGNMENT' )

for x in samples_info_vars[1:]:
	sample = x[0]
	x.append( pca_dictionary[sample] )

outPut_NODATE( samples_info_vars , '%s/%s' % ( sample_info_dir , 'all_cohorts_HS_AFFECTED_captured_samples_withPCAinfo.txt') )
