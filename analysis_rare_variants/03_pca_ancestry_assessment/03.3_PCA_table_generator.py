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
##  03.3_PCA_TABLE_GENERATOR.PY
##
##  This script creates a genotype table to be used for geographic ancestry principal component analysis.  Genotypes for individuals
##  for each variant in a variant list will be coded as 0 (ref/ref), 1 (het), or 2 (alt/alt).  In addition to recoding the HS patient
##  genotype calls, the script simulates genotypes for individuals for the 5 subpopulation groups pulled from gnomad (nfe, afr, amr, eas, sas).
##  For a simulated individual, for each variant a random number between 0 and 1 is selected. That number is compared to the geographic ancestry
##  aaf for the individual being simulated.  If the random number is less than or equal to the ancestry aaf, the individual accrues an alt allele; if it is greater
##  than the ancestry aaf, the individual accrues a reference allele. Ultimately the simulated individuals will be used in the PCA analysis with
##  the HS patient samples in order to determine what geographic ancestry group the HS patients look most similar to.
##
##  This script generates one 0/1/2 coded genotype file with 3 categories of data.
##         - HS capture genotype data
##         - simulated genotype data to be used for training kNN analysis to categorize HS patients. these data are labeled as 'training'.
##         - simulated genotype data used to validate kNN algorithm and optimize K value. these data are labeled as 'validation'. they will
##           be used as a test case for the kNN algorithm. we will use kNN to categorize the validation genotypes using different K values.
##           The K value that produces the fewest ancestry misclassifications of the validation dataset will then be used for categorizing the HS cohort.
##

pca_dir = '../../output/pca_results'
sample_info_dir = '../../output/samples_info_output'

###
### generate PCA table
###
capture_file = '%s/%s' % ( pca_dir , 'PCA_capture_PASSED_variants_capturedHSpatients_only__gnomad_aaf_appended.txt' )
capture_vars=extractVariants(capture_file)

AAF_COMMON_THRESHHOLD_FILTER = 0.0

aaf_filtered_variant_list = []
for x in capture_vars[1:]:
	temp = VARIANT(x,capture_vars[0])

	if float(temp.info['aaf_gnomad_all_subsets']) > AAF_COMMON_THRESHHOLD_FILTER:
		aaf_filtered_variant_list.append(temp)
		#temp.printVariantHeaderInfo()
	else:
		None
		#temp.printVariantHeaderInfo()

variant_set = VARIANT_SET(aaf_filtered_variant_list)

simulated_training_genotypes = 1000
simulated_validation_genotypes = 1000

principal_component_table = variant_set.prCompTable_w_simulations(simulated_training_genotypes, simulated_validation_genotypes)

###
### Filter non-HS-captured individuals out of PCA table
###

cohort_file = '%s/%s' % ( sample_info_dir , 'all_cohorts_HS_AFFECTED_captured_samples.txt' )
cohort_vars = extractVariants( cohort_file )
cohort_list = [ x[0] for x in cohort_vars[1:] ]

filtered_PCA_table = []
filtered_PCA_table.append( principal_component_table[0] )
for x in principal_component_table[1:]:
	if 'training' in x[0]:
		filtered_PCA_table.append( x )
	elif 'validation' in x[0]:
		filtered_PCA_table.append( x )
	elif x[0] in cohort_list:
		filtered_PCA_table.append( x )
	else:
		None

outPut_NODATE( filtered_PCA_table , '%s/%s' % ( pca_dir , 'PCA_%i-%i_simulated-validation_and_HS_captured_pca_table.txt' % (simulated_training_genotypes , simulated_validation_genotypes) ) )
