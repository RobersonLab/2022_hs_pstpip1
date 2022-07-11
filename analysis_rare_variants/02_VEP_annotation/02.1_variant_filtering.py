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
##  02.1_VARIANT_FILTERING.PY
##
##  All of the previous analysis steps have included the genetic variants for the entire captured
##  cohort.
##
##  However not all of the individuals that were captured should be in the final variant burden analysis.
##  Some samples are non-affected family members from Peter's cohort. Some samples are non-affected
##  controls from the Kappenberger cohort. For the purposes of reporting variant totals counts for missense,
##  synonymous, etc., we need to present only the variants from individuals included in the burdens analysis.
##  The individuals who were included were affected probands from Peter's cohort, all samples from Irish cohort,
##  and the HS affected individuals from the Kaffenberger cohort.
##
##  This code just takes the annotated variants generated from 02.1_VARIANT_EVALUATION and filters the lists
##  to keep only the variants that were found in the individuals that will be analyzed for burden analysis.
##
##
##

def filter_out_nonHS_patients( variant_list, cohortSet ):
	passed_variants = []
	filtered_variants = []

	for x in variant_list:
		affected_individuals = set(x.showAltIndividuals())

		if len( affected_individuals.intersection( cohortSet ) ) > 0:
			passed_variants.append(x)
		else:
			filtered_variants.append(x)


	PASS_SET = VARIANT_SET( passed_variants )
	PASS_RETURN = PASS_SET.returnSet()

	FILTER_SET = VARIANT_SET( filtered_variants )
	FILTER_RETURN = FILTER_SET.returnSet()

	return ( PASS_RETURN , FILTER_RETURN )

#
# Read in file list of analysis cohort.
#

info_dir = '../../info'
output_samples_info_dir = '../../output/samples_info_output'
vep_variant_dir = '../../output/VEP_variant_files'
sites_dir = '../../output/sites_files'

cohort_file = '%s/%s' % ( output_samples_info_dir , 'all_cohorts_HS_AFFECTED_captured_samples.txt' )
cohort_vars = extractVariants(cohort_file)
cohort_dict = {}
for x in cohort_vars[1:]:
	cohort_dict[x[0]] = 1

cohort_set = set(cohort_dict.keys())

allQC_filtered_sites_file = '%s/%s' % ( sites_dir , 'CAPTURE_QC_passed_sites_gnoAAF_appended.txt' )
allQC_vars = extractVariants( allQC_filtered_sites_file )
allQC_VARIANTS = [VARIANT( x , allQC_vars[0] ) for x in allQC_vars[1:]]

allQC_PASSED_SITES , allQC_FILTERED_SITES = filter_out_nonHS_patients( allQC_VARIANTS , cohort_set )

outPut_NODATE( allQC_PASSED_SITES , '%s/%s' % ( sites_dir , 'HS_CAPTURE_allQC_variants_capturedHSpatients_only.txt') )
