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
##  03.1_PCA_ANALYSIS_capture_exome_QC.pyplot
##
##  This script removes positions/sites from the gnomAD EXOME,and HS_CAPTURE
##  data sets based on:
##        - coverage depth at the site
##        - number of site alleles called as unknown for HS CAPTURE data set
##
##  Site filtering is reciprocal in order to maintain that both the gnomAD and
##  the HS_CAPTURE are only analyzed for variant burden on the same set of
##  possible genomic positions/sites.  So, low quality sites in HS_CAPTURE are
##  removed from both HS_CAPTURE and the gnomAD EXOME data set.  And likewise, sites
##  that are low quality in both gnomAD EXOME data set are removed from the gnomAD
##  set and the HS_CAPTURE.
##
##  This is very similar to 01.4_sites_quality_control.py, but we are not using
##  the gnomad GENOME data.  The GENOME data does not include SAS genomes.  Therefore,
##  the GENOMES data cannot be used for PCA purposes to identify SAS individuals.
##  So for the PCA ancestry analysis, we will use only the gnomAD EXOME data.
##
##  The variants that pass this quality control are then filtered so that the only variants that
##  are kept are ones for which a captured HS-affected individual (to be included in the burdens analysis)
##  has an alt allele.
##
##

##########################
##                      ##
##  Filter funciton     ##
##                      ##
##############################################################
def filter_out_nonHS_patients( variant_list , cohortSet ):
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

#######################################################

sites_dir = '../../output/sites_files'
variants_dir = '../../output/variant_files'
pca_dir = '../../output/pca_results'

gnomad_coverages_dir = '../../output/gnomad_coverages_files'
capture_coverages_dir = '../../data/capture_coverage_file'
info_dir = '../../info'
output_samples_info_dir = '../../output/samples_info_output'

print 'Reading in files...'
print

########################################################
##  Read in gnomad EXOME sites COVERAGE file.
##  Read in gnomad EXOME missense/synonymous files, and create VARIANT_SETs for sites.

gnoEXOME_coverages_file = '%s/%s' % ( gnomad_coverages_dir , 'gnomad_exomes_coverages_summary_table_CAPTURE_REGION.txt' )

gnoEXOME_SITES_file = '%s/%s' % ( sites_dir , 'gnomad_exomes_sites_CAPTURE_REGION.txt' )
gnoEXOME_SITES_vars = extractVariants( gnoEXOME_SITES_file )
gnoEXOME_SITES_variants = [ VARIANT( x , gnoEXOME_SITES_vars[0] ) for x in gnoEXOME_SITES_vars[1:] ]

######################################################
##  Read in HS CAPTURE sites COVERAGE file.
##  Read in HS CAPTURE missense/synonymous files, and create VARIANT_SETs for sites.

capture_coverage_file = '%s/%s' % ( capture_coverages_dir , 'HS_capture_all_samples_fosmid_coverages.txt' )

capture_sites_file = '%s/%s' % ( sites_dir , 'DJMH_016_029_joint_call_sites.txt' )
capture_sites_vars = extractVariants(capture_sites_file)
capture_sites_variants = [ VARIANT( x , capture_sites_vars[0] ) for x in capture_sites_vars[1:] ]


#######################################################
## Establish BED capture region dictionary for filtering purposes
##          key = chromosome# ; value = list of lists start/end of regions on the chromosome that are in the capture.  ex) key = 14 ; value = [ [ 73599404 , 73701440 ] , [ 105601962 , 105642819 ] ]

bed_dict = {}

capture_bed_file = '%s/%s' % ( info_dir , 'HS_capture_fosmid_regions_no_overlaps.bed' )
capture_bed_vars = extractVariants( capture_bed_file )

for x in capture_bed_vars:
	bed_chrom = x[0]
	bed_start = float(x[1])
	bed_end = float(x[2])

	if bed_chrom in bed_dict:
		bed_dict[bed_chrom].append([ bed_start , bed_end ])
	else:
		bed_dict[bed_chrom] = []
		bed_dict[bed_chrom].append([ bed_start , bed_end ])

##########################################################
#########  Determine positions to filter from EXOMES    ##
##########################################################

#
#  gnomAD coverage files are structured as tab separated vaules:
#      chrom   pos     mean    median  over_1  over_5  over_10 over_15 over_20 over_25 over_30 over_50 over_100
#
#  Will filter sites where the MEDIAN depth was < 20.
#

MEDIAN_DEPTH_COVERAGE_CUTOFF = 20.0

##########################
#  gnomad EXOME FILTERING
#

#print 'Determining gnomad EXOMES positions to filter...'
#print

gnomad_EXOME_positions_to_filter = []
gnomad_EXOME_PASS_positions = []

is_first_line = True
covg_index_dictionary = {}

file_handle = open(gnoEXOME_coverages_file, 'rb')

for line in file_handle:
	if is_first_line is True:
		is_first_line = False

		tmp_header_line = line.strip('\n').split('\t')
		index_counter = 0

		for header_field in tmp_header_line:
			covg_index_dictionary[header_field] = index_counter
			index_counter += 1
	else:
		tmp_line = line.strip('\n').split('\t')
		covg_chrom  = tmp_line[ covg_index_dictionary['chrom'] ]
		covg_pos = float( tmp_line[ covg_index_dictionary['pos'] ])
		covg_20fold = float( tmp_line[ covg_index_dictionary['over_20'] ])
		covg_median = float( tmp_line[ covg_index_dictionary['median'] ])

		if covg_chrom in bed_dict:
			for region in bed_dict[covg_chrom]:
				if covg_pos >= region[0] and covg_pos <= region[1]:
					if covg_median < MEDIAN_DEPTH_COVERAGE_CUTOFF:
						gnomad_EXOME_positions_to_filter.append( '%s_%s' % ( str(int(covg_chrom)) , str(int(covg_pos)) ) )
					else:
						gnomad_EXOME_PASS_positions.append( '%s_%s' % ( str(int(covg_chrom)) , str(int(covg_pos)) ) )
				else:
					None
		else:
			None

file_handle.close()

#############################################################################################
##     Identify positions from HS_CAPTURE that are low coverage and need to be filtered.    #
#############################################################################################

print 'Determining capture low_depth positions to filter...'
print

capture_low_depth_filter_sites = []
capture_low_depth_PASS_sites = []

is_first_line = True
file_handle = open(capture_coverage_file, 'rb')

sample_count = 0

for line in file_handle:
	if is_first_line is True:
		is_first_line = False

		tmp_header_line = line.strip('\n').split('\t')
		sample_count = ( len(tmp_header_line) - 2 )         # This is kind of an ugly non-universal way for handling this.  Capture coverages file was generated
		                                                    # as:  'CHROM\tPOSITION\tSAMPLE_1\tSAMPLE_2\tSAMPLE_3\t......\tLAST_SAMPLE'

	else:
		tmp_covgs = numpy.zeros(sample_count)

		tmp_line = line.strip('\n').split('\t')

		covg_chrom  = tmp_line[0]
		covg_pos = float(tmp_line[1])

		index_tracker = 0
		for sample_covg in tmp_line[2:]:
			tmp_covgs[index_tracker] = float(sample_covg)
			index_tracker += 1

		covg_median = float(numpy.median(tmp_covgs))

		if covg_chrom in bed_dict:
			for region in bed_dict[covg_chrom]:
				if covg_pos >= region[0] and covg_pos <= region[1]:
					if covg_median < MEDIAN_DEPTH_COVERAGE_CUTOFF:
						capture_low_depth_filter_sites.append('%s_%s' % ( str(int(covg_chrom)) , str(int(covg_pos)) ) )
					else:
						capture_low_depth_PASS_sites.append('%s_%s' % ( str(int(covg_chrom)) , str(int(covg_pos)) ) )
				else:
					None
		else:
			None

file_handle.close()

####################################################################################
#########  Determine uncalled_genotype positions to filter from HS_CAPTURE         #
####################################################################################

print 'Determining capture high uncalled GT positions to filter...'
print

UNKNOWN_CALL_MAXIMUM = 0.1    #Percent of alleles uncalled... i.e., percent of alleles called as '.', rather than ref(0) or alt(1)

capture_uncalled_filter_sites = []

for x in capture_sites_variants:
	chrm = x.info['CHROM']
	position = x.info['POS']
	chrm_pos = '%s_%s' % ( chrm , position )

	sample_number = len(x.samples_list)
	uncalled_alleles = x.all_genotypes_string.count('.')

	if float(uncalled_alleles)/float(sample_number) >= UNKNOWN_CALL_MAXIMUM:
		capture_uncalled_filter_sites.append(chrm_pos)
	else:
		None

##########################################
####                                    ##
####    COMMENCE FILTERING              ##
####								    ##
##########################################

print 'Filtering gnomad and capture files...'
print

#filtering lists:
#
#	capture_uncalled_filter_sites
#	capture_low_depth_filter_sites
#	gnomad_EXOME_positions_to_filter
#
#

#gnomad_EXOME_positions_to_filter
#gnomad_EXOME_PASS_positions

####
####  filter gnomad EXOMES

gnoEXOME_FILTERED_variants = []
gnoEXOME_PASSED_variants = []

for x in gnoEXOME_SITES_variants:
	chrom_pos = '%s_%s' % ( x.info['CHROM'] , x.info['POS'] )

	if chrom_pos in gnomad_EXOME_positions_to_filter or chrom_pos in capture_low_depth_filter_sites or chrom_pos in capture_uncalled_filter_sites:
		gnoEXOME_FILTERED_variants.append(x)
	else:
		gnoEXOME_PASSED_variants.append(x)

####
#### filter capture

capture_FILTERED_variants = []
capture_PASSED_variants = []

for x in capture_sites_variants:
	chrom_pos = '%s_%s' % ( x.info['CHROM'] , x.info['POS'] )

	if chrom_pos in capture_low_depth_filter_sites or chrom_pos in capture_uncalled_filter_sites or chrom_pos in gnomad_EXOME_positions_to_filter:
		capture_FILTERED_variants.append(x)
	else:
		capture_PASSED_variants.append(x)

#print 'Returning filtered variant data to output files...'

#
#  Return filtered variants to output files.
#

gnoEXOME_PASSED_variants_SET = VARIANT_SET( gnoEXOME_PASSED_variants )
gnoEXOME_PASSED_RETURN = gnoEXOME_PASSED_variants_SET.returnSet()
outPut_NODATE( gnoEXOME_PASSED_RETURN , '%s/%s' % ( pca_dir , 'PCA_gnomad_EXOME_PASSED_variants.txt')	)

###
### Filter capture variants so that the list only contains variants from HS affected
### individuals whose samples were actually captured and analyzed.
###

cohort_file = '%s/%s' % ( output_samples_info_dir , 'all_cohorts_HS_AFFECTED_captured_samples.txt' )
cohort_vars = extractVariants(cohort_file)
cohort_dict = {}
for x in cohort_vars[1:]:
	cohort_dict[x[0]] = 1

cohort_set = set(cohort_dict.keys())

HS_captured_samples_PASSED_variants , HS_captured_samples_FILTERED_variants = filter_out_nonHS_patients( capture_PASSED_variants , cohort_set )

outPut_NODATE( HS_captured_samples_PASSED_variants , '%s/%s' % ( pca_dir , 'PCA_capture_PASSED_variants_capturedHSpatients_only.txt') )

print 'Number of sites filtered from EXOME database:\t%i' % ( len( gnoEXOME_FILTERED_variants ) )
print 'Number of sites passed   from EXOME database:\t%i' % ( len( gnoEXOME_PASSED_variants ) )
print
print 'Number of sites filtered from CAPTURE database:\t%i' % ( len( capture_FILTERED_variants ) )
print 'Number of sites passed   from CAPTURE database:\t%i' % ( len( capture_PASSED_variants ) )
