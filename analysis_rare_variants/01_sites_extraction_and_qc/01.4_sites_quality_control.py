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
##  01.4_SITES_QUALITY_CONTROL.PY
##
##  This script removes positions/sites from the gnomAD EXOME, gnomad GENOME,
##  and HS_CAPTURE data sets based on:
##        - coverage depth at the site
##        - number of alleles called as unknown for HS CAPTURE data set
##
##  Site filtering is reciprocal in order to maintain that both the gnomAD and
##  the HS_CAPTURE are only analyzed for variant burden on the same set of
##  possible genomic positions/sites.  So, low quality sites in HS_CAPTURE are
##  removed from both HS_CAPTURE and the gnomAD data sets.  And likewise, sites
##  that are low quality in both gnomAD data sets (exome and genome) are removed
##  from the gnomAD sets and the HS_CAPTURE.
##

gnomad_coverages_dir = '../../output/gnomad_coverages_files'
capture_coverages_dir = '../../data/capture_coverage_file'
sites_dir = '../../output/sites_files'
info_dir = '../../info'

########################################################################################
##  Read in gnomad GENOME sites COVERAGE file.
##  Read in gnomad GENOME missense/synonymous files, and create VARIANT_SETs for sites.

print 'Reading in files...'
print

gnoGENOME_coverages_file = '%s/%s' % ( gnomad_coverages_dir , 'gnomad_genomes_coverages_summary_table_CAPTURE_REGION.txt' )

gnoGENOME_SITES_file = '%s/%s' % ( sites_dir , 'gnomad_genomes_sites_CAPTURE_REGION_wSAS.txt' )
gnoGENOME_SITES_vars = extractVariants( gnoGENOME_SITES_file )
gnoGENOME_SITES_variants = [ VARIANT( x , gnoGENOME_SITES_vars[0] ) for x in gnoGENOME_SITES_vars[1:] ]

########################################################################################
##  Read in gnomad EXOME sites COVERAGE file.
##  Read in gnomad EXOME missense/synonymous files, and create VARIANT_SETs for sites.

gnoEXOME_coverages_file = '%s/%s' % ( gnomad_coverages_dir , 'gnomad_exomes_coverages_summary_table_CAPTURE_REGION.txt' )

gnoEXOME_SITES_file = '%s/%s' % ( sites_dir , 'gnomad_exomes_sites_CAPTURE_REGION.txt' )
gnoEXOME_SITES_vars = extractVariants( gnoEXOME_SITES_file )
gnoEXOME_SITES_variants = [ VARIANT( x , gnoEXOME_SITES_vars[0] ) for x in gnoEXOME_SITES_vars[1:] ]

########################################################################################
##  Read in HS CAPTURE sites COVERAGE file.
##  Read in HS CAPTURE missense/synonymous files, and create VARIANT_SETs for sites.

capture_coverage_file = '%s/%s' % ( capture_coverages_dir , 'HS_capture_all_samples_fosmid_coverages.txt' )

capture_sites_file = '%s/%s' % ( sites_dir , 'DJMH_016_029_joint_call_sites.txt' )
capture_sites_vars = extractVariants(capture_sites_file)
capture_sites_variants = [ VARIANT( x , capture_sites_vars[0] ) for x in capture_sites_vars[1:] ]

##################################################################
## Establish BED capture region dictionary for filtering purposes
##          key = chromosome+# ; value = list of lists start/end of regions on the chromosome that are in the capture.  ex) key = 14 ; value = [ [ 73599404 , 73701440 ] , [ 105601962 , 105642819 ] ]

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

###########################
#  gnomad GENOME FILTERING
#

#print 'Determining gnomad GENOME positions to filter...'
#print

gnomad_GENOME_positions_to_filter = []
gnomad_GENOME_PASS_positions = []

is_first_line = True

file_handle = open(gnoGENOME_coverages_file, 'rb')

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

		if covg_chrom in bed_dict:
			covg_pos = float( tmp_line[ covg_index_dictionary['pos'] ])
			covg_20fold = float( tmp_line[ covg_index_dictionary['over_20'] ])
			covg_median = float( tmp_line[ covg_index_dictionary['median'] ])

			for region in bed_dict[covg_chrom]:
				if covg_pos >= region[0] and covg_pos <= region[1]:
					if covg_median < MEDIAN_DEPTH_COVERAGE_CUTOFF:
						gnomad_GENOME_positions_to_filter.append('%s_%s' % ( str(int(covg_chrom)) , str(int(covg_pos)) ) )
					else:
						gnomad_GENOME_PASS_positions.append('%s_%s' % ( str(int(covg_chrom)) , str(int(covg_pos)) ) )
				else:
					None
		else:
			None

file_handle.close()

#########################################################################################################################################
###		                                                                                                                                #
###	 MERGE  GNOMAD  EXOME  AND  GENOME  FILTER  POSITIONS																				#
###        Compile filtered sites from gnomad to filter from HS_CAPTURE.									                            #
###       (If a site is low_depth filtered in EXOMEs but passed in GENOMES (or vice versa), then do NOT filter site from HS_CAPTURE)    #
### 																												                    #
#########################################################################################################################################

print 'Determining gnomad positions that are low depth in both genomes and exomes...'
print

#gnomad_GENOME_positions_to_filter
#gnomad_GENOME_PASS_positions

#gnomad_EXOME_positions_to_filter
#gnomad_EXOME_PASS_positions

gnomad_positions_to_filter_from_HS_CAPTURE = []

for position in gnomad_GENOME_positions_to_filter:
	if position in gnomad_EXOME_PASS_positions:
		None
	else:
		gnomad_positions_to_filter_from_HS_CAPTURE.append(position)

for position in gnomad_EXOME_positions_to_filter:
	if position in gnomad_GENOME_PASS_positions:
		None
	else:
		gnomad_positions_to_filter_from_HS_CAPTURE.append(position)

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
#	gnomad_positions_to_filter_from_HS_CAPTURE
#	gnomad_GENOME_positions_to_filter
#	gnomad_EXOME_positions_to_filter
#
#
#output lists:
#
#	gnoGENOME_MISSENSE_filtered_variants			gnoGENOME_MISSENSE_PASS_variants
#	gnoGENOME_SYNON_filtered_variants            	gnoGENOME_SYNON_PASS_variants
#
#
#	gnoEXOME_MISSENSE_filtered_variants				gnoEXOME_MISSENSE_PASS_variants
#	gnoEXOME_SYNON_filtered_variants				gnoEXOME_SYNON_PASS_variants
#
#
#	capture_MISSENSE_filtered_variants				capture_MISSENSE_PASS_variants
#	capture_SYNON_filtered_variants					capture_SYNON_PASS_variants
#

####
####  filter gnomad GENOMES


####################################################
##   HOW  TO  FILTER??
##
##   In looking back through code while cleaning and editing, I started wondering if I was using the right filtering
##      strategy.  I had to re-convince myself that I think I am.  I thought I would explain things more here.
##
##   Right now I am filtering as follows:
##         GENOME_SITES:
##         - filter gnomad GENOME positions from the GENOME sites set if they are below MIN_MEDIAN_CUTOFF.
##               - do not automatically filter this position from EXOME or HS_CAPTURE SET
##         - filter gnomad GENOME positions if the position is below MIN_MEDIAN_CUTOFF in HS_CAPTURE
##         - filter gnomad GENOME positions if the position is above MAX_UNKNOWN_ALLELE_CUTOFF in HS_CAPTURE
##
##
##         EXOME_SITES:
##         - filter gnomad EXOME positions from the EXOME sites set if they are below MIN_MEDIAN_CUTOFF.
##               - do not automatically filter this position from GENOME or HS_CAPTURE SET
##         - filter gnomad EXOME positions if the position is below MIN_MEDIAN_CUTOFF in HS_CAPTURE
##         - filter gnomad EXOME positions if the position is above MAX_UNKNOWN_ALLELE_CUTOFF in HS_CAPTURE
##
##
##         HS_CAPTURE_SITES:
##         - filter gnomad HS_CAPTURE positions from the HS_CAPTURE sites set if they are below MIN_MEDIAN_CUTOFF.
##         - filter gnomad HS_CAPTURE positions if the position is above MAX_UNKNOWN_ALLELE_CUTOFF
##         - filter gnomad HS_CAPTURE positions if the position NOT above MIN_MEDIAN_CUTOFF in BOTH the gnomad GENOME and EXOME data sets
##
##
##         Thoughts:
##            The idea is to have a reciprocal depth filter.  If the site is bad depth in the capture, then we should
##				 also filter it from the gnomad data set.  The same is true the other way around EXCEPT that a single
##               'BELOW_MINIMUM_COVERAGE' flag does not result in filtering the site from all data sets.  Rather, a site
##               from gnomAD has to be either unpresent or below coverage cutoff for both data sets in order to remove the position
##               from the capture dataset.
##
##               I think that this makes sense.  One thing I have struggled with is deciding if a single 'ABOVE_MINUM_COVERAGE'
##               hit in the gnomAD data set should allow for both data sets to be kept and combined.  Or if a single 'BELOW_MINUMUM_COVERAGE'
##               hit should trigger removal from both EXOME and GENOME data sets.
##
##               My reasoning for the current filtering logic is as follows:

##               Should we require two coverage fails before filtering a site from gnomAD? The idea behind the depth coverage filter is that low
##               coverage reflects a lower confidence in genotype calls.  If GENOME is below coverage minimum, but EXOME is above,
##               summing the site data from GENOME and EXOME does not increase confidence in the original GENOME genotype call... thus
##               a coverage PASS for one data set should not rescue the site for the other data set.
##
##               Should we filter site from both gnomAD sets if coverage threshold isn't met in only one of the sets?
##               The GENOME and EXOME data sets have very different properties with regards to what regions of the genome will be
##               sequenced with good coverage.  A site that is intronic near the end of an exon will likely have a EXOME coverage that
##               is low.  But the GENOME data set might have sufficient coverage for the same intronic site.  Allowing to keep sites
##               that fail depth coverage cutoff in one set but pass the threshold in the other allows us to maximize the
##               amount of useful data we extract from the data sets.
##
##
##

gnoGENOME_FILTERED_variants = []
gnoGENOME_PASSED_variants = []

for x in gnoGENOME_SITES_variants:
	chrom_pos = '%s_%s' % ( x.info['CHROM'] , x.info['POS'] )

	if chrom_pos in gnomad_GENOME_positions_to_filter or chrom_pos in capture_low_depth_filter_sites or chrom_pos in capture_uncalled_filter_sites:
		gnoGENOME_FILTERED_variants.append(x)
	else:
		gnoGENOME_PASSED_variants.append(x)

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

	if chrom_pos in capture_low_depth_filter_sites or chrom_pos in capture_uncalled_filter_sites or chrom_pos in gnomad_positions_to_filter_from_HS_CAPTURE:
		capture_FILTERED_variants.append(x)
	else:
		capture_PASSED_variants.append(x)

#print 'Returning filtered variant data to output files...'

#
#  Return filtered variants to output files.
#

gnoGENOME_FILTERED_variants_SET = VARIANT_SET( gnoGENOME_FILTERED_variants )
gnoGENOME_FILTERED_RETURN = gnoGENOME_FILTERED_variants_SET.returnSet()
outPut_NODATE( gnoGENOME_FILTERED_RETURN , '%s/%s' % ( sites_dir , 'gnomad_GENOME_QC_filtered_sites.txt' ) )

gnoGENOME_PASSED_variants_SET = VARIANT_SET( gnoGENOME_PASSED_variants )
gnoGENOME_PASSED_RETURN = gnoGENOME_PASSED_variants_SET.returnSet()
outPut_NODATE( gnoGENOME_PASSED_RETURN , '%s/%s' % ( sites_dir , 'gnomad_GENOME_QC_passed_sites.txt' ) )

##
##

gnoEXOME_FILTERED_variants_SET = VARIANT_SET( gnoEXOME_FILTERED_variants )
gnoEXOME_FILTERED_RETURN = gnoEXOME_FILTERED_variants_SET.returnSet()
outPut_NODATE( gnoEXOME_FILTERED_RETURN , '%s/%s' % ( sites_dir , 'gnomad_EXOME_QC_filtered_sites.txt' ) )

gnoEXOME_PASSED_variants_SET = VARIANT_SET( gnoEXOME_PASSED_variants )
gnoEXOME_PASSED_RETURN = gnoEXOME_PASSED_variants_SET.returnSet()
outPut_NODATE( gnoEXOME_PASSED_RETURN , '%s/%s' % ( sites_dir , 'gnomad_EXOME_QC_passed_sites.txt' ) )

##
##

capture_FILTERED_variants_SET = VARIANT_SET( capture_FILTERED_variants )
capture_FILTERED_RETURN = capture_FILTERED_variants_SET.returnSet()
outPut_NODATE( capture_FILTERED_RETURN , '%s/%s' % ( sites_dir , 'CAPTURE_QC_filtered_sites.txt' ) )

capture_PASSED_variants_SET = VARIANT_SET( capture_PASSED_variants )
capture_PASSED_RETURN = capture_PASSED_variants_SET.returnSet()
outPut_NODATE( capture_PASSED_RETURN , '%s/%s' % ( sites_dir , 'CAPTURE_QC_passed_sites.txt' ) )
