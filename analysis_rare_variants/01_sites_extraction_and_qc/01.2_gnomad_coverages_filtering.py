#!/usr/bin/env python

##########
# Import #
##########
import gzip
import numpy
import sys

from burdens_code import *

##############################################################
##############################################################
##
##  01.2_GNOMAD_COVERAGES_FILTERING.PY
##
##  We will quality filter genomic positions based on coverage depth at the position.  To do this we will use both our own capture
##         data and gnomad coverage information.  Gnomad only reports population coverage metrics, so we will end up having to
##         use gnomad median coverage at positions to filter out low coverage sites from gnomad.  Gnomad coverage tables are very large
##         This script filters the gnomAD coverage table to keep only information for positions that fall within our targeted capture
##         regions.  This will give us a smaller gnomAD coverages file that will be easier to work with when we perform the actual
##         quality filtering of variants.
##

def filterCoverageFiles( bed_file , coverages_file ):
	####
	#### Establish BED capture region dictionary for filtering purposes

	bed_dict = {}
	covg_index_dictionary = {}

	bed_vars = extractVariants( bed_file )

	for x in bed_vars:

		bed_chrom = x[0]
		bed_start = float(x[1])
		bed_end = float(x[2])

		if bed_chrom in bed_dict:
			bed_dict[bed_chrom].append([ bed_start , bed_end ])
		else:
			bed_dict[bed_chrom] = []
			bed_dict[bed_chrom].append([ bed_start , bed_end ])

	##################################################################
	#########  Determine positions to filter from gnomad GENOMES    ##
	##################################################################

	positions_in_capture_region = []

	is_first_line = True
	file_handle = gzip.open(coverages_file, 'rb')
	for line in file_handle:
		if is_first_line is True:

			is_first_line = False

			tmp_header_line = line.strip('\n').split('\t')

			positions_in_capture_region.append(tmp_header_line)

			index_counter = 0
			for header_field in tmp_header_line:
				#print header_field
				covg_index_dictionary[header_field] = index_counter
				index_counter += 1

		else:
			tmp_line = line.strip('\n').split('\t')
			covg_chrom  = tmp_line[ covg_index_dictionary['chrom'] ]

			if covg_chrom in bed_dict:

				covg_pos = float( tmp_line[ covg_index_dictionary['pos'] ])

				for region in bed_dict[covg_chrom]:
					if covg_pos >= region[0] and covg_pos <= region[1]:
						positions_in_capture_region.append( tmp_line )
					else:
						None
			else:
				None

	file_handle.close()

	return( positions_in_capture_region )

##########################
##########################

coverage_file_dir = '../../data/gnomad_coverages_files'
info_dir = '../../info'
out_dir = '../../output/gnomad_coverages_files'

gnomad_genomes_coverages_file = '%s/%s' % ( coverage_file_dir , 'gnomad.genomes.coverage.summary.tsv.gz' )
gnomad_exomes_coverages_file = '%s/%s' % ( coverage_file_dir , 'gnomad.exomes.coverage.summary.tsv.gz' )
capture_bed_file = '%s/%s' % ( info_dir , 'HS_capture_fosmid_regions_no_overlaps.bed' )

print 'starting genome coverages filtering...'
gnomad_genomes_capture_region = filterCoverageFiles( capture_bed_file , gnomad_genomes_coverages_file )
gnomad_exomes_capture_region = filterCoverageFiles( capture_bed_file , gnomad_exomes_coverages_file )

outPut_NODATE( gnomad_genomes_capture_region , '%s/%s' % ( out_dir , 'gnomad_genomes_coverages_summary_table_CAPTURE_REGION.txt') )
outPut_NODATE( gnomad_exomes_capture_region , '%s/%s' % ( out_dir , 'gnomad_exomes_coverages_summary_table_CAPTURE_REGION.txt' ) )
