#!/usr/bin/env python

##########
# Import #
##########
import numpy
import sys

from burdens_code import *

##############################################################
##############################################################
##
##  01.3_GNOMAD_GENOMIC_SAS_APPEND.PY
##
##  gnomAD genome data does not include SAS population subgroup because the sample size for SAS individuals was too small.
##         This script adds an SAS subpopulation category to the genomic gnomAD data.  It adds alt allele counts and total allele number
##         of 0.  This is just a place holder to make it easier to merge the exome and genome sites data from gnomAD.  The SAS
##         group information was intentionally ordered as last subgroup in the exome vcf sites extraction because SAS was going to be appended
##         as the last category for the genomic data, and it was just easier to merge the data with all columns in the same order for both data sets.
##         This is all just to say that the subpopulation group extraction order for the exome and genome data is important.
##

sites_dir = '../../output/sites_files'

genomic_sites_file = '%s/%s' % ( sites_dir , 'gnomad_genomes_sites_CAPTURE_REGION.txt')
gen_sites_vars = extractVariants( genomic_sites_file )

GEN_variants = [ VARIANT( x , gen_sites_vars[0] ) for x in gen_sites_vars[1:] ]

for x in GEN_variants:
	x.appendInfo( 'AC_sas' , '0' )
	x.appendInfo( 'AN_sas' , '0' )

GEN_set = VARIANT_SET( GEN_variants )

GEN_return = GEN_set.returnSet()

outPut_NODATE( GEN_return , '%s/%s' % ( sites_dir , 'gnomad_genomes_sites_CAPTURE_REGION_wSAS.txt' ) )
