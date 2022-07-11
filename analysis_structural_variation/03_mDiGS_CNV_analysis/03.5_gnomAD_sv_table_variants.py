#!/usr/bin/env python

##########
# Import #
##########
import numpy
import re
import sys

from burdens_code import *

###########################################################################################
###########################################################################################
###
###    03.5_gnomAD_sv_table_variants.py
###
###    This script is very similar to 03.2.  It prints out variants from the gnomAD SV database,
###    but it prints out variants for the paper's final HS-capture/gnomAD CNV comparison table.  There were two CNV identified in our
###	   HS cohort that were not polymorphisms.  One is not in the gnomAD variant databse as far as I can tell. The other looks like
###	   it is in the gnomAD SV database, but with a slightly different CNV positional window.  This script prints out both the
###    gnomAD polymorphisms and the non-polymorphic capture CNV.
###
###
###

gnomad_sv_vars_dir = '../../output/gnomad_structural_variants'
results_dir = '../../results'
sv_file = '%s/%s' % ( gnomad_sv_vars_dir , 'gnomad_structural_variants_in_capture_region.txt' )
sv_vars = extractVariants( sv_file )

##
## Parse header returned by bcftools from original gnomAD SV extraction in 03.1
##

temp_list = []
header = []
for x in sv_vars[0]:
	col = x.strip('# ').split(']')[-1]
	header.append(col)

temp_list.append(header)

##
## Extract polymorphisms...
##

header_index_dictionary = {}
index_tracker = 0
for x in header:
	header_index_dictionary[x] = index_tracker
	index_tracker += 1

for struc_var in sv_vars[1:]:

	frequencies_list = ['EUR_AF', 'AFR_AF', 'AMR_AF', 'EAS_AF']

	var_length = float( struc_var[ header_index_dictionary['SVLEN']])

	if all( [ float(struc_var[ header_index_dictionary[ x ] ]) >= 0.01 for x in frequencies_list ] ) and var_length >= 500:

		temp_list.append( struc_var )

	elif int(struc_var[header_index_dictionary['CHROM']]) == 4 and int(struc_var[header_index_dictionary['POS']] ) == 26255264:
		temp_list.append( struc_var )
	else:
		None

out_HEADER = ['CHROM' , 'START' , 'END', 'TYPE' , 'LENGTH', 'all_AF', 'afr_AF', 'amr_AF' , 'eas_AF', 'eur_AF' , 'AC', 'AN']
out_list = []
out_list.append( out_HEADER )

for x in temp_list[1:]:
	temp_CHROM = x[ header_index_dictionary['CHROM'] ]
	temp_START = int(x[ header_index_dictionary['POS'] ])
	temp_END = temp_START + int(x[ header_index_dictionary['SVLEN'] ])
	temp_TYPE = x[ header_index_dictionary['ALT'] ]
	temp_LENGTH = x[ header_index_dictionary['SVLEN'] ]
	temp_ALL_AF  = x[ header_index_dictionary['AF'] ]
	temp_afr_AF  = x[ header_index_dictionary['AFR_AF'] ]
	temp_amr_AF  = x[ header_index_dictionary['AMR_AF'] ]
	temp_eas_AF  = x[ header_index_dictionary['EAS_AF'] ]
	temp_eur_AF  = x[ header_index_dictionary['EUR_AF'] ]
	temp_AC  = x[ header_index_dictionary['AC'] ]
	temp_AN  = x[ header_index_dictionary['AN'] ]

	temp_entry = []
	temp_entry.append(temp_CHROM)
	temp_entry.append(temp_START)
	temp_entry.append(temp_END)
	temp_entry.append(temp_TYPE)
	temp_entry.append(temp_LENGTH)
	temp_entry.append(temp_ALL_AF)
	temp_entry.append(temp_afr_AF)
	temp_entry.append(temp_amr_AF)
	temp_entry.append(temp_eas_AF)
	temp_entry.append(temp_eur_AF)
	temp_entry.append(temp_AC)
	temp_entry.append(temp_AN)

	out_list.append( temp_entry )

for x in out_list:
	print x

outPut_NODATE( out_list , '%s/%s' % ( results_dir, 'gnomad_variants_for_paper_table.txt') )
