#!/usr/bin/env python

##########
# Import #
##########
import numpy
import sys

from burdens_code import *

##################
##################

########################################################################
###
###   PRE.1_ORGANIZE_BIOMART.PY
###   This script takes fasta output from biomart and organizes it into
###   a tab separated table where each line is exon information for the
###   transcripts to be analyzed for the variant burdens analysis.
###

info_dir = '../../info'
transcript_dir = '../../output/transcript_file'

biomart_file = '%s/%s' % ( info_dir , 'hs_capture_primary_transcript_info_biomart_unorganized.txt' )

out_list = []
InFile = open( biomart_file , 'r' )
NO_SEQUENCES_YET = True
for LINE in InFile:
	if LINE[0] == '#':
		HEADER = LINE.strip('#').strip('\r').strip('\n').split('|')
	elif LINE[0] == '>':
		if NO_SEQUENCES_YET:
			transcript_info = LINE.strip('>').strip('\r').strip('\n').split('|')
			SEQUENCE = ''
			NO_SEQUENCES_YET = False
		else:
			transcript_info.append(SEQUENCE)
			out_list.append(transcript_info)

			transcript_info = LINE.strip('>').strip('\r').strip('\n').split('|')
			SEQUENCE = ''
	else:
		SEQUENCE += LINE.strip('>').strip('\r').strip('\n')

transcript_info.append(SEQUENCE)
out_list.append(transcript_info)

InFile.close()

out_list.sort(key=lambda z: ( int(z[0]) , int(z[2]) ) )

final_out = []
final_out.append(HEADER)

for x in out_list:
	final_out.append(x)

outPut_NODATE( final_out , '%s/%s' % ( transcript_dir , 'HS_capture_primary_transcript_exon_sequences.txt' ) )
