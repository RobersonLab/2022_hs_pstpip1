#!/usr/bin/env python

###########
# Imports #
###########
import numpy
import random
import re
import sys

###########
# Classes #
###########
class VARIANT:
	def __init__ (self,vcf_info_list,vcf_info_header):
		#bcftools header seems to be returned as '# [1]FIELD\t[2]FIELD\t[3]FIELD\n

		self.header_dictionary = {}                    # key = FIELD (i.e. column name) --> value = INDEX number ; i.e. the column number where the information is located.  Thus given a list you can use the INDEX from the header dictionary to hash to the entry with the correct information
		self.header_field_list = []                    # list of column names for fields that are not SAMPLE specific; i.e. CHROM, POS, REF, ALT, FILTER (filter status)
		self.sample_fields_list = []                   # list of raw column names for sample information columns; i.e. column header [5]SAMPLE1:GT --> sample_field_list entry = SAMPLE1:GT
		self.samples_dictionary = {}                   # dictionary of dictionaries: key = SAMPLE_NAME { key = INFORMATION_CATEGORY (i.e. GT for genotype, DP for depth) --> value = INDEX number ; i.e. the column number where the information is located.  Thus given a list you can use the INDEX from the samples dictionary to hash to the entry with the correct information}
		self.samples_list = []                         # list of sample names/IDs in the data set
		self.depths_list = []                          # list of the coverage depths of all samples for the variant;  this will help for variant filtering; will allow for quickly determining if mean/median depth at this variant meets filtering threshold
		self.gtquals_list = []                         # list of the gt_quality of all samples for the variant;  this will help for variant filtering; will allow for quickly determining if mean/median gt_quality at this variant meets filtering threshold
		self.depth_dictionary = {}                     # key = sample_name/ID --> value = coverage depth for that sample for this variant
		self.gtquals_dictionary = {}                   # key = sample_name/ID --> value = gt_quality for that sample for this variant
		self.tag_order_list = []                       # when reading in variant header information, keeps track of the order of sample information; i.e., the header contains [GT,DP] genotypes and then depths

		self.info = {} 								   #dictionary, keys = header columns of gemini search --> value = respective information for category requested ;  EXAMPLE:  VARIANT.self.info["CHROM"] = 1, VARIANT.self.info["POS"] = 105260, VARIANT.self.info["REF"] = A, VARIANT.self.info["ALT"] = T
		                                               #            ** for this dictionary, using sample_name/ID as a key will return genotype, i.e. VARIANT.self.info[SAMPLE_NAME]  = "0/1"
													   #            ** also can access other sample information with the more complete request, i.e. VARIANT.self.info[SAMPLE_NAME:GT] = "0/1" , VARIANT.self.info[SAMPLE_NAME:DP] = 80
													   #
													   #            ** this will be kind of the heart of this entire class; other aspects will be used for variant filtering;
													   #               but this will be used for the burden analysis... for instance: alt_frequency = VARIANT.self.info["aaf_gnomad_nfe"] = 0.00134

		self.all_genotypes_string = ''                 # compiled string of all sample genotypes for variant; i.e. "0/10/00/00/1./.0/0./00/1" ; this was just used a simple way of figuring out alt, ref and uncalled counts for the variant.  genotypes are just smushed into one string and then count number of "0"s, "1"s or "."s

	    #self.total_uncalled_alleles    (see below for instantiation)     # number of alleles for this variant that were reported as "."
	    #self.total_ref_alleles         (see below for instantiation)     # number of alleles for this variant that were reported as "0"
	    #self.total_alt_alleles         (see below for instantiation)     # number of alleles for this variant that were reported as "1"
	    #self.alt_affected_individuals  (see below for instantiation)     # number of samples that have at least one alt allele (i.e. at least one "1" in genotype)

		#######################
		# FUNCTIONS AVAILABLE #
		#######################

		#	def appendInfo(self,header_column_name,value)     # append information to the variant; ex) variant is created with sequencing information from our targeted capture, and then later gnomAD alt allele frequencies are added
		#	def patientAltCount(self,patient)			      # provide patient/sample_id , returns the number of alt alleles observed in that individual
		#	def returnVariant(self)                           # return variant information ; should look identical to information list used for creating VARIANT unless information has been appended
		#	def returnVariantHeaderInfo(self)				  # return variant summary information (i.e., chromosome, position, reference, alt, etc...; but doesn't print sample specific information)
		#	def printVariant(self)                            # prints variant information
		#	def printVariantHeaderInfo(self)                  # prints variant summary information
		#	def showAltIndividuals(self)                      # returns list of sample_names/ID that have alternate allele for this variant

		#read in header information for variant and sample cohort

		for each in vcf_info_header:
			cleaned = each.strip('#').strip(' ').strip('[').split(']')
			field_index = int(cleaned[0])-1
			field = cleaned[1]

			if ':' in field:

			## handle sample information fields
			## genotype information is extracted as "[COLUMN#]SAMPLE_NAME:INFORMATION_ABBREVIATION"  (ex: [5]sample_1:GT , [7]sample_2:DP)

				self.sample_fields_list.append(field)
				sample_field = field.split(':')[0]
				tmp_format_tag = field.split(':')[1]
				if tmp_format_tag in self.tag_order_list:
					None
				else:
					self.tag_order_list.append(tmp_format_tag)

				if sample_field in self.samples_list:
					self.samples_dictionary[sample_field][tmp_format_tag] = field_index
				else:
					self.samples_list.append(sample_field)
					self.samples_dictionary[sample_field] = {}
					self.samples_dictionary[sample_field][tmp_format_tag] = field_index

				self.header_dictionary[sample_field] = field_index                               ## I don't think that this is needed
			else:
				self.header_dictionary[field] = field_index
				self.header_field_list.append(field)


		# fill in formation for self.info = {} dictionary

		for each_field in self.header_field_list:
			self.info[each_field] = vcf_info_list[self.header_dictionary[each_field]]


		# fill in all genotypes string
		# fill in self.info sample genotype information
		# fill in depth and gt_qual lists and dictionaries

		for each_sample in self.samples_list:
			self.all_genotypes_string += vcf_info_list[self.samples_dictionary[each_sample]['GT']]   # I think this can probably be simplified to ' += self.info

			self.info[each_sample] = vcf_info_list[self.samples_dictionary[each_sample]['GT']]

			self.depths_list.append(vcf_info_list[self.samples_dictionary[each_sample]['DP']])
			self.gtquals_list.append(vcf_info_list[self.samples_dictionary[each_sample]['GQ']])
			self.depth_dictionary[each_sample] = vcf_info_list[self.samples_dictionary[each_sample]['DP']]
			self.gtquals_dictionary[each_sample] = vcf_info_list[self.samples_dictionary[each_sample]['GQ']]


		# fill in more self.info information

		for tag in self.tag_order_list:
			for sample in self.samples_list:
				self.info['%s:%s' % ( sample , tag ) ] = vcf_info_list[self.samples_dictionary[sample][tag]]



		# calculate summary information for number of uncalled, ref, alt alleles; determine number of samples that have an alt allele

		self.total_uncalled_alleles = self.all_genotypes_string.count('.')
		self.total_ref_alleles = self.all_genotypes_string.count('0')
		self.total_alt_alleles = self.all_genotypes_string.count('1')

		self.alt_affected_individuals = 0
		for each_sample in self.samples_list:
			if '1' in self.info[each_sample]:
				self.alt_affected_individuals += 1
			else:
				None

	def appendInfo(self,header_column_name,value):
		if header_column_name in self.header_field_list:
			self.info[header_column_name] = '%s/%s' % ( self.info[header_column_name] , value )
			#print 'sorry, \'%s\' information already exists in this gemini variant'
			#print 'adding additional info to existing header field'
		else:
			self.header_field_list.append(header_column_name)
			self.info[header_column_name] = value
			self.header_dictionary[header_column_name] = len(self.header_field_list) - 1

			for each in self.sample_fields_list:
				tmp_field = each.split(':')[0]
				tmp_tag = each.split(':')[1]
				self.samples_dictionary[tmp_field][tmp_tag] += 1

	def patientAltCount(self,patient):
		alt_alleles = self.info[patient].count('1')
		return(alt_alleles)

	def returnVariant(self):
		return_variant = []
		for x in self.header_field_list:
			return_variant.append(self.info[x])

		for tag in self.tag_order_list:
			for sample in self.samples_list:
				return_variant.append(self.info['%s:%s' % ( sample , tag )])

		return (return_variant)

	def returnVariantHeaderInfo(self):
		return_variant_header_info = []
		for x in self.header_field_list:
			return_variant_header_info.append(self.info[x])

		return (return_variant_header_info)

	def printVariant(self):
		variant_to_print = self.returnVariant()
		print variant_to_print

	def printVariantHeaderInfo(self):
		variant_to_print = self.returnVariantHeaderInfo()
		print variant_to_print


	def showAltIndividuals(self):
		alt_affected_samples = []
		for sample in self.samples_list:
			if '1' in self.info[sample]:
				alt_affected_samples.append(sample)
			else:
				None

		return (alt_affected_samples)

	def add_AAFs(self):

		##
		##  This finds AC and AN values for same populations and then calculates and appends
		##  the aaf_gnomad_subpopulation allele frequency.  This is used for generating AAFs
		##  for PCA analysis.  It is not used for the broader variant burdens analysis, because
		##  for the burdens analysis I compiled the gnomad GENOMEs and EXOMEs data.  To do that I
		##  had to go in manually to add the genomes AC and exomes AC, so I just also calculate the
		##  aaf in that code while I have the total AC and total AN on hand.
		##

		aaf_groups = []
		for x in self.header_field_list:
			if 'AC' in x:
				aaf_temp = x.split('_')
				if aaf_temp in aaf_groups:
					None
				else:
					aaf_groups.append(aaf_temp)
			else:
				None
		for x in aaf_groups:

			## Subgroups taken in as 'AC_nfe' split into ['AC','nfe'].
			## There was a time when I was playing around with the data and generated combined groups
			## like 'AC_nfe_afr', so code below determines if it is an individual or combined population group.

			if len(x) == 3:
				subgroup = '%s_%s' % ( x[1] ,  x[2] )
			elif len(x) == 2:
				subgroup = x[1]
			else:
				None

			alt_count = float(self.info['AC_%s' % (subgroup)])
			allele_count = float(self.info['AN_%s' % (subgroup)])
			aaf = '%.15f' % (alt_count/allele_count)

			self.appendInfo('aaf_gnomad_%s' % (x[-1]) , aaf)

class VARIANT_SET:
	##
	## Creates a set of VARIANTs.  Allows for easier filtering/processing/looping over VARIANTS
	##
	## VARIANTS must have identical information in the identical order;  the hashing information is simply
	##          taken from the first variant in the set and used for all samples;  identical field information/order
	##          is tested for upon VARIANT_SET instantiation.
	##

	def __init__ ( self , VARIANT_list ):

		if len(VARIANT_list) == 0:
			print 'there are no variants to SET'

		else:
			None
		###
		###  Verify that all variants have same header fields;  if not, throw error
		self.header_field_list = VARIANT_list[0].header_field_list
		self.sample_fields_list = VARIANT_list[0].sample_fields_list
		self.samples_list = VARIANT_list[0].samples_list

		self.header_dictionary = VARIANT_list[0].header_dictionary
		self.samples_dictionary = VARIANT_list[0].samples_dictionary
		self.depths_list = VARIANT_list[0].depths_list
		self.gtquals_list = VARIANT_list[0].gtquals_list
		self.depth_dictionary = VARIANT_list[0].depth_dictionary
		self.gtquals_dictionary = VARIANT_list[0].gtquals_dictionary
		# self.variants         # list of VARIANTS

		all_variant_fields_are_identical = True
		for x in VARIANT_list:
			if x.header_field_list == self.header_field_list and x.sample_fields_list == self.sample_fields_list and x.samples_list == self.samples_list:
				None
			else:
				all_variant_fields_are_identical = False
				print 'First variant in VARIANT list:'
				VARIANT_list[0].printVariantHeaderInfo()
				print 'Mismatched variant in list:'
				x.printVariantHeaderInfo()

		if all_variant_fields_are_identical is False:

			print 'variants in the set you are trying to create do not contain identical header info'
			sys.exit(1)
		else:
			None

		self.variants = [x for x in VARIANT_list]

	#######################
	# FUNCTIONS AVAILABLE #
	#######################
	#
	#	def addVariant(self , NEW_VARIANT)       # adds additional variant
	#	def returnSet(self)						 # return the information in the VARIANT_SET as a list of list of all the VARIANT information
	#

	def addVariant(self , NEW_VARIANT):
		if NEW_VARIANT.header_field_list == self.header_field_list and NEW_VARIANT.sample_fields_list == self.sample_fields_list and NEW_VARIANT.samples_list == self.samples_list :
			self.variants.append(NEW_VARIANT)
		else:
			print 'could not append VARIANT to VARIANT_set: header fields are not the same.'

	def generate_AAF_INFO(self):
		for x in self.variants:
			x.add_AAFs()

	def returnSet(self):
		return_list = []
		HEADER = []
		for x in self.header_field_list:
			tmp_field = '[%i]%s' % ((self.header_dictionary[x]+1),x)
			HEADER.append(tmp_field)

		for x in self.sample_fields_list:
			sample = x.split(':')[0]
			tag = x.split(':')[1]
			index = self.samples_dictionary[sample][tag]
			tmp_field= '[%i]%s:%s' % ( (index+1) , sample , tag )

			HEADER.append(tmp_field)

		return_list.append(HEADER)

		for each_variant in self.variants:
			tmp_var = each_variant.returnVariant()
			return_list.append(tmp_var)

		return (return_list)

	def prCompTable_w_simulations( self , NUMBER_OF_SIMULATED_training_INDIVIDUALS , NUMBER_OF_SIMULATED_validation_INDIVIDUALS ):

		##
		##  Create a variant PCA table for analyzing cohort samples to categorize individuals in to PCA ancestry populations.
		##  This converts each individuals genotype into 0,1, or 2 to represent the number of alt alleles the individual has
		##  for that variant.  Additionally, each variant is simulated, using a random number generator, 200 times for 5 gnomad
		##  subpopulations (NFE, AFR, AMR, SAS, and EAS).
		##
		##

		compiled_samples_list = []  # [x for x in self.variants[0].samples_list]
		var_counter  = 0
		allele_dictionary = {}

		for variant in self.variants:                      #  LOOP  OVER  VARIANTS
			var_counter += 1
			if var_counter % 100 == 0:
				print '%i variants processed' % (var_counter)
			else:
				None

			# track variants and individuals' genotypes at that variant with a dictionary using variant information as keys.

			var_id = '%s_%s_%s_%s' % (variant.info['CHROM'],variant.info['POS'],variant.info['REF'],variant.info['ALT'])
			allele_dictionary[var_id] = {}


			# SIMULATE variant genotypes for various gnomad_subpopulation groups for PCA TRAINING SET.

			allele_dictionary, compiled_samples_list = pca_subPopulation_simulation( 'training', variant, var_id, 'NFE', 'aaf_gnomad_nfe' , NUMBER_OF_SIMULATED_training_INDIVIDUALS, allele_dictionary, compiled_samples_list)
			allele_dictionary, compiled_samples_list = pca_subPopulation_simulation( 'training', variant, var_id, 'AFR', 'aaf_gnomad_afr' , NUMBER_OF_SIMULATED_training_INDIVIDUALS, allele_dictionary, compiled_samples_list)
			allele_dictionary, compiled_samples_list = pca_subPopulation_simulation( 'training', variant, var_id, 'AMR', 'aaf_gnomad_amr' , NUMBER_OF_SIMULATED_training_INDIVIDUALS, allele_dictionary, compiled_samples_list)
			allele_dictionary, compiled_samples_list = pca_subPopulation_simulation( 'training', variant, var_id, 'SAS', 'aaf_gnomad_sas' , NUMBER_OF_SIMULATED_training_INDIVIDUALS, allele_dictionary, compiled_samples_list)
			allele_dictionary, compiled_samples_list = pca_subPopulation_simulation( 'training', variant, var_id, 'EAS', 'aaf_gnomad_eas' , NUMBER_OF_SIMULATED_training_INDIVIDUALS, allele_dictionary, compiled_samples_list)

			# SIMULATE variant genotypes for various gnomad_subpopulation groups for PCA VALIDATION SET.

			allele_dictionary, compiled_samples_list = pca_subPopulation_simulation( 'validation', variant, var_id, 'NFE', 'aaf_gnomad_nfe' , NUMBER_OF_SIMULATED_validation_INDIVIDUALS, allele_dictionary, compiled_samples_list)
			allele_dictionary, compiled_samples_list = pca_subPopulation_simulation( 'validation', variant, var_id, 'AFR', 'aaf_gnomad_afr' , NUMBER_OF_SIMULATED_validation_INDIVIDUALS, allele_dictionary, compiled_samples_list)
			allele_dictionary, compiled_samples_list = pca_subPopulation_simulation( 'validation', variant, var_id, 'AMR', 'aaf_gnomad_amr' , NUMBER_OF_SIMULATED_validation_INDIVIDUALS, allele_dictionary, compiled_samples_list)
			allele_dictionary, compiled_samples_list = pca_subPopulation_simulation( 'validation', variant, var_id, 'SAS', 'aaf_gnomad_sas' , NUMBER_OF_SIMULATED_validation_INDIVIDUALS, allele_dictionary, compiled_samples_list)
			allele_dictionary, compiled_samples_list = pca_subPopulation_simulation( 'validation', variant, var_id, 'EAS', 'aaf_gnomad_eas' , NUMBER_OF_SIMULATED_validation_INDIVIDUALS, allele_dictionary, compiled_samples_list)


			# collect genotype information from HS targeted capture cohort.

			for individual in variant.samples_list:

				if individual in compiled_samples_list:
					None

				else:
					compiled_samples_list.append(individual)

				sample_alt_count = variant.patientAltCount(individual)
				allele_dictionary[var_id][individual] = sample_alt_count


		##
		## Return simulated and observed variant genotype information in PCA-friendly-table.
		## First column will be sample names.  Each other column is genotype (0,1, or 2) for an
		## individual... ie, each row is a patient or simulated individuals genotypes for all positions.
		##

		variants_list = allele_dictionary.keys()

		print 'number of individuals:\t', len(compiled_samples_list)

		out_list = []
		header = []
		header.append('SAMPLE')

		for x in variants_list:

			header.append(x)

		out_list.append(header)
		individual_counter = 0

		for individual in compiled_samples_list:

			individual_counter += 1
			tmp = []
			tmp.append(individual)

			for var in variants_list:

				tmp.append(allele_dictionary[var][individual])

			out_list.append(tmp)

		return(out_list)

class Exon:
	def __init__ (self,exonVar):
		self.chrom = exonVar[0]
		self.tx_id = exonVar[1]
		self.exon_start = int(exonVar[2])
		self.exon_end = int(exonVar[3])
		try:
			self.coding_start = int(exonVar[4])
		except:
			self.coding_start = 'None'
		try:
			self.coding_end = int(exonVar[5])
		except:
			self.coding_end = 'None'
		self.strand = int(exonVar[6])
		self.gene_name = exonVar[7]
		self.seq = exonVar[8]

##################
# Constant setup #
##################
codon_dictionary = {}

codon_dictionary['TTT'] = 'F' #'Phenylalanine'
codon_dictionary['TTC'] = 'F' #'Phenylalanine'
codon_dictionary['TTA'] = 'L' #'Leucine'
codon_dictionary['TTG'] = 'L' #'Leucine'

codon_dictionary['CTT'] = 'L' #'Leucine'
codon_dictionary['CTC'] = 'L' #'Leucine'
codon_dictionary['CTA'] = 'L' #'Leucine'
codon_dictionary['CTG'] = 'L' #'Leucine'

codon_dictionary['ATT'] = 'I' #'Isoleucine'
codon_dictionary['ATC'] = 'I' #'Isoleucine'
codon_dictionary['ATA'] = 'I' #'Isoleucine'
codon_dictionary['ATG'] = 'M' #'Methionine'

codon_dictionary['GTT'] = 'V' #'Valine'
codon_dictionary['GTC'] = 'V' #'Valine'
codon_dictionary['GTA'] = 'V' #'Valine'
codon_dictionary['GTG'] = 'V' #'Valine'

codon_dictionary['TCT'] = 'S' #'Serine'
codon_dictionary['TCC'] = 'S' #'Serine'
codon_dictionary['TCA'] = 'S' #'Serine'
codon_dictionary['TCG'] = 'S' #'Serine'

codon_dictionary['CCT'] = 'P' #'Proline'
codon_dictionary['CCC'] = 'P' #'Proline'
codon_dictionary['CCA'] = 'P' #'Proline'
codon_dictionary['CCG'] = 'P' #'Proline'

codon_dictionary['ACT'] = 'T' #'Threonine'
codon_dictionary['ACC'] = 'T' #'Threonine'
codon_dictionary['ACA'] = 'T' #'Threonine'
codon_dictionary['ACG'] = 'T' #'Threonine'

codon_dictionary['GCT'] = 'A' #'Alanine'
codon_dictionary['GCC'] = 'A' #'Alanine'
codon_dictionary['GCA'] = 'A' #'Alanine'
codon_dictionary['GCG'] = 'A' #'Alanine'

codon_dictionary['TAT'] = 'Y' #'Tyrosine'
codon_dictionary['TAC'] = 'Y' #'Tyrosine'
codon_dictionary['TAA'] = '*' #'STOP'
codon_dictionary['TAG'] = '*' #'STOP'

codon_dictionary['CAT'] = 'H' #'Histidine'
codon_dictionary['CAC'] = 'H' #'Histidine'
codon_dictionary['CAA'] = 'Q' #'Glutamine'
codon_dictionary['CAG'] = 'Q' #'Glutamine'

codon_dictionary['AAT'] = 'N' #'Asparagine'
codon_dictionary['AAC'] = 'N' #'Asparagine'
codon_dictionary['AAA'] = 'K' #'Lysine'
codon_dictionary['AAG'] = 'K' #'Lysine'

codon_dictionary['GAT'] = 'D' #'Aspartic_Acid'
codon_dictionary['GAC'] = 'D' #'Aspartic_Acid'
codon_dictionary['GAA'] = 'E' #'Glutamic_Acid'
codon_dictionary['GAG'] = 'E' #'Glutamic_Acid'

codon_dictionary['TGT'] = 'C' #'Cysteine'
codon_dictionary['TGC'] = 'C' #'Cysteine'
codon_dictionary['TGA'] = '*' #'STOP'
codon_dictionary['TGG'] = 'W' #'Tryptophan'

codon_dictionary['CGT'] = 'R' #'Arginine'
codon_dictionary['CGC'] = 'R' #'Arginine'
codon_dictionary['CGA'] = 'R' #'Arginine'
codon_dictionary['CGG'] = 'R' #'Arginine'

codon_dictionary['AGT'] = 'S' #'Serine'
codon_dictionary['AGC'] = 'S' #'Serine'
codon_dictionary['AGA'] = 'R' #'Arginine'
codon_dictionary['AGG'] = 'R' #'Arginine'

codon_dictionary['GGT'] = 'G' #'Glycine'
codon_dictionary['GGC'] = 'G' #'Glycine'
codon_dictionary['GGA'] = 'G' #'Glycine'
codon_dictionary['GGG'] = 'G' #'Glycine'

#############
# Functions #
#############
def pca_subPopulation_simulation(simulation_id, VARIANT, var_id, population_name, sub_population_code , NUMBER_OF_SIMULATED_INDIVIDUALS, allele_dictionary, compiled_samples_list):

	counter=0

	for x in xrange(0, NUMBER_OF_SIMULATED_INDIVIDUALS):
		counter += 1

		indv = '%s_%s_%i' % ( simulation_id , population_name , counter )

		if indv in compiled_samples_list:
			None
		else:
			compiled_samples_list.append(indv)

		indv_allele_count = 0

		##
		## Generate random numbers to simulate alleles.
		##

		indv_random_1 = random.random()
		indv_random_2 = random.random()

		##
		## Determine if alleles are 'alt' or 'ref' by whether the random number is '<' or '>='
		## gnomad subpopulation alt allele frequency.
		##

		if float(indv_random_1) < float(VARIANT.info[sub_population_code]):
			indv_allele_count += 1
		else:
			None
		if float(indv_random_2) < float(VARIANT.info[sub_population_code]):
			indv_allele_count += 1
		else:
			None

		allele_dictionary[var_id][indv] = indv_allele_count

	return (allele_dictionary, compiled_samples_list)

def extractVariants(FileName):

	##  just a simple function for converting a file of tab separated values into a list of lists
	##  each list is a line from the original file, each item in the list is a tab separated value

	InFileName = FileName
	varList = []
	
	with open( InFileName, 'r' ) as InFile:
		for Line in InFile:
			temp = Line.strip('\n')
			temp = temp.strip('\r')
			tempLine = temp.split('\t')
			cleanlist = []
			for x in tempLine:
				if x == '':
					cleanlist.append('None')
				else:
					cleanlist.append(x)
			varList.append(cleanlist)
			
	return (varList)

def outPut_NODATE(varList, fileNameNoDate):
	OutFileName = fileNameNoDate
	
	with open( OutFileName, 'a' ) as OutFile:
		for x in varList:
			Line = ''
			elementCounter = 1
			for y in x:
				if elementCounter == len(x):
					Line = Line+str(y)+'\n'
				else:
					Line = Line+str(y)+'\t'
					elementCounter += 1
			OutFile.write(Line)
	return

def nucleotide_comp(base):
	if base.upper() == 'A':
		return('T')
	elif base.upper() == 'T':
		return('A')
	elif base.upper() == 'G':
		return('C')
	elif base.upper() == 'C':
		return('G')

def GenerateTranscriptGenomicPositionDictionaries(transcript_info_list, gene_to_analyze, REGULATORY_BUFFER_WINDOW):

	#  This code is responsible for parsing and storing transcript/protein sequence information in order to later evaluate
	#  whether a variant is in a particular gene region, whether it is in a coding region, etc. The dictionaries produced
	#  from this function also allow one to ask what effect a variant has on protein sequence.
	#
	#  Currently this function requires a list of exon sequences and information about their genomic positions, strandedness, etc.
	#  See EXON class to see what information is used.  Also, this is currently set up to require the exon information be
	#  provided in a particular order.  This could be changed to require the information in any order but with particular
	#  field headers so as to make it easy to parse what information is in what column.

	exon_positions_dictionary = {}         # exon sites...     keys = genomic_position in chr#_position form ex.) 'chr5_1252342'; value = nucleotide  }
	coding_positions_dictionary = {}       # coding sites      same as exon_positions_dictionary but only contains coding positions
	gene_position_dictionary = {}          # key = chr_position ; value = gene for which that chr_position is associated with

	#################################################
	#################################################
	## ESTABLISH NUCLEOTIDE/CHR:POS SEQUENCE FOR GENE
	##

	# Read in exon information

	for x in transcript_info_list[1:]:
		curr_exon = Exon(x)

		if curr_exon.gene_name == gene_to_analyze:
			strandedness = curr_exon.strand

			# First sets up a list of all exonic genomic positions for the gene of interest.

			exon_positions_list = range(curr_exon.exon_start,curr_exon.exon_end+1)
			directional_exon_positions_list=[]
			if curr_exon.strand == 1:
				for nt in exon_positions_list:
					directional_exon_positions_list.append(nt)
			elif curr_exon.strand == -1:
				for nt in reversed(exon_positions_list):
					directional_exon_positions_list.append(nt)

			# Create a dictionary to be able to access any nucleotide value for any exonic genomic position.
			# Keys are established as ex) 'chr5_1252342'.  Originally keys were set up so that one could create
			# one big dictionary that contained information for all transcripts.  Having the chromosome identifier
			# in the key prevented numeric positional collisions between chromosomes.  However, I think that I ended
			# up just analyzing one gene at a time.  So I don't know that using the chromosome in the key is
			# actually essential to the process.
			#
			# *** MUST EVENTUALLY KEEP STRANDEDNESS/DIRECTIONALITY IN MIND.  Biomart data returns exon sequences
			#     that are in correct strand/order for direct translation.  So for a negative stranded gene that is encoded
			#     between genomic positions 100-500, the exon base sequence string will be reverse complement. (i.e., exon
			#     sequence nucleotide position 1 is the complement the positive strand nucleotide at genomic position 500).
			#     HOWEVER, when you are evaluating a variant on a negative stranded gene, you have to keep in mind that
			#     the genotype caller is only calling variants with respect to the positive strand of the reference genome.
			#     So for a variant from the HS_CAPTURE (or the gnomAD data sets), if the variant is ref: A , alt: G at genomic
			#     position 400, this is actually a T --> C variant at position 100 in the mRNA.
			#

			for (position,nucleotide) in zip(directional_exon_positions_list,curr_exon.seq):
				chrm_pos = 'chr%s_%i' % (curr_exon.chrom,position)
				exon_positions_dictionary[chrm_pos] = nucleotide.upper()

			# Create a dictionary as above, but now for coding positions.
			# Key = chrom/position ex) 'chr5_1252342' --> value = A (nucleotide)

			if curr_exon.coding_start != 'None' and curr_exon.coding_end != 'None':
				for cod_position in xrange(curr_exon.coding_start,curr_exon.coding_end+1):
					cod_chrm_pos = 'chr%s_%i' % (curr_exon.chrom,cod_position)
					coding_positions_dictionary[cod_chrm_pos] = exon_positions_dictionary[cod_chrm_pos]

			else:
				None

		else:
			None

	# This creates a broader 'positions list' and figures out the maximum and minimum genomic positions that encompass the exonic region
	# of the gene. The idea here is to be able to adjust the window that categorizes a variant as being 'associated' with or near a gene.
	# after establishing the exonic max and min positions, a user input for regulatory buffer window can be used to create new bigger region
	# around the exonic coding region to serve as an area of interest.  So if you wanted to include variants that are within 500 bases
	# of the exonic start/end of NCSTN as being "NCSTN associated variants", you can do that.

	pos_list = [int(each.split('_')[1]) for each in exon_positions_dictionary.keys()]
	pos_list.sort()
	min_pos = min(pos_list) - REGULATORY_BUFFER_WINDOW
	max_pos = max(pos_list) + REGULATORY_BUFFER_WINDOW

	current_chrom = exon_positions_dictionary.keys()[0].split('_')[0]

	# This then provides the mechanism for categorizing positions as associated with genes.
	# For this gene_position_dictionary:  key = ex.) 'chr5_1252342' --> value = 'NCSTN'
	# a

	for y in xrange(min_pos,max_pos +1):
		chrm_pos = '%s_%i' % (current_chrom, y)
		gene_position_dictionary[chrm_pos] = gene_to_analyze

	###########################################
	###########################################
	## ESTABLISH AA_SEQUENCE FOR PROTEIN GENE
	## (FACTOR IN STRANDEDNESS OF ENCODED GENE)

	codon_positions_dictionary = {}        # {keys = chr_position (ex. 'chr5_1252342') ; values = [[list of the three chr_position's that make up the codon to which the 'key' chr_position belongs], amino acid number within protein]}
	                                       # ex) transcript: chr_1, positions 1-300, strand = + ;  codon_positions_dictionary[chr1_1] = [[chr1_1 , chr1_2 , chr1_3 ], 1]... or  codon_positions_dictionary[chr1_8] =  [[chr1_7 , chr1_8 , chr1_9 ], 3]

	# Extract coding positions into a list and then order list depending on transcript strandedness.

	coding_positions_list = [int(each.split('_')[1]) for each in coding_positions_dictionary.keys()]
	coding_positions_list.sort()
	directional_coding_list = []

	if strandedness == 1:
		for base in coding_positions_list:
			directional_coding_list.append('%s_%i' % (current_chrom,base) )
	elif strandedness == -1:
		for base in reversed(coding_positions_list):
			directional_coding_list.append('%s_%i' % (current_chrom,base) )

	# Loop through coding positions collecting 3 positions (and then 3 nucleotides by querying what nucleotide is at each given position)
	# Fill in codon_positions_dictionary described above.
	# Codons are translated using a codon dictionary that I set up in the actual code, outside of this function.
	# Code is written to determine full amino acid sequence that is produced by the transcript.  This was used primarily for
	#      trouble shooting purposes.  The code was left in to make it easy to go back and trigger it to print out amino
	#      acid sequence string if desired.

	codon_string = ''
	aa_sequence = ''
	aa_counter = 1
	codon_positions_list = []
	codon_positions_dictionary = {}
	gene_sequence = ''

	for base in directional_coding_list:

		codon_string += exon_positions_dictionary[base]
		codon_positions_list.append(base)

		if len(codon_string)>0 and len(codon_string)%3 == 0:
			for nuc in codon_positions_list:
				codon_positions_dictionary[nuc] = [codon_positions_list , aa_counter]

			aa_sequence+= codon_dictionary[codon_string]
			aa_counter += 1
			gene_sequence += codon_string
			codon_string = ''
			codon_positions_list = []
		else:
			None

	#print aa_sequence

	return(exon_positions_dictionary,coding_positions_dictionary,gene_position_dictionary, codon_positions_dictionary, strandedness)

def evaluateVariants( TRACKER_dictionary, VARIANT_dictionary, exon_positions_dictionary, coding_positions_dictionary, gene_position_dictionary, codon_positions_dictionary, strandedness):

	# This function takes in the output of GenerateTranscriptGenomicPositionDictionaries(). Given a variant from the targeted capture
	#    this function will determine if the variant is in the gene being analyzed, determine if it is a non-coding/coding variant, evaluate the
	#    effect on amino acid sequence if it is a coding variant, and sort variants into lists of synonymous and missense variants.

	# **** READ ****
	#
	# You cannot simply look at position of variant to determine if it is coding.
	# For example if there is a 3 base multibase substitution (ie. ATG --> GTA) or deletion (ATGA --> T),
	# The variant position is reported as the position of reference A.
	# If A is in an intron, but TG is in an exon, then this is a coding variant.
	# Currently this code does not make use of intron sequences, so handling exon/intron border
	# straddling variants cant be evaluated for what bases would take the place of the deleted variants.
	# Even had we incorporated intron sequences into the analysis, one cannot say what the coding consequences would be.
	# Odds are that splicing would be dramatically affected, so I don't think that analyzing impact of the
	# deletion variant would be as simple as retranslating the protein with the bases that fall into the
	# positions of the deleted bases.  Right now I am collecting these variants separately for manual curation.

	## *** as described in transcript function above...
	##     MUST KEEP STRANDEDNESS/DIRECTIONALITY IN MIND.  Biomart data returns exon sequences
	##     that are in correct strand/order for direct translation.  So for a negative stranded gene that is encoded
	##     between genomic positions 100-500, the exon sequence string will be reverse complement. (i.e., exon
	##     sequence nucleotide position 1 is the complement the positive strand nucleotide at genomic position 500).
	##
	##     HOWEVER, the genotype callers only call variants with respect to the positive strand of the reference genome.
	##     For a variant from the HS_CAPTURE (or the gnomAD data sets), if the variant is ref = A , alt = G at genomic
	##     position 400, this is actually a T --> C variant at position 100 in the mRNA for a negative stranded gene.
	##     So if you were to evaluate the new ALT codon in a negative stranded gene produced by a point substitution
	##     variant, the 2 bases from the biomart data are already correct complementarity, but the VCF variant data is not
	##     and you need to use the complementary base to what is reported in the VCF.
	#

	# Establish which chromosome the current analysis is taking place on.
	# Kind of a silly way that I set this up... might be a smarter cleaner way of doing it.
	# Basically just sample one key from codon_positions_dictionary and extract the chromosome.

	key_counter = 0
	for key in coding_positions_dictionary.keys():
		if key_counter > 0:
			break
		else:
			CURRENT_CHROM_BEING_ANALYZED = key.split('_')[0].strip('chr')
		key_counter += 1

	# Set up lists for output...

	missense_coding_vars = []
	synonymous_coding_vars = []
	manual_curation_coding_vars = []
	non_coding_vars = []
	off_target_vars = []

	# Start looping through variants on the chromosome of current interest...

	for tempVariant in VARIANT_dictionary[CURRENT_CHROM_BEING_ANALYZED]:

		# Set or re-set switches/variables...

		mutation = 'None'             # I originally just set things up to refer to everything as 'mutation'
		                              # which obviously isn't accurate; might consider doing find/replace if we
									  # being genetically precise about coding variable names...  I also use 'mutPos','mutChrom', etc. below
		curr_gene = 'None'
		isOnTargetGene = False
		isFrameshift = False
		isIndel = False
		isExonic = False
		isNoncoding = False
		needsManualCurating = False

		mutChrom = tempVariant.info['CHROM']
		mutPos = int(tempVariant.info['POS'])

		mutPosition = 'chr%s_%i' % ( mutChrom , mutPos )
		mutRef = tempVariant.info['REF']
		mutAlt = tempVariant.info['ALT']

		####
		####  Double-check to see if variant is even on the same chromosome as gene being analyzed; break if not;
		####  If everything is ok, proceed with analysis...

		if mutChrom != CURRENT_CHROM_BEING_ANALYZED:
			None
		else:

			#######################################
			#######################################
			##  Variant   AA_CHANGE  Evaluation  ##
			#######################################
			##
			## Use boolean switches from above to classify what type of effect the variant will have on the gene.

			isPointMutation , isMultiBaseSubstitution , isIndel , isFrameshift , needsManualCurating, isNoncoding , isOFFtarget= variantCategorizer( tempVariant , strandedness , exon_positions_dictionary , coding_positions_dictionary, codon_positions_dictionary , codon_dictionary , gene_position_dictionary)

			#  Figure out coding impact of variant.  If variant is in coding region, put the variant
			#  through a function to determine if the variant is a single base substitution, an indel,
			#  or a frameshift variant.  Then determine functional impact on coding region.

			if isIndel or isFrameshift:

				curr_gene , mutation = evaluate_Frameshift_Indel_Mutation(tempVariant, exon_positions_dictionary, coding_positions_dictionary, gene_position_dictionary, codon_positions_dictionary, strandedness)

			elif isPointMutation or isMultiBaseSubstitution:

				curr_gene , mutation = evaluate_Substitution_Mutation(tempVariant, exon_positions_dictionary, coding_positions_dictionary, gene_position_dictionary, codon_positions_dictionary, strandedness)

			## cycle through genomic positions of ref/alt bases to determine which gene (if any) the non-coding or manual_inspection variant resides in
			##
			## 10-31-2019  I have tested this code on some previously annotated indel/substitution variants
			## but I don't yet know of any straddling variants to use as test cases.  So need to watch this when
			## I run on bigger variant sets to see if it is working.

			elif needsManualCurating or isNoncoding:

				curr_gene = manualCuration_nonCoding_Handling( tempVariant , strandedness , exon_positions_dictionary , coding_positions_dictionary, codon_positions_dictionary , codon_dictionary , gene_position_dictionary )

				if needsManualCurating:

					mutation = 'Manual_Inspection'

				elif isNoncoding:

					mutation = 'None'

			elif isOFFtarget:

				mutation = 'None'
				curr_gene = 'None'

			##
			## Add gene/mutation information to VARIANT
			##

			try:
				tempVariant.info['GENE']

				if curr_gene == 'None':
					None

				else:

					# When a variant is re-assessed within a domain of a gene that it has already been analyzed for,
					# the mutation AA number is evaluated with respect to the domain.  i.e., if a gene has 500 amino acids
					# and the domain is amino acids 350-450, the current evaluator tries to annotate the variant
					# as GENE: PSTPIP1/PSTPIP1_DOMAIN ; AA_CHANGE: L400R/L50R.   So currently I am handling this below by trying to
					# evaluate if the current analysis is for a domain and the duplicating the AA_CHANGE with respect to the entire gene.
					# So it should report as GENE: PSTPIP1/PSTPIP1_DOMAIN ; AA_CHANGE: L400R/L400R.

					# The reason this is set up this way is to be able to track overlapping genes.  So if GENE_Y were a gene on the opposite
					# strand as PSTPIP1, there could be a variant that then gets reported as:
					# GENE: PSTPIP1/PSTPIP1_DOMAIN/GENE_Y ; AA_CHANGE: L400R/L400R/Y23W

					if tempVariant.info['GENE'] in curr_gene:                                  ## *** NOTE: this really only works if domains are listed after full genes in the transcript file. Ought to tidy this up at some point.
						tempVariant.appendInfo('AA_CHANGE',tempVariant.info['AA_CHANGE'])

					else:
						tempVariant.appendInfo('AA_CHANGE', mutation)

					tempVariant.appendInfo('GENE',curr_gene)
			except:

				if curr_gene == 'None':
					None

				else:
					tempVariant.appendInfo('GENE',curr_gene)
					tempVariant.appendInfo('AA_CHANGE',mutation)

				######################################
				######################################
				#####  Variant Return Handling
				######################################

				##
				##  Sort variants into 'No mutation', 'Synonymous', or 'Missense' variant bins.
				##  Indels/frameshifts/aa_changing_substitutions are all filtered to missense.
				##
				##  TRACKER_dictionary prevents duplication of variants in output list.  This could happen from having genes
				##  that overlap genomic positions, or by analyzing a gene domain after analyzing the full gene.  If a variant
				##  is in two genes, this is tracked by having multiple gene names in the variant.info['GENE'] information slot.
				##  But I didn't want to have the same exact variant appended to the output missense list twice.  So the TRACKER_dictionary
				##  just keeps track of what variant ID's have been put into the output lists so far.  Then if a the code analyzes
				##  another missense variant, it checks to see if that variant has already has already been analyzed before adding
				##  it to the output list.
				##

				var_id = '%s_%s_%s_%s' % ( tempVariant.info['CHROM'] , tempVariant.info['POS'] , tempVariant.info['REF'] , tempVariant.info['ALT'])

				''' TROUBLESHOOTING
				print 'variant_id:\t', var_id
				print 'mutation:\t', mutation
				print 'aa_change_list:\t', aa_change_list
				'''

				if mutation == 'None' and curr_gene != 'None':
					if var_id in TRACKER_dictionary['noncoding']:
						None
					else:
						TRACKER_dictionary['noncoding'][var_id] = 1
						non_coding_vars.append(tempVariant)

				elif mutation != 'None' and mutation != 'Manual_Inspection':

					aa_change_list = re.split( '[\d]+', mutation )

					if aa_change_list[0] == aa_change_list[1]:

						if var_id in TRACKER_dictionary['synonymous']:
							None
						else:
							TRACKER_dictionary['synonymous'][var_id] = 1
							synonymous_coding_vars.append(tempVariant)

					elif aa_change_list[0] != aa_change_list[1]:

						if var_id in TRACKER_dictionary['missense']:
							None
						else:
							TRACKER_dictionary['missense'][var_id] = 1
							missense_coding_vars.append(tempVariant)

				elif mutation == 'Manual_Inspection':
					if var_id in TRACKER_dictionary['manual_inspection']:
						None
					else:
						TRACKER_dictionary['manual_inspection'][var_id] = 1
						manual_curation_coding_vars.append(tempVariant)

				else:
					# variant is not on target for this gene; cant categorize as OFF TARGET though because it might be
					# ON TARGET for a different gene, since we evaluate genes one at a time... if we were to analyze genes
					# simultaneously, then we could create an OFF TARGET collection here

					# print 'Variant was not placed into one of the following groups: missense, synonymous, noncoding, off_target.'
					None
			else:
				None

	return( missense_coding_vars , synonymous_coding_vars , manual_curation_coding_vars, non_coding_vars )

def variantCategorizer( current_variant , strandedness , exon_positions_dictionary , coding_positions_dictionary, codon_positions_dictionary , codon_dictionary , gene_position_dictionary ):

	# Is it an indel, or a frameshift, or a point mutation, or multi_base_substitution ??
	# This code returns a list of booleans to identify what type of variant it is.
	#
	# Currently it is not flagging anything more descriptive than 'isNoncoding', if it does not
	# affect a coding region.  We could change this to try to keep track of intronic vs intergenic,
	# or to catch variants that are a particular distance from a splice site.
	#
	# Categorization logic:
	#
	# len(REF) == len(ALT) :  if any position in alt/ref is coding, this will be either point mutation or multibase substitution
	# len(REF) != len(ALT)
	#       if:   ABS( ( len(REF) - len(ALT) ) ) % 3 == 0 :  amino acid insertion or deletion   ***Must check for exon border straddling case --> funnel to manual_curation_variants
	#       if:   ABS( ( len(REF) - len(ALT) ) ) % 3 != 0 :  frameshift variant                 ***Must check for exon border straddling case --> funnel to manual_curation_variants


	# set/re-set return boolean switches

	isPointMutation = False
	isMultiBaseSubstitution = False
	isIndel = False
	isFrameshift = False
	needsManualCurating = False
	isNoncoding = False
	isOFFtarget = False

	STRAND = strandedness

	current_chrom = current_variant.info['CHROM']
	var_pos = int(current_variant.info['POS'])
	chrom_position = 'chr%s_%i' % (current_chrom,var_pos)

	ref_bases = current_variant.info['REF']
	alt_bases = current_variant.info['ALT']


	# Collect all coding positions for gene into a list.  Collect just integer value from 'chr5_1252342' chrm_position.

	coding_positions_list = []

	for each_key in codon_positions_dictionary.keys():
		coding_positions_list.append(int(each_key.split('_')[1]))


	###
	### Handle substitutions
	###

	if len(ref_bases) == len(alt_bases):

		if len(ref_bases) == 1:

			if chrom_position in coding_positions_dictionary:
				isPointMutation = True

			else:
				None												  #  isNONcodingPointMutation = True

		elif len(ref_bases) > 1:

			##
			## Exon border testing. Determine how many of ref/alt genomic positions are coding positions.
			## If only one base of 3 base substitution is coding, then this is really just a point mutation... and
			## actually is more likely a splice affecting variant, or affects START or STOP codons, but it will just
			## be tracked as a point mutation.

			var_positions_list = []
			temp_var_pos = var_pos

			for each_base in ref_bases:

				curr_chr_pos = 'chr%s_%i' % ( current_chrom , temp_var_pos )
				var_positions_list.append(curr_chr_pos)

				temp_var_pos += 1


			var_coding_positions_counter = 0

			for each_var_pos in var_positions_list:

				if each_var_pos in coding_positions_dictionary:
					var_coding_positions_counter += 1

				else:
					None

			if var_coding_positions_counter == 0:
				None                                                  #  isNONcodingSubstitution = True

			elif var_coding_positions_counter == 1:
				isPointMutation = True                                #  might be better to funnel to manual curation

			elif var_coding_positions_counter > 1:
				isMultiBaseSubstitution = True


	###
	### Handle insertions.  Insertions are easier to handle than deletions.
	###

	elif len(ref_bases) < len(alt_bases):

	#                                    1     12345
	# 3 base insertion is referenced as: T --> TCAG  in a VCF.
	# Another way to look at it is:
	#
	#            123     1abc23
	#            TXX --> TCAGXX
	#
	# For this type of variant we are interested in whether T and (T + 1), positions #1 and #2 are in a coding region.
	# If they are both coding, then the inserted CAG is also coding.  If only one of the positions is coding, then the
	# insertion is intron adjacent... it will likely affect splicing, but theoretically not a coding variant... filter
	# these intron adjacent insertions for manual inspection.

		left_boundary_pos = var_pos
		right_boundary_pos = var_pos + 1

		left_chrm_pos  = 'chr%s_%i' % ( current_chrom ,  left_boundary_pos )
		right_chrm_pos = 'chr%s_%i' % ( current_chrom , right_boundary_pos )

		if left_chrm_pos in coding_positions_dictionary and right_chrm_pos in coding_positions_dictionary:

			if abs( len(ref_bases) - len(alt_bases)) % 3 == 0:

				isIndel = True      # INSERTION

			elif abs( len(ref_bases) - len(alt_bases)) % 3 != 0:
				isFrameshift = True

		elif left_chrm_pos in coding_positions_dictionary or right_chrm_pos in coding_positions_dictionary:

			needsManualCurating = True

		else:

			None                                                     #  isNONcodingINDEL = True


	###
	###  Handle deletions
	###

	elif len(ref_bases) > len(alt_bases):

		# LOGIC:  deletions are trickier...
		#                         1234     1      <-- position
		#   example variant #1:   TCAG --> T      <-- nucleotide
		#
		#		  ccccccc      cccc     c = coding position
		#         1234abc      1abc
		#         TCAGTXX  --> TTXX
		#
		# If a deletion is flanked, as above, by genomic positions that are coding (i.e. position #1 and position 'a'), then the deletion is a simple coding deletion.
		# In the case of an extremely big deletion, you could have a situation where the deletion is flanked by coding positions, but those coding
		# positions are in different exons.  We can check for this by asking whether every genomic position between the flanking site is a coding position; if
		# not then the deletion is spanning an intron.  These variants may actually be easy to analyze amino_acid change because there is no splicing at the
		# deleted intron to worry about... but I will just filter these variants for manual inspection.
		#
		#   example variant #2:
		#
		#         cccnnnn      cnnn    <--  c = coding position ; n = non-coding position
		#         1234abc      1abc    <--  position
		#         TCAGTXX  --> TTXX    <--  nucleotide
		#
		# If a deletion is not flanked by coding positions, that does not necessarily mean that it is not a coding deletion. In the
		# example above, base at position 'a' is not a coding position, but some of the deleted bases were in coding bases.  This kind of variant
		# will likely be a splice affector and there will be no easy way to assess impact on amino acid sequence.  To find these variants,
		# we need to ask whether any of the deleted positions are coding positions. If there are deleted coding positions, but they
		# are not flanked by coding positions, then filter variant to inspect manually.
		# I think that this strategy should also catch variants where both flanking positions are non-coding, but the deletion spans an entire exon.
		# In that case all of the positions in the deletion would be sampled, we would detect a coding site, and it sould filter for manual inspection.
		#

		left_boundary_pos = var_pos
		right_boundary_pos = ( var_pos + (abs( len(ref_bases) - len(alt_bases)) + 1 ) )

		left_chrm_pos  = 'chr%s_%i' % ( current_chrom , left_boundary_pos )
		right_chrm_pos = 'chr%s_%i' % ( current_chrom , right_boundary_pos )

		#
		# Handle cases where deletion-flanking positions are coding positions.
		#

		if left_chrm_pos in coding_positions_dictionary and right_chrm_pos in coding_positions_dictionary:

			# Check that all genomic positions between deletion-flanking positions are coding positions.

			all_positions_are_contiguous = True
			for all_positions in xrange( left_boundary_pos , (right_boundary_pos+1) ):
				tmp_chrm_pos = 'chr%s_%i' % ( current_chrom , all_positions )

				if tmp_chrm_pos in coding_positions_dictionary:
					None

				else:
					all_positions_are_contiguous = False

			if all_positions_are_contiguous is False:
				needsManualCurating = True

			else:

				if abs( len(ref_bases) - len(alt_bases)) % 3 == 0:
					isIndel = True                                     # DELETION

				else:
					isFrameshift = True

		#
		# Handle cases where NOT BOTH of the deletion-flanking positions are coding positions.
		#

		else:

			var_positions_list = []
			current_var_pos = var_pos

			base_counter = 1
			for each_position in ref_bases:

				temp_position = current_var_pos + base_counter
				temp_chr_pos = 'chr%s_%i' % ( current_chrom , temp_position )
				var_positions_list.append( temp_chr_pos )

				base_counter += 1

			deletion_coding_positions_counter = 0

			for each_del_pos in var_positions_list:

				if each_del_pos in coding_positions_dictionary:


					deletion_coding_positions_counter += 1

				else:
					None

			if deletion_coding_positions_counter > 0:

				needsManualCurating = True

			else:

				None                                        # isNONcodingINDEL = True

	# If variant has not fallen into any of our categories, then classify it as Noncoding

	variant_bool_list = [ isPointMutation , isMultiBaseSubstitution , isIndel , isFrameshift , needsManualCurating ]

	if any(variant_bool_list):
		isNoncoding = False
		isOFFtarget = False

	else:

		exonic_positions_list = [int(x.split('_')[1]) for x in exon_positions_dictionary ]
		min_exon_position = min( exonic_positions_list )
		max_exon_position = max( exonic_positions_list )

		count_tracker = 0
		current_variant_positions_list = []
		if len(ref_bases) > len(alt_bases):
			for each_position in ref_bases:
				current_variant_positions_list.append( var_pos + count_tracker )
				count_tracker += 1

		else:
			for each_position in alt_bases:
				current_variant_positions_list.append( var_pos + count_tracker )
				count_tracker += 1

		if any( [x >= min_exon_position and x <= max_exon_position for x in current_variant_positions_list ]):
			isNoncoding = True
			isOFFtarget = False

		else:
			isNoncoding = False
			isOFFtarget = True


	variant_bool_list.append( isNoncoding )
	variant_bool_list.append( isOFFtarget )

	variant_classification_counts = variant_bool_list.count( True )

	if variant_classification_counts > 1:

		print 'oops, variant has been classified as more than one variant category: PointMutation, MultiBaseSubstitution, Indel, Frameshift, Manual_Curation, Noncoding'
	else:

		return ( variant_bool_list )

def manualCuration_nonCoding_Handling( current_variant , strandedness , exon_positions_dictionary , coding_positions_dictionary, codon_positions_dictionary , codon_dictionary , gene_position_dictionary ):
	current_chrom = current_variant.info['CHROM']
	var_pos = int(current_variant.info['POS'])
	chrom_position = 'chr%s_%i' % (current_chrom,var_pos)

	ref_bases = current_variant.info['REF']
	alt_bases = current_variant.info['ALT']

	# Establish exon min and max genomic positions
	exonic_positions_list = [int(x.split('_')[1]) for x in exon_positions_dictionary ]
	min_exon_position = min( exonic_positions_list )
	max_exon_position = max( exonic_positions_list )

	count_tracker = 0
	current_variant_positions_list = []
	if len(ref_bases) > len(alt_bases):
		for each_position in ref_bases:
			current_variant_positions_list.append( var_pos + count_tracker )
			count_tracker += 1
	else:
		for each_position in alt_bases:
			current_variant_positions_list.append( var_pos + count_tracker )
			count_tracker += 1

	if any( [x >= min_exon_position and x <= max_exon_position for x in current_variant_positions_list ]):

		current_gene = gene_position_dictionary[ exon_positions_dictionary.keys()[0] ]  # This just regurgitates whatever gene we are currently analyzing.

	else:
		current_gene = 'None'

	return( current_gene )

def evaluate_Frameshift_Indel_Mutation(current_variant, exon_positions_dictionary, coding_positions_dictionary, gene_position_dictionary, codon_positions_dictionary, strandedness):

	# Function to evaluate the amino acid sequence effect of frameshift or indel variants.
	#
	# Premise here is to collect all of the coding positions in sorted order, and convert those numeric positions to nucleotides.  All of this is done
	# using the different transcript sequence dictionaries that are created prior to variant analysis.  mRNA sequences are built both for the reference
	# transcript and the alternate transcript.  As the sequence mRNA nucleotide sequence is being built, the function tracks whether you have arrived
	# at the variant position, at which point the ALT variant is incorporated into the alternate mRNA sequence.  After mRNA sequences are built, the
	# sequences are compared and the amino acid differences are established.
	#
	# **** NOTE ****
	# Keep in mind that transcript dictionaries are mRNA complementarity for both + and - strand genes, but that the REF and ALT bases from the VCF VARIANT
	# are all + strand.
	#
	#

	#define some variables

	STRAND = strandedness

	current_chrom = current_variant.info['CHROM']
	var_pos = int(current_variant.info['POS'])
	chrom_position = 'chr%s_%i' % ( current_chrom , var_pos )

	current_eval_gene = gene_position_dictionary[chrom_position]

	# collect all coding positions

	coding_positions_list = []
	for each_key in codon_positions_dictionary.keys():
		coding_positions_list.append(int(each_key.split('_')[1]))

	#figure out how coding sequence is changed with variant
	ref_bases = current_variant.info["REF"]
	alt_bases = current_variant.info["ALT"]

	if len(ref_bases) > len(alt_bases):
		mut_type = 'deletion'
		base_difference = len(ref_bases) - len(alt_bases)

	elif len(alt_bases) > len(ref_bases):
		mut_type = 'insertion'
		base_difference = len(alt_bases) - len(ref_bases)

	#create new ALT transcript sequence
	new_alt_transcript = ''
	original_transcript = ''

	BASE_DIFF_COUNTER = base_difference
	START_VARIATION = False

	coding_positions_list.sort()       # Regardless of gene strandedness, it is easiest at this point to work from left to right on + strand.  This is
	                                   # because variants from VCF are all + strand order perspective.  So consider an insertion of G --> GTT.  For a negative
									   # strand mRNA, the 'TT' insertion after a 'G' is really an insertion of 'AA' before a 'C'.  If you work in the
									   # correct mRNA order for - strand gene, you would run across the C at genomic position 10, and then you would
                                       # have to go back and insert AA before the C, which is a pain to code for.  The Biomart transcript data is all
									   # correct strand.  So if you just work with - strand gene from the + strand perspective, as you build the mRNA
                                       # sequence, the nucleotides will be correct complementarity, but just in the reverse order.  This will allow you
                                       # to simply add the AA after the C.  Then when you want to evaluate amino acid sequence below, you can just
									   # new_sequence.sort(reverse = True) to flip the mRNA back to the correct order.
                                       #

	for current_position in coding_positions_list:

		temp_chr_pos = 'chr%s_%i' % ( current_chrom , current_position )
		current_reference_base = exon_positions_dictionary[temp_chr_pos]

		if current_position < var_pos:
			new_alt_transcript += current_reference_base
			original_transcript += current_reference_base

		elif current_position == var_pos:
			new_alt_transcript += current_reference_base
			original_transcript += current_reference_base

			if mut_type == 'insertion':

				for base in alt_bases[1:]: # variants reported as REF= A, ALT = ACCG; so what you need to actually insert is everything in ALT after the first base

					if STRAND == 1:
						new_alt_transcript += base

					elif STRAND == -1:
						new_alt_transcript += nucleotide_comp(base)   # Stored Biomart data is already in the proper strand.  Variant data is + strand, so need to add complementary base here.

		elif current_position > var_pos:

			if mut_type == 'deletion':

				if BASE_DIFF_COUNTER > 0:

					BASE_DIFF_COUNTER -= 1
					original_transcript += current_reference_base
					None

				else:

					new_alt_transcript += current_reference_base
					original_transcript += current_reference_base

			elif mut_type == 'insertion':

					new_alt_transcript += current_reference_base
					original_transcript += current_reference_base
	#
	# Rearrange transcripts to appropriate strand order.
	#

	if STRAND == -1:

		original_transcript = original_transcript[::-1]
		new_alt_transcript = new_alt_transcript[::-1]

	else:
		None

	#
	# translate original transcript
	#
	original_translation_string = ''
	temp_codon_string = ''
	for current_nucleotide in original_transcript:

		temp_codon_string += current_nucleotide

		if len(temp_codon_string) != 0 and len(temp_codon_string) % 3 == 0:

			temp_AA = codon_dictionary[temp_codon_string]
			original_translation_string += temp_AA
			temp_codon_string = ''

		else:
			None

	#translate new transcript until * (stop/nonsense) base is hit
	new_translation_string = ''
	temp_codon_string = ''

	for current_nucleotide in new_alt_transcript:

		temp_codon_string += current_nucleotide

		if len(temp_codon_string) != 0 and len(temp_codon_string) % 3 == 0:

			temp_AA = codon_dictionary[temp_codon_string]
			new_translation_string += temp_AA
			temp_codon_string = ''

		else:
			None

	# Determine how many coding positions are changed...
	#
	#
	# The variantCategorizer() should have filtered exon/intron boundary straddling variants for 'manual_inspection'.
	# So, all variants evaluated for coding effect that have reached this point, should be 'simple' frameshift/indels
	# where the genomic positions flanking the isertion/deletion are coding positions.
	#


	##
	## Determine codon insertion/deletion AA coding effect
	##
	##

	'''
	## We can't simply compare reference and variant translation strings here...
	##
	## edge case example:
	##
	##                                    *            <-- deleted amino acid
	## Reference aa sequence:  MLPRGAVADLAQQQQQRTPL
	## Variant aa sequence:    MLPRGAVADLAQQQQRTPL
	##                                        ^        <-- first difference in AA strings
	## If we just compared sequences until we found a difference in amino acid sequence, it might look like a Q-->R substitution, or we would
	## report the last Q as being deleted, when technically the first Q was deleted.  This might be irrelevant with regards to impact on protein,
	## but there is a difference with regards to variant genomic position.
	##


	#
	# Determine number of codons affected.  Note that this will depend on placement of the insert:
	#
	#       CODON #:       1   2   3   4   5   6
	#   Ref sequence:     ATG CTA ATA TTC GAT TCC
	#   Alt sequence 1:   ATG Cxx xTA ATA TTC GAT TCC     (TWO CODONS AFFECTED)
	#   Alt sequence 2:   ATG CTA xxx ATA TTC GAT TCC     (ONE CODON  AFFECTED)
	#
	#           REPORT:   [CODON1-AA]___[CODON1-#]___[AA-OF-ALL-AFFECTED-CODONS]
    #	#  Strategy:
	#     Determine # of codons affected. Report last reference AA, and then
	#     report subsequent AA of alt_AA_sequence.  Figure out how many codon positions
	#     the indel occupies. Then just read out alt_AA_sequence[last_ref_position + x]
	#     where x = xrange(1,(#_codons_occupied + 1))
	#
	'''

	## Update 4-30-20:  In the process of touching up figures/tables for manuscripts I needed to convert protein amino acid change nomenclature
	## to abide by HGVS nomenclature.  It turns out that my strategy described above is not what HGVS uses.  They reference amino acid change
	## at the most 3' position possible.  So in the case of polyQ deletions, even if the most 5' Q codon is affected with respect to the genomic
	## nucleotide variant, you would in fact slide the new protein sequence down in the alignment to the reference amino acid sequence, and then
	## call the 3' most positions as the deleted Q's.
	##
	##

	if len(ref_bases) > len(alt_bases):
		aa_change_type = 'deletion'
		ref_or_alt = ref_bases
	else:
		aa_change_type = 'insertion'
		ref_or_alt = alt_bases

	##
	##  Figure out INDEL AA impact.
	##

	if 	( len(ref_or_alt) -1 ) !=0 and ( ( len(ref_or_alt) -1 ) % 3 ) == 0:

		#
		# Verify that we are dealing with a 'simple' indel. Flanking positions should be coding.
		# mRNA length difference should be equal to (indel_base_#)/3.
		#
		left_flank_position = var_pos
		right_flank_position = ( var_pos + len(ref_or_alt) )

		LF_chrm_pos = 'chr%s_%i' % ( current_chrom , left_flank_position )
		RF_chrm_pos = 'chr%s_%i' % ( current_chrom , right_flank_position )

		if LF_chrm_pos in coding_positions_dictionary and RF_chrm_pos in coding_positions_dictionary:
			None

		else:
			sys.exit( 'Something went wrong:  Evaluating simple INDEL, but flanking genomic positions are NOT coding positions.')

		if abs( len(original_translation_string) - len(new_translation_string) ) == ( len(ref_or_alt) - 1)/3 :
			None

		else:
			print original_translation_string
			print 'original_protein_length:\t%s' % ( len(original_translation_string) )
			print new_translation_string
			print 'new_protein_length:\t%s' % ( len(new_translation_string) )
			print
			print '%s_%i_%s_%s' % (current_chrom , var_pos, ref_bases, alt_bases)
			sys.exit( 'Something went wrong:  Evaluating simple INDEL, but number of INDEL bases does not agree with aa_sequence lengths.')

		if aa_change_type == 'insertion':
			LONGER_AA_SEQUENCE = new_translation_string
			SHORTER_AA_SEQUENCE = original_translation_string
		elif aa_change_type == 'deletion':
			LONGER_AA_SEQUENCE = original_translation_string
			SHORTER_AA_SEQUENCE = new_translation_string

		number_of_indel_AAs = (len(ref_or_alt)-1)/3

		position_counter = 1
		temp_translation_string = SHORTER_AA_SEQUENCE[:position_counter]
		while temp_translation_string in LONGER_AA_SEQUENCE:
			position_counter += 1
			temp_translation_string = SHORTER_AA_SEQUENCE[:position_counter]

		temp_last_aa_position = position_counter - 1
		temp_last_ref_aa = LONGER_AA_SEQUENCE[position_counter -1 ]

		current_aa_position = temp_last_aa_position
		indel_string = temp_last_ref_aa
		for indel_AA in xrange(1,(number_of_indel_AAs+1)):
			current_aa_position += 1
			indel_string += LONGER_AA_SEQUENCE[current_aa_position - 1]

		last_indel_aa_position = current_aa_position

		''' TROUBLESHOOTING
		print temp_last_ref_aa
		print temp_last_aa_position
		print indel_string
		print
		print last_indel_aa_position
		print 'original_protein_length:\t%s' % ( len(original_translation_string) )
		print 'new_protein_length:\t%s' % ( len(new_translation_string) )
		'''

		if aa_change_type == 'insertion':

			indel_aa_change = '%s%i%s' % ( temp_last_ref_aa , temp_last_aa_position , indel_string )

		elif aa_change_type == 'deletion':

			indel_aa_change = '%s%i%s' % ( indel_string , temp_last_aa_position , temp_last_ref_aa )

		'''
		## This code was to annotate the indel by the number of positions affected by the nucleotide variant, not the HGVS nomenclature.
		#
		# Find the last codon NOT affected by insertion.
		#

		first_indel_base = left_flank_position + 1
		first_indel_pos = 'chr%s_%i' % ( current_chrom , first_indel_base )
		first_indel_codon = codon_positions_dictionary[ first_indel_pos ][1]
		last_reference_codon = first_indel_codon - 1

		#
		# Figure out how many coding position, indel bases occupy.
		#

		base_count_tracker = 1
		affected_positions_list = []
		affected_coding_positions_list = []
		affected_CODON_positions_list = []

		for x in ref_or_alt[1:]:
			affected_positions_list.append( var_pos + base_count_tracker )
			base_count_tracker += 1

		for x in affected_positions_list:
			tmp_pos_key = 'chr%s_%i' % ( current_chrom , x )

			if tmp_pos_key in coding_positions_dictionary:
				affected_coding_positions_list.append(x)

				if codon_positions_dictionary[ tmp_pos_key ] in affected_CODON_positions_list:
					None
				else:
					affected_CODON_positions_list.append( codon_positions_dictionary[ tmp_pos_key ] )

			else:
				None

		if aa_change_type == 'insertion':

			last_ref_aa = original_translation_string[ last_reference_codon - 1 ]
			insertion_aa_string = ''
			insertion_aa_string += last_ref_aa

			codon_position_tracker = last_reference_codon

			for each_codon in affected_CODON_positions_list :
				codon_position_tracker += 1
				insertion_aa_string += new_translation_string[ codon_position_tracker - 1 ]

			indel_aa_change = '%s%i%s' % ( last_ref_aa , last_reference_codon , insertion_aa_string )

		elif aa_change_type == 'deletion':

			last_ref_aa = original_translation_string[ last_reference_codon - 1 ]
			deletion_aa_string = ''
			deletion_aa_string += last_ref_aa

			codon_position_tracker = last_reference_codon

			for each_codon in affected_CODON_positions_list :
				codon_position_tracker += 1
				deletion_aa_string += original_translation_string[ codon_position_tracker - 1 ]

			indel_aa_change = '%s%i%s' % ( deletion_aa_string , last_reference_codon , last_ref_aa )

		'''

	else:
		##
		## Determine frameshift AA coding effect
		##
		aa_position_counter = 1

		for amino_acid in zip(original_translation_string , new_translation_string):

			if amino_acid[0] == amino_acid[1]:
				aa_position_counter += 1

			else:
				frameshift_aa_effect_position = aa_position_counter
				break

		last_ref_amino_acid = new_translation_string[ frameshift_aa_effect_position - 2]  # minus 2 here... subtract 1 for moving back to get last accurate amino acid... subtract 1 to put into base 0 for accessing correct amino acid from sequence string
		FRAMESHIFT_translation_string = ''

		for aa_position in new_translation_string[(frameshift_aa_effect_position - 2):]:
			FRAMESHIFT_translation_string += aa_position
			if aa_position == '*':
				break
			else:
				None


		# A frameshift close to the end of a protein might not result in a stop codon within the coding sequence of the mRNA.
		# The transcript information may lack the sequence required to assess frameshift translation past the coding region.
		# These variants will be labeled with '...' at the end of the frameshift sequence to denote that the frameshift would
		# continue translating.
		# ex) R415RFSGSDFSDF...


		if '*' in FRAMESHIFT_translation_string:
			None
		else:
			FRAMESHIFT_translation_string +='...'

		indel_aa_change = '%s%i%s' % ( last_ref_amino_acid , (frameshift_aa_effect_position - 1) , FRAMESHIFT_translation_string )

		#print frmshft_aa_change						# for troubleshooting

	return( current_eval_gene , indel_aa_change )

def evaluate_Substitution_Mutation(current_variant, exon_positions_dictionary, coding_positions_dictionary, gene_position_dictionary, codon_positions_dictionary, strandedness):

	# Function to evaluate the amino acid sequence effect of point mutations or multibase substitutions.
	#
	# Premise here is to collect all of the coding positions in sorted order, and convert those numeric positions to nucleotides.  All of this is done
	# using the different transcript sequence dictionaries that are created prior to variant analysis.  mRNA sequences are built both for the reference
	# transcript and the alternate transcript.  As the sequence mRNA nucleotide sequence is being built, the function tracks whether you have arrived
	# at the variant position, at which point the ALT variant is incorporated into the alternate mRNA sequence.  After mRNA sequences are built, the
	# sequences are compared and the amino acid differences are established.
	#
	# **** NOTE ****
	# Keep in mind that transcript dictionaries are mRNA complementarity for both + and - strand genes, but that the REF and ALT bases from the VCF VARIANT
	# are all + strand.
	#
	#
	#define some variables

	STRAND = strandedness

	current_chrom = current_variant.info['CHROM']
	var_pos = int(current_variant.info['POS'])
	chrom_position = 'chr%s_%i' % ( current_chrom , var_pos )

	current_eval_gene = gene_position_dictionary[chrom_position]

	# collect all coding positions

	coding_positions_list = []
	for each_key in codon_positions_dictionary.keys():
		coding_positions_list.append(int(each_key.split('_')[1]))

	#figure out how coding sequence is changed with variant
	ref_bases = current_variant.info["REF"]
	alt_bases = current_variant.info["ALT"]

	if len(ref_bases) == len(alt_bases) > 1:
		mut_type = 'multi_base_substitution'

	elif len(ref_bases) == len(alt_bases) == 1:
		mut_type = 'point_mutation'

	elif len(ref_bases) != len(alt_bases):
		print 'oops, this variant is labeled as substitution, but REF and ALT are different lengths.'


	#
	# Set up an ALT dictionary where key is numeric genomic position and value is the ALT nucletide at that position.
	# This will allow for us to loop through coding positions to build mRNA, and if we hit a value that is a key in
	# the ALT dictionary, we will insert the ALT base instead of the REF base.
	#

	pos_counter = 0
	alt_base_dictionary = {}

	for each_var_pos in	ref_bases:
		alt_base_dictionary[ (var_pos + pos_counter) ] = alt_bases[ pos_counter ]
		pos_counter += 1

	#create REF and ALT transcript sequences
	new_alt_transcript = ''
	original_transcript = ''

	coding_positions_list.sort()

	for current_position in coding_positions_list:

		temp_chrm_pos = 'chr%s_%i' % ( current_chrom , current_position )
		current_reference_base = exon_positions_dictionary[ temp_chrm_pos ]

		if current_position in alt_base_dictionary:

			if STRAND == 1:

				new_alt_transcript += alt_base_dictionary[ current_position ]
				original_transcript += current_reference_base

			elif STRAND == -1:

				new_alt_transcript += nucleotide_comp( alt_base_dictionary[ current_position ] )
				original_transcript += current_reference_base

		else:

			new_alt_transcript += current_reference_base
			original_transcript += current_reference_base

	#
	# Rearrange transcripts to appropriate strand order.
	#

	if STRAND == -1:

		original_transcript = original_transcript[::-1]
		new_alt_transcript = new_alt_transcript[::-1]

	else:

		None

	#
	# translate REF transcript
	#

	original_translation_string = ''
	temp_codon_string = ''
	for current_nucleotide in original_transcript:

		temp_codon_string += current_nucleotide

		if len(temp_codon_string) != 0 and len(temp_codon_string) % 3 == 0:

			temp_AA = codon_dictionary[temp_codon_string]
			original_translation_string += temp_AA
			temp_codon_string = ''

		else:
			None

	#
	# translate ALT transcript
	#

	new_translation_string = ''
	temp_codon_string = ''

	for current_nucleotide in new_alt_transcript:

		temp_codon_string += current_nucleotide

		if len(temp_codon_string) != 0 and len(temp_codon_string) % 3 == 0:

			temp_AA = codon_dictionary[temp_codon_string]
			new_translation_string += temp_AA
			temp_codon_string = ''

		else:
			None

	#
	# Identify amino acid change.
	#

	aa_position_counter = 1
	aa_change_switch = False
	ref_aa = ''
	alt_aa = ''
	affected_aa_positions_list = []

	for amino_acid in zip(original_translation_string , new_translation_string):

		if amino_acid[0] == amino_acid[1]:
			aa_position_counter += 1

		elif amino_acid[0] != amino_acid[1]:

			affected_aa_positions_list.append( aa_position_counter )
			aa_position_counter += 1

			ref_aa += amino_acid[0]
			alt_aa += amino_acid[1]

	## If amino acid sequences are identical (SYNONYMOUS VARIANT), just use position dictionaries to figure out which amino acids were 'synonymously affected'

	if len(affected_aa_positions_list) == 0:

		for each_key_position in alt_base_dictionary.keys():

			tmp_chrm_pos_key = 'chr%s_%i' % ( current_chrom , each_key_position )

			if tmp_chrm_pos_key in codon_positions_dictionary:

				tmp_aa_position = codon_positions_dictionary[tmp_chrm_pos_key][1]

				if tmp_aa_position in affected_aa_positions_list:
					None

				else:
					affected_aa_positions_list.append( tmp_aa_position )

			else:

				None

		reference_affected_aa_string = ''
		alt_affected_aa_string = ''
		affected_aa_positions_list.sort()

		for each_aa_position in affected_aa_positions_list:

			aa_position_index = ( each_aa_position - 1 )

			reference_affected_aa_string += original_translation_string[ aa_position_index ]
			alt_affected_aa_string += new_translation_string[ aa_position_index ]

		substitution_aa_change = '%s%i%s' % ( reference_affected_aa_string , min(affected_aa_positions_list) , alt_affected_aa_string )

	else:

		substitution_aa_change = '%s%i%s' % ( ref_aa , min(affected_aa_positions_list) , alt_aa )

	#
	# Report variant AA change.
	#

	return( current_eval_gene , substitution_aa_change )

if __name__ == "__main__":
	pass