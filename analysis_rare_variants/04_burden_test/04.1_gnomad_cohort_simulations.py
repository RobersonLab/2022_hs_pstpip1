#!/usr/bin/env python

##########
# Import #
##########
import numpy
import random
import sys

from burdens_code import *

###################################################################
###################################################################
##
##  04.1_GNOMAD_COHORT_SIMULATIONS.PY
##
##  This script performs ALT allele counts for targeted capture cohort, and it performs gnomAD individual simulations to predict
##  the expected ALT allele observations for a cohort of similar PCA ancestry to targeted capture cohort.
##
##  This requires that you perform a PCA ancestry analysis.  I did this in R.
##
##  The number of simulations can be changed.  As can the cohort make-up.  This takes in a "SAMPLES_FOR_ANALYSIS" file
##  which has sample_names and PCA_ancestry information.  This file is important because the capture was performed both on HS_AFFECTED
##  and NON_AFFECTED individuals, so the script needs to know which individuals to analyze and which to pass over.
##
##  Variant filtering strategy:
##  We wanted to use AAF cutoffs to perform multiple burdens analyses for different levels of variant rare-ness.
##  However, this becomes a problem when some variants have different AAF subpopulation frequencies that straddle
##  the cutoff frequency. (ex, cutoff_AAF = 0.01, variant: aaf_gnomad_nfe = 0.000001, aaf_gnomad_afr = 0.015).
##  The strategy I ended up choosing was to allow a variant to pass the cutoff filter if at least one subpopulation had an
##  AAF that was non-zero and < AAF_THRESHOLD.  The idea was that a rare variant could be rare and disease contributing in a
##  subpopulation, and that that same variant might be more common in another population and explain an increased risk for
##  disease within that subpopulation.
##
##  Potential gnomAD multi-allelic site handling:
##  gnomAD variant data is compiled from a variant_AAF list to genomic_position_AAF list.  The idea behind this was to deal with potential multi-allelic sites
##  Simulations are run by using a random number generator to select two values (between 0 and 1) to represent 'alleles'.  If a value/allele is less than the
##  gnomad_subpopulation_AAF being simulated, the allele is scored as ALT, otherwise it is scored as REF. We decided to pool the alt_allele_frequencies for all
##  variants at a genomic position, and then challenge that combined_alt_frequency during simulation.  REASON BEHIND DOING IT THIS WAY: Consider a genomic site
##  that is multiallelic, lets say three different variants (A->T, A->G, A->C).  If you didn't compile these, the simulator would challenge each variant.  Essentially
##  you would be drawing 6 random numbers, or 6 alleles for the individual.  Instead, I compiled the variant frequencies to figure out the total AAF for
##  any ALT (ie, A-> T or G or C).  Then the simulator chooses 2 alleles, and determines if an alt of any kind was observed.
##  POTENTIAL ISSUES:  I think this strategy should work because theoretically the frequencies of A, T, G, C at these alleles will total to 1.0, so we should just
##  be able to add the alt frequencies to get a compiled frequency for single_base substitutions. However, there is theoretically an unlimited number of possibilities
##  for alt alleles at a position (for the way we track variant positions, where indels are attached to the position of the last reference base position).  ie, A -> TGC, A-> TG, A->CAC, etc.
##  So it would technically be possible for a genomic_position_AAF to be > 1.0.  I don't know if that is an issue or not.   Also, while it makes more intuitive sense
##  to just choose 2 alleles (random numbers) I haven't sat down to figure out if it is mathematically the same probability as making 6 choices and challenging
##  against the smaller individual variant_AAFs.
##
##  GNOMAD SIMULATION CHALLENGE THRESHOLD:
##  I wasn't sure what the best way would be to handle variants observed in one population but not in others.  Specifically,
##  what challenge threshold should be used during simulation for a subpopulation in which the variant being analyzed was not seen.
##  If a variant isn't seen in a subpopulation technically we only can say that, if the allele it circulates, it circulates at
##  a frequency below the minumum_allele_frequency for that cohort [ie, 0.01 (or 1 in 100) if the total number of individuals
##  sampled from the population was 50 ]  So, saying the alt_frequency is 0.0 is not necessarilly accurate. However,
##  using the minimum_frequency resolution also does not seem fair.  Consider an allele that is seen only once in NFE
##  where total allele counts is ~100,000. If AMR theoretical maximum_allele_count is 20,000, then the min_AMR_frequency
##  would be 1 in 20,000.  So the AMR sub_cohort individuals would have a 5X greater chance of simulating the variant than
##  the NFE sub_population, even though that variant wasn't actually seen in the AMR population.
##
##  Almost half of the gnomAD variants are variants seen only once in all populations.  The simulation strategy I decided
##  to use in the end was to say that a variant/alt_allele seen only one time in all subpopulations combined is a singleton
##  that is not actually a circulating allele.  If a variant/alt_allele is observed more than once, even if seen only once in
##  2 different subpopulations, then we consider that variant to be circulating in those populations at those allele frequencies,
##  BUT we also assume that the variant/allele is NOT ciruclating in the populations where the variant was not observed. When simulating
##  the latter populations, 'alleles' generated from the random number generator will be challenged against 0.0; ie, those populations
##  will have no chance of obtaining the variant.
##
##  Also, variants defined as 'singletons' will not be individually challenged in the simulator.  Rather, the pool of singletons
##  will be used to generate a gene specific singleton frequency --> number of singletons observed per total number of alleles sampled.
##  As every individual hs two copies of the gene, the simulator will draw two random number challenges for the possibility of seeing
##  a singleton in the individual, and the numbers will be challenged against the gene specific singleton rate.

##
##
##  Onwards to some code and less rambling...
##################
##################

########################
########################
##
## Simulation Function
##

def altAlleleSimulation(CAUCASIAN_GROUP, simulation_result, GNOMAD_ALTS_DICTIONARY, cutoff_aaf, aaf_cohort, gene_to_analyze, GNOMAD_compiled_variant_AAF_list, pca_AAF_DICTIONARY, minAF_dictionary, SINGLETON_DICTIONARY):

	## individuals are separated into pca groups based on gnomad aaf pca analysis

	pca_aaf_group = aaf_cohort

	#
	# Determine the sample size that needs to be simulated
	#

	cohort_samples = pca_AAF_DICTIONARY[pca_aaf_group]           # dictionary is: key = pca_group ; value = list of individuals categorized as beloning to that group


	if aaf_cohort == 'aaf_gnomad_nfe':                           # This is just a switch I put in to allow you to easily change the whole analysis to use aaf_NWE for the caucasian group, instead of aaf_NFE.  I might just remove this.
		aaf_comparison_category = CAUCASIAN_GROUP
	else:
		aaf_comparison_category = aaf_cohort



	## Remember that the 'simulation_result' passed into the function is: [CURRENT_GENE, cut_off_value, simulation_counter].
	## Each simulation will append simulated_individual_alt_allele_counts to this list; ex) [ NCSTN , 0.01 , 1 , SIM_INDV_#1_ALT_COUNTS, SIM_INDV_#2_ALT_COUNTS, ...]


	#
	# Print pca_ancestry group that is being simulated and the number of individuals in that group
	#

	if simulation_result[2] == 1:

		# This provides a documented record of what the simulator did.  It prints out at the first simulation of a particular AAF Cutoff analysis.  There are some
		# other print functions throughout this code that tracks progress, but this is nice because if you save the print out as a log file you can go back and see:
		#
		# simulating: 30 individuals using aaf_gnomad_afr
		# simulating: 70 individuals using aaf_gnomad_nfe
		# simulating: 17 individuals using aaf_gnomad_amr
		#
		# ... and if everything worked correctly, this should be the same sample pca_group break down as the HS_COHORT, which is also printed out.

		print 'simulating %i %s individuals using %s as threshhold_aaf for cutoff' % (len(cohort_samples), aaf_cohort, aaf_comparison_category)

	else:
		None


	min_AF = minAF_dictionary[aaf_cohort]                                          #  Not currently using, see notes in the 06_gnomad_control_simulations.py script for description


	for individual in range(len(cohort_samples)):

		sim_indv_TOTAL_alts = 0

		for position in GNOMAD_compiled_variant_AAF_list:

			var_pos = position

			gnomad_sub_population_aaf = float(GNOMAD_compiled_variant_AAF_list[position][aaf_comparison_category])

			###
			### Don't need to check whether gnomad is below a cutoff frequency.  We are dealing with compiled genomic position alt frequencies
			### where the frequency is the accumulated frequency of seeing any alt at the position.  Common variants were already filtered
			### during the construction of the positional frequency dictionary.  It is possible that some polymorphic positions
			### will have accumulated frequencies above aaf_cutoff.  For example, if there are 3 variants at a position each with
			### 0.005 population allele frequency.  The accumulated position frequency will be 0.015, above a 1% cutoff; however all of the
			### contributing variants have below 1% aaf, so the simulation should procede using the position_alt_frequency challenge threshold
			### even though it is above the variant_alt_frequency cut off

			if gnomad_sub_population_aaf == -1.0 or gnomad_sub_population_aaf == 0.0:
				pca_comparison_aaf = 0.0
				#pca_comparison_aaf = min_AF

			else:

				pca_comparison_aaf = gnomad_sub_population_aaf

			#
			# Simulate alleles with random number generator to determine if alleles are ALT or REF
			#

			sim_indv_alt_count = 0

			for allele in range(2):

				random_draw = random.random()

				if float(random_draw) < float(pca_comparison_aaf):
					sim_indv_alt_count += 1

				else:
					None

			if sim_indv_alt_count > 0:                                                                 # Track which positions are accumulating ALT alleles in simulation
				if var_pos in GNOMAD_ALTS_DICTIONARY[cutoff_aaf][gene_to_analyze]:
					GNOMAD_ALTS_DICTIONARY[cutoff_aaf][gene_to_analyze][var_pos] += sim_indv_alt_count
				else:
					GNOMAD_ALTS_DICTIONARY[cutoff_aaf][gene_to_analyze][var_pos] = sim_indv_alt_count
			else:
				None



			sim_indv_TOTAL_alts += sim_indv_alt_count


		### Singleton simulation...
		current_gene_SINGLETON_FREQUENCY = SINGLETON_DICTIONARY[gene_to_analyze]['gene_specific_singleton_frequency']

		for allele in range(2):
			random_draw = random.random()

			if float(random_draw) < float(current_gene_SINGLETON_FREQUENCY):
				sim_indv_alt_count += 1

			else:
				None

		simulation_result.append(sim_indv_TOTAL_alts)

	return(simulation_result, GNOMAD_ALTS_DICTIONARY)

########################
########################
##
## Burdens Test Function
##

def burdens_analysis( capture_variant_file , gnomad_variant_file , samples_info_file , GENE_LIST , cut_off_list , NUMBER_OF_SIMULATIONS):

	##################################
	##################################
	##  READ  IN  VARIANT  INFORMATION

	#hs_capture_missense_variants_file = args.CAPTURE_VARIANT_FILE

	hs_capture_vars = extractVariants(capture_variant_file)
	capture_gem_vars_list = [ VARIANT( x , hs_capture_vars[0] ) for x in hs_capture_vars[1:] ]

	#gnomad_missense_variants_file = args.GNOMAD_VARIANT_FILE

	gnomad_vars = extractVariants( gnomad_variant_file )
	gnomad_gem_vars_list = [ VARIANT( x , gnomad_vars[0] ) for x in gnomad_vars[1:] ]

	###########################################
	###########################################
	##  READ  IN  PATIENT  COHORT  INFORMATION
	##
	##  Sample information file is a list of all of the samples to include in the analysis.  It has sample_name, affected_status,
	##  and pca_ancestry_group.  (This file should be only the samples you want to analyze.  ie, affected samples.) The targted
	##  capture variant file has information for ALL SAMPLES CAPTURED, that is to say, both affected and unaffected.  Samples were
	##  captured in 3 batches, and sequenced in 2 batches.  I don't know the distribution of affected samples throughout the batches.
	##  I figured that keeping all samples for all steps prior to the burden analysis would give us more data for site filtering and
	##  for PCA analysis/classification (EDIT: though in the end I removed the unaffected samples for the PCA analysis).
	##
	##  So, when doing ALT allele counts for the HS COHORT we cant just evaluate ALT COUNT at each variant.  We need to loop over
	##  our sample list and add up the ALT COUNTS for the samples in the affected list.
	##

	#########################
	# ALL HS AFFECTED SAMPLES
	sample_info_vars = extractVariants( samples_info_file )

	####################################################
	####################################################
	##
	## ASSESS  COHORT  PCA  MAKEUP
	##
	## Evaluate pca_ancestry break down of our affected sample cohort.  Put this information into reciprocal dictionaries.
	## sample_to_pca_DICTIONARY{}: for any sample_name key, provide pca_group for that sample
	## pca_cohort_dictionary{}: for any pca_group key, provide list of samples beloning to that pca_group
	##
	## This information is going to help us when we need to do pca_ancestry specific gnomAD_simulations. For a simulation,
	## we can loop over the set() of pca_groups.  For every group, simulate len(pca_cohort_dictionary[aaf_gnomad_nfe]) number
	## of individuals.
	##
	###
	###  SEPARATE ANALYSIS SAMPLES INTO PCA COHORT GROUPS

	print 'assessing sample pca distributions...'
	print

	sample_list = []
	pca_list = []
	pca_cohort_dictionary = {}                   # key = pca_group ; value = list of samples that belong to pca_group
	sample_to_pca_DICTIONARY = {}                # key = sample_name ; value = pca_group for that sample

	for x in sample_info_vars[1:]:
		sample = x[0]                            #  These indices are just based on cohort_information_file format.
		pca_group = x[-1]                        #  Could be more elegant by analyzing header names to figure out which index is sample_name and which is PCA_ancestry

		pca_list.append(pca_group)



		if sample in sample_list:
			None

		else:

			sample_list.append(sample)
			sample_to_pca_DICTIONARY[sample] = pca_group

		if pca_group in pca_cohort_dictionary:
			pca_cohort_dictionary[pca_group].append(sample)
		else:
			pca_cohort_dictionary[pca_group] = []
			pca_cohort_dictionary[pca_group].append(sample)

	################################
	################################
	##
	##  Print out some summary information about cohort.
	##

	print 'sample summary:'
	print
	total_samples = 0
	for x in pca_cohort_dictionary.keys():
		print 'number of %s:\t%i' % ( x , len(pca_cohort_dictionary[x]))
		total_samples += len(pca_cohort_dictionary[x])

	aaf_groups_to_analyze = set(pca_cohort_dictionary.keys())

	COHORT_SIZE = len(sample_list)
	print 'analyzing following aaf pca groups:\t',aaf_groups_to_analyze
	print 'total cohort size:\t', COHORT_SIZE

	print
	print

	########################################
	########################################
	##
	##  ORGANIZE VARIANTS
	##
	##  Separate out GNOMAD and CAPTURE variants into VARIANT dictionaries where KEYS are gene names and
	##  values are lists of VARIANTS that affect the KEY gene.  This speeds processing up by not having to
	##  loop over variants in all genes in order to find NCSTN variants.
	##
	##  ALSO, handle 'singleton variants'.  Singletons are separated inton a separate GNOMAD_singleton dictionary. For our analysis
	##  we will assume that a variant seen only once in all population subgroups is a singleton variant.  We will therefore also assume that:
	##
	## 		- a variant/allele seen more than once (which could be once in population_A and once in
	##         population_B) is actually ciruclating in the population(s) at that frequency.
	##       - a variant/allele that is not a 'singleton' and was not observed in a sub-population is
	##         not present/circulating in that sub-population;  thus, when simulating individuals from
	##         that sub-population we will challenge a random number against 0.0, giving no chance for
	##         the individual to observe that variant.

	CAPTURE_VARIANT_DICTIONARY = {}
	for x in GENE_LIST:
		CAPTURE_VARIANT_DICTIONARY[x] = []

	for variant in capture_gem_vars_list:

		genes_associated = variant.info['GENE'].split('/')

		for cur_gene in genes_associated:
			CAPTURE_VARIANT_DICTIONARY[cur_gene].append(variant)

	GNOMAD_VARIANT_DICTIONARY = {}
	for x in GENE_LIST:
		GNOMAD_VARIANT_DICTIONARY[x] = []

	gnomad_singleton_dictionary = {}

	for x in GENE_LIST:
		gnomad_singleton_dictionary[x] = {}
		gnomad_singleton_dictionary[x]['singleton_counts'] = 0
		gnomad_singleton_dictionary[x]['allele_depth_counts'] = []

	for variant in gnomad_gem_vars_list:

		genes_associated = variant.info['GENE'].split('/')

		ac_nfe = float(variant.info['AC_nfe'])
		ac_afr = float(variant.info['AC_afr'])
		ac_amr = float(variant.info['AC_amr'])
		ac_sas = float(variant.info['AC_sas'])
		ac_eas = float(variant.info['AC_eas'])

		an_nfe = float(variant.info['AN_nfe'])
		an_afr = float(variant.info['AN_afr'])
		an_amr = float(variant.info['AN_amr'])
		an_sas = float(variant.info['AN_sas'])
		an_eas = float(variant.info['AN_eas'])

		ac_list = [ ac_nfe , ac_afr, ac_amr, ac_sas , ac_eas ]
		ac_total = sum( ac_list )

		an_list = [ an_nfe , an_afr, an_amr, an_sas , an_eas ]
		an_total = sum( an_list )

		if ac_total == 1:

			for each_gene in genes_associated:
				gnomad_singleton_dictionary[each_gene]['singleton_counts'] += ac_total
				gnomad_singleton_dictionary[each_gene]['allele_depth_counts'].append( an_total )

		else:

			for each_gene in genes_associated:

				GNOMAD_VARIANT_DICTIONARY[each_gene].append(variant)

	for each_gene in GENE_LIST:

		total_singletons = gnomad_singleton_dictionary[each_gene]['singleton_counts']

		gnomad_singleton_dictionary[each_gene]['allele_depth_ARRAY'] = numpy.zeros( int(total_singletons) )

		for x in xrange(0 , int(total_singletons) ):

			gnomad_singleton_dictionary[each_gene]['allele_depth_ARRAY'][x] = gnomad_singleton_dictionary[each_gene]['allele_depth_counts'][x]

		gnomad_singleton_dictionary[each_gene]['mean_allele_depth'] = gnomad_singleton_dictionary[each_gene]['allele_depth_ARRAY'].mean()
		gnomad_singleton_dictionary[each_gene]['mean_allele_stdev'] = gnomad_singleton_dictionary[each_gene]['allele_depth_ARRAY'].std()

		gnomad_singleton_dictionary[each_gene]['gene_specific_singleton_frequency'] = total_singletons/(gnomad_singleton_dictionary[each_gene]['mean_allele_depth'])

	##########################################################
	##########################################################
	##
	##  Set up an alt allele frequency threshold cutoff.  The simulations analysis is set up to perform simulations
	##  for a list of different alt allele frequency thresholds, ie: what does variant burdens look like for variants <0.01 or < 0.001.
	##  However, if a subpopulation had only 500 exomes, alt allele frequencies for variants in this cohort can only be defined down to 0.001.
	##  So there is no point in doing a burdens analysis looking at variants with alt allele frequencies of < 0.0001.
	##
	##	min_aaf_dictionary is used to establish the theoretical maximum allele counts for each subpopulation, and therefore the theoretical
	##  subpopulation minimum alt allele frequencies for variants.  The CUT_OFF_MIN_AAF_RESOLUTION is established, and the code ignores any
	##  requests to perform a burden analysis below the CUT_OFF_MIN_AAF_RESOLUTION.
	##

	# total sample number data collected from https://macarthurlab.org/2018/10/17/gnomad-v2-1/
	min_aaf_dictionary = {}
	min_aaf_dictionary['aaf_gnomad_amr'] = float(1)/float(2*( 17296 + 424 ))  #[exomes + genomes]
	min_aaf_dictionary['aaf_gnomad_afr'] = float(1)/float(2*( 8128 + 4359 ))  #[exomes + genomes]
	min_aaf_dictionary['aaf_gnomad_nfe'] = float(1)/float(2*( 56885 + 7718 ))  #[exomes + genomes]
	min_aaf_dictionary['aaf_gnomad_eas'] = float(1)/float(2*( 9197 + 780 ))   #[exomes + genomes]
	min_aaf_dictionary['aaf_gnomad_sas'] = float(1)/float(2*( 15308 + 0 ))   #[exomes + genomes]
	min_aaf_dictionary['aaf_gnomad_all_subsets'] = float(1)/float(2*( 106814 + 13281 ))   #[exomes + genomes]

	CUT_OFF_MIN_AAF_RESOLUTION = 1.0
	for x in aaf_groups_to_analyze:
		if min_aaf_dictionary[x] < CUT_OFF_MIN_AAF_RESOLUTION:
			CUT_OFF_MIN_AAF_RESOLUTION = min_aaf_dictionary[x]
		else:
			None

	print 'begining variant analysis...'
	print
	print

	#####
	#####
	#####  START  SIMULATION  SCRIPT
	#####
	#####

	#CAUC_GROUP = args.CAUCASIAN_GROUP.lower() # this is a switch I had built in that allowed me to easily change the analysis between using NFE or NWE for the caucasian group
	CAUC_GROUP = 'aaf_gnomad_nfe'

	#
	# Prep for simulations output. Create output lists and setup output headers.
	#

	CAPTURE_RESULTS = []
	capture_results_header = ['gene','cut_off']
	COHORT_TO_ANALYZE = sample_list
	for x in COHORT_TO_ANALYZE:
		capture_results_header.append(x)

	SIMULATION_FINAL_RESULTS = []
	simulation_results_header = ['gene','cut_off','simulation']
	counter = 1
	for x in aaf_groups_to_analyze:
		for y in pca_cohort_dictionary[x]:
			simulation_results_header.append('sim_indv_%i' % (counter))
			counter += 1

	CAPTURE_RESULTS.append(capture_results_header)
	SIMULATION_FINAL_RESULTS.append(simulation_results_header)

	#
	# Set up dictionaries to catch simulation results for evaluating which positions are being frequently chosen as ALT alleles.  Dictionaries of dictionaries of dictionaries.
	# key = SIMULATION_AAF_CUTOFF ; value = { key = gene_name ; value = { key = genomic_position;  value = alt_counts (observed or simulated) }
	#
	# Note that the alt_counts values are the summed values for all simulations.
	# This isn't necessary for any of the analysis.  I just added this in at some point because some genes were outputting large expected alt counts.  I was curious
	# about which variants were being frequently 'hit' as alt alleles in the simulations.  These dictionaries allow you to start to get a feel for that.  So if a genomic
	# position 103250 in NOTCH1 was being chosen as an ALT allele in every other simulation, out of 200 simulations, then:
	#
	# GNOMAD_ALTS_DICTIONARY[ 0.01][ NOTCH_1 ][ 103250 ] = 100
	# You can write a script to output/plot positions linearly and get a feel for whether a particular region is being hit a lot or not.
	#

	GNOMAD_ALTS_DICTIONARY = {}
	CAPTURE_ALTS_DICTIONARY= {}
	for x in cut_off_list:
		GNOMAD_ALTS_DICTIONARY[x] = {}
		CAPTURE_ALTS_DICTIONARY[x] = {}
		for y in GENE_LIST:
			GNOMAD_ALTS_DICTIONARY[x][y] = {}
			CAPTURE_ALTS_DICTIONARY[x][y] = {}


	##############################
	##############################
	##
	## ANALYZE already!
	##
	## ALT counting and simulation analyzes one gene at a time. Then one individual at a time, within that gene.  So you end up tracking
	## (and output-ing) the number of alt alleles for each individual (or simulated individual) within each gene.  This sets you up to be
	## able to easily calculate observed alt alleles for each gene by summing the alt allele counts for all individualas. It also allows
	## you to ask how many individuals are affected by an alt allele for a particular gene within the cohort.
	##

	for CURRENT_GENE in GENE_LIST:                                                # LOOP OVER GENES

		print CURRENT_GENE


		for cut_off_value in cut_off_list:                                          # LOOP OVER CUTOFFS;  WHILE ANALYZING EACH GENE, ANALYZE THE GENE AT ALL CUTOFFS

																					# If the user is asking for an analysis at an AAF threshold below the resolution for the subpopulation with
			if cut_off_value > CUT_OFF_MIN_AAF_RESOLUTION :                         # the fewest samples, do not perform the analysis.

				tmp_result = [CURRENT_GENE , cut_off_value]                         # set up for collecting results for this particular gene/cutoff.  will be a list of [ GENE , CURRENT_CUTOFF, #ALTS_INDV_1 , #ALTS_INDV_2, #ALTS_INDV_3, ...]


				##########################################################
				## COUNT  ALT  ALLELES  IN  TARGETED_CAPTURE  COHORT
				##

				for indv in COHORT_TO_ANALYZE:                                      # LOOP OVER INDIVIDUALS; DETERMINE NUMBER OF ALT ALLELES FOR EACH INDIVIDUAL WITHIN THIS GENE AT THIS CUTOFF

					patient_TOTAL_alts = 0

					for cap_variant in CAPTURE_VARIANT_DICTIONARY[CURRENT_GENE]:    # LOOP OVER VARIANTS

						variant_min_freq = 1.0                                      # Identify gnomAD_subpopulation_aaf to use for AAF threshold filtering.
						THRESHHOLD_AAF_TO_USE = 'None'                              # Find the smallest subpopulation AAF that is NON-ZERO and NOT -1
						temp_freq_list = []

						for aaf_group in aaf_groups_to_analyze:

							temp_freq_list.append( float( cap_variant.info[aaf_group] ) )

							if float(cap_variant.info[aaf_group]) < variant_min_freq and float(cap_variant.info[aaf_group])!= 0.0 and float(cap_variant.info[aaf_group]) != -1.0:

								variant_min_freq = float(cap_variant.info[aaf_group])
								THRESHHOLD_AAF_TO_USE = aaf_group

							else:
								None

						if all([x <=0.0 for x in temp_freq_list]):
							variant_min_freq = 0.0
						else:
							None


						patient_alt_count = 0                                            #  Patient alt count for this one variant being analyzed

						if variant_min_freq < cut_off_value:

							patient_alt_count = cap_variant.info[indv].count('1')

							var_chrom = cap_variant.info['CHROM']
							var_pos = cap_variant.info['POS']
							cap_gene = cap_variant.info['GENE']
							cap_alt = cap_variant.info['ALT']
							cap_ref = cap_variant.info['REF']
							variant_ID = '%s_%s_%s_%s' % (var_chrom,var_pos,cap_ref,cap_alt)

						else:
							None

						patient_TOTAL_alts += patient_alt_count


						if patient_alt_count > 0:                                                                       # Track at which genomic positions alt alleles are being found
							if var_pos in CAPTURE_ALTS_DICTIONARY[cut_off_value][CURRENT_GENE]:
								CAPTURE_ALTS_DICTIONARY[cut_off_value][CURRENT_GENE][var_pos] += patient_alt_count
							else:
								CAPTURE_ALTS_DICTIONARY[cut_off_value][CURRENT_GENE][var_pos] = patient_alt_count
						else:
							None

					tmp_result.append(patient_TOTAL_alts)

				CAPTURE_RESULTS.append(tmp_result)

				#############################
				## RUN  GNOMAD  SIMULATIONS
				##

				#
				# Create multiallelic genome position ALT frequencies data from gnomAD variants
				#

				GNOMAD_COMPILED_VARIANT_POSITIONS = {}                                          # key = genomic position ; value = { key = gnomad_subpopulation ; value = combined multiallelic ALT allele frequency}

				for variant in GNOMAD_VARIANT_DICTIONARY[CURRENT_GENE]:

					gnom_variant_min_freq = 1.0
																								# Filter variants below AAF cutoff threshold
					for aaf_group in aaf_groups_to_analyze:
						if float(variant.info[aaf_group]) < gnom_variant_min_freq and float(variant.info[aaf_group])!= 0.0 and float(variant.info[aaf_group]) != -1.0:
							gnom_variant_min_freq = float(variant.info[aaf_group])
						else:
							None

					variant_position = variant.info['POS']

					if gnom_variant_min_freq < cut_off_value:

						if variant_position in GNOMAD_COMPILED_VARIANT_POSITIONS:              # Compile gnomAD variant list to merge variant aaf for variants at same position

							for sub_group in pca_cohort_dictionary.keys():

								if sub_group == 'aaf_gnomad_nfe':
									sub_group = CAUC_GROUP

								else:
									None

								if float(variant.info[sub_group]) == -1:
									GNOMAD_COMPILED_VARIANT_POSITIONS[variant_position][sub_group] += 0.0
								else:
									GNOMAD_COMPILED_VARIANT_POSITIONS[variant_position][sub_group] += float(variant.info[sub_group])

						else:

							GNOMAD_COMPILED_VARIANT_POSITIONS[variant_position] = {}

							for sub_group in pca_cohort_dictionary.keys():

								if sub_group == 'aaf_gnomad_nfe':
									sub_group = CAUC_GROUP

								else:
									None

								GNOMAD_COMPILED_VARIANT_POSITIONS[variant_position][sub_group] = 0

								if float(variant.info[sub_group]) == -1:
									GNOMAD_COMPILED_VARIANT_POSITIONS[variant_position][sub_group] += 0.0
								else:
									GNOMAD_COMPILED_VARIANT_POSITIONS[variant_position][sub_group] += float(variant.info[sub_group])

					else:
						None

				#
				# Run simulations
				#

				simulation_counter = 0

				for simulation in range(NUMBER_OF_SIMULATIONS):

					simulation_counter += 1

					simulation_result = [CURRENT_GENE, cut_off_value, simulation_counter]

					if simulation_result[2] == 1:
						print 'Running %i simulations for aaf cutoff: %f' % (NUMBER_OF_SIMULATIONS,cut_off_value)

					else:
						None

					for aaf_cohort in pca_cohort_dictionary:

						simulation_result , GNOMAD_ALTS_DICTIONARY = altAlleleSimulation(CAUC_GROUP , simulation_result, GNOMAD_ALTS_DICTIONARY, cut_off_value, aaf_cohort , CURRENT_GENE, GNOMAD_COMPILED_VARIANT_POSITIONS, pca_cohort_dictionary, min_aaf_dictionary , gnomad_singleton_dictionary)


					SIMULATION_FINAL_RESULTS.append(simulation_result)

					if simulation_counter % 20 == 0:
						print 'finished simulation number %i' % (simulation_counter)

					else:
						None
			else:
				None

			sys.stdout.flush()

	return ( SIMULATION_FINAL_RESULTS , CAPTURE_RESULTS )

###########################################################################################################
###########################################################################################################
###########################################################################################################
###
###   OK, NOW ACTUALLY IMPLEMENT ALL THAT BORING STUFF CODED UP THERE.
###
####  RUN VARIANT ENRICHMENT ANALYSIS
###
###

VEP_variants_dir = '../../output/VEP_variant_files'
samples_info_dir = '../../output/samples_info_output'
burdens_output_dir = '../../output/burdens_results'
transcript_dir = '../../output/transcript_file'

hs_capture_MISSENSE_variants_file = '%s/%s' % ( VEP_variants_dir , 'HS_CAPTURE_missense_variants_capturedHSpatients_only.txt' )
gnomad_MISSENSE_variants_file = '%s/%s' % ( VEP_variants_dir , 'gnomad_MISSENSE_variants.txt' )

hs_capture_SYNONYMOUS_variants_file = '%s/%s' % ( VEP_variants_dir , 'HS_CAPTURE_synonymous_variants_capturedHSpatients_only.txt' )
gnomad_SYNONYMOUS_variants_file = '%s/%s' % ( VEP_variants_dir , 'gnomad_SYNONYMOUS_variants.txt' )

samples_info_file = '%s/%s' % ( samples_info_dir , 'all_cohorts_HS_AFFECTED_captured_samples_withPCAinfo.txt' )

transcript_file = '%s/%s' % ( transcript_dir , 'HS_capture_primary_transcript_exon_sequences_WITH_DOMAINS.txt' )
transcript_vars = extractVariants(transcript_file)
#
# Collect unique list of genes to analyze from sequence file.
#

GENE_LIST = []
for x in transcript_vars[1:]:
	temp_exon = Exon(x)
	if temp_exon.gene_name in GENE_LIST:
		None
	else:
		GENE_LIST.append(temp_exon.gene_name)

cut_off_list = [ 0.05 , 0.01 , 0.005 , 0.001 ]
NUMBER_OF_SIMULATIONS = 200

MISSENSE_simulation_results , MISSENSE_capture_results = burdens_analysis( hs_capture_MISSENSE_variants_file , gnomad_MISSENSE_variants_file , samples_info_file , GENE_LIST , cut_off_list , NUMBER_OF_SIMULATIONS)

SYNONYMOUS_simulation_results , SYNONYMOUS_capture_results = burdens_analysis( hs_capture_SYNONYMOUS_variants_file , gnomad_SYNONYMOUS_variants_file , samples_info_file , GENE_LIST , cut_off_list , NUMBER_OF_SIMULATIONS)

outPut_NODATE( MISSENSE_simulation_results , '%s/%s' % ( burdens_output_dir , 'MISSENSE_simulation_burdens_results.txt' ) )
outPut_NODATE( MISSENSE_capture_results , '%s/%s' % ( burdens_output_dir , 'MISSENSE_capture_burdens_results.txt' )	)

outPut_NODATE( SYNONYMOUS_simulation_results , '%s/%s' % ( burdens_output_dir , 'SYNONYMOUS_simulation_burdens_results.txt' ) )
outPut_NODATE( SYNONYMOUS_capture_results , '%s/%s' % ( burdens_output_dir , 'SYNONYMOUS_capture_burdens_results.txt' ) )

##################################################################################################################

'''
##
##  This is some old code that I had used to print out the data within GNOMAD_ALTS_OUT_LIST to get an idea of
##  which positions were frequent contributors to the gnomAD simulation ALT counts.  I have just kept in here just
##  in case I want to go back to it.  It is not necessary for the analysis.
##
##  I think that I was able to use this clean up noise in the MAML gene analyses.  I think this is what led to me
##  doing the reciprocal quality control sites filtering.  I think there were a lot of variants like indels in the
##  MAML polyQ tracts that were low depth, or had a high percentage of '.' genotype calls in the capture.
##
##

GNOMAD_ALTS_OUT_LIST = []
tmp_header = ['gene', 'aaf_cutoff', 'position','counts','avg','pct_of_vars']
GNOMAD_ALTS_OUT_LIST.append(tmp_header)

for x in GNOMAD_ALTS_DICTIONARY.keys():
	unsorted_gene_list = []
	aaf = x
	for y in GNOMAD_ALTS_DICTIONARY[x].keys():
		gene = y
		for z in GNOMAD_ALTS_DICTIONARY[x][y].keys():
			position = z
			counts = GNOMAD_ALTS_DICTIONARY[x][y][z]
			nfe_aaf = GNOMAD_COMPILED_VARIANT_POSITIONS[position]['aaf_gnomad_nfe']
			print nfe_aaf
			tmp = [gene ,aaf ,position, counts,nfe_aaf]
			unsorted_gene_list.append(tmp)

	total_variants = 0
	for each in unsorted_gene_list:
		avg_var_observed = float(each[3])/float(NUMBER_OF_SIMULATIONS)
		each.append('%.4f' % (avg_var_observed))
		total_variants += avg_var_observed

	for each in unsorted_gene_list:
		pct_contr = (float(each[5])/float(total_variants))*100
		each.append('%.4f' % (pct_contr))

	unsorted_gene_list.sort(key = lambda z: (-float(z[-1]),float(z[2])))
	for each in unsorted_gene_list:
		GNOMAD_ALTS_OUT_LIST.append(each)

outPut_NODATE(GNOMAD_ALTS_OUT_LIST,'gnomad_alt_positions_list.txt')

singleton_contribution = float(0)
sim_totals = []
tracker = 3
for x in range(2):
	tmp_total = 0
	for y in GNOMAD_ALTS_OUT_LIST[1:]:
		tmp_total += float(y[tracker])
	sim_totals.append(tmp_total)
	tracker += 2

for x in GNOMAD_ALTS_OUT_LIST[1:]:
	if float(x[4]) < 0.00001:
		singleton_contribution += float(x[5])
	else:
		None

print
print sim_totals
print
print 'avg number of singletons:\t', singleton_contribution

CAPTURE_ALTS_OUT_LIST = []
tmp_header = ['gene', 'aaf_cutoff', 'position','counts']
CAPTURE_ALTS_OUT_LIST.append(tmp_header)

for x in CAPTURE_ALTS_DICTIONARY.keys():
	capture_unsorted_list = []
	aaf = x
	for y in CAPTURE_ALTS_DICTIONARY[x].keys():
		gene = y
		for z in CAPTURE_ALTS_DICTIONARY[x][y].keys():
			position = z
			counts = CAPTURE_ALTS_DICTIONARY[x][y][z]
			tmp = [gene,aaf,position,counts]
			capture_unsorted_list.append(tmp)

	capture_unsorted_list.sort(key = lambda z: (z[0],-float(z[1]),float(z[2])))
	for each in capture_unsorted_list:
		CAPTURE_ALTS_OUT_LIST.append(each)

outPut_NODATE(CAPTURE_ALTS_OUT_LIST,'capture_alt_positions_list.txt')
'''
