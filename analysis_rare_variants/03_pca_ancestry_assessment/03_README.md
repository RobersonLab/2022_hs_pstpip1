The goal of the code in 03_pca_ancestry_assessment is to unbiasedly categorize the patients in our HS cohort into
   gnomAD ancestry subpopulations (AFR, AMR, EAS, NFE, or SAS). gnomad variant data is used to simulate genotypes
   of these gnomAD ancestries across the targeted capture region.  A PCA analysis is used to evaluate variance between
   these ancestries and our HS patient genotypes.  The resulting principal component values for the simulated gnomad
   genotypes are used as input to train a kNN analysis, which then estimates the most likely ancestry for each patient.
   
   The PCA ancestry analysis has an R markdown that has to be run as a step of the pipeline. Python code is used to 
   generate PCA genotype tables for HS patients, and to generate PCA genotypes for subpopulations simulated
   from gnomAD data. Then a PCA analysis and kNN analysis has to be run in the 03.4_HS_capture_PCA_assessment.RMd.
   This generates a file with the assigned estimated gnomAD ancestry for each HS patient, and then 03.5_PCA_append_assigned_aaf.py
   can be run.
   

1 - 03.1_PCA_ANALYSIS_capture_exome_QC.py : This script generates a new set of variants to use for gnomad genotype simulation.
    As in the 01_sites_extraction section, this script performs quality control filteringg on gnomAD exome and HS targeted capture 
	sites to generate a refined pool of HS capture variants. The genome sites are not used for the PCA analyis because they do not 
	contain SAS genomic information.
	
2 - 03.2_PCA_ANALYSIS_exome_aaf_append.py : For each variant in the QC HS capture variant list.  The corresponding gnomad exome 
    alt allele frequency is appended to the variant information.  This will allow us to simulate this HS capture variant set in 
	subsequent steps. 
	
3 - 03.3_PCA_table_generator.py : Converts HS targeted capture individual genotypes into table for PCA analysis where genotype
    is converted from ALLELE/ALLELE to 0,1, or 2 to indicate the number of alternate alleles the individual has at a particular
	genomic site. In addition simulations are performed to generate PCA table converted genotypes for simulated gnomad individuals
	of various genetic ancestry backgrounds.

4 - 03.4_HS_capture_PCA_assessment.RMd : R markdown file that performs PCA analysis on the data generated from 03.3_PCA_table_generator.py.
    Then the PCA data from a cohort of simulated gnomad individuals is used to train a kNN algorithm.  Data from a second cohort of 
	simulated gnomad individuals is used to identify the optimal K value for kNN.  Then kNN assessment is performed on the targeted
	capture cohort.
	
5 - 03.5_PCA_append_assigned_aaf.py : The assigned PCA ancestry from the kNN analysis is appended to the patient information file. 
