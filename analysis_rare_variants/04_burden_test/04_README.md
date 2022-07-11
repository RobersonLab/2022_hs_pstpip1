The code in 04_burden_test finally performs the variant burden enrichment analysis.

1 - 04.1_gnomad_cohort_simulations.py : This script does two things.  First, it tallies up the number of alt alleles for 
    each gene of interest observed in the HS targeted capture cohort. Second, it performs simulations over the gnomad variant list
    to identify how many alt alleles are oberved for each gene in a cohort of gnomad simulated individuals. The results are returned
    as tables where each column is an individual, and each row is a gene. The value at the intersection is the number of alt alleles
	observed for that individual in that gene (either in the targeted capture or for a particular gnomad simulation).
	
2 - 04.2_burdens_analysis.Rmd : This R markdown analyzes the simulations performed in step #1. It calculates the total number of alt
    alleles for each gene in the entire HS or simulated cohort.  It calculates the mean alt allele count observed for each gene over
	all gnomad simulations and performs statistical testing to determine if the number of variants observed in our targeted capture
	is statistically significant. It also produces figures of the burdens testing results for the manuscript.
