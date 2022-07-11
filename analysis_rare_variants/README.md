The bin directory contains the code for performing the variant enrichment analysis. The code should be executed
in the following order:

1 - PRE_analyses
2 - 01_sites_extraction_and_qc
3 - 02_VEP_annotation
4 - 03_pca_ancestry_assessment
5 - 04_burden_test
6 - supplementary_analyses
7 - code in 'markdowns' directory will further process data generated from bin code (particularly the supplementary_analyses)
    and will generate figures. for example the markdown files in the markdowns directory generate the figures describing
	the targeted capture performance.
	
More details about the executing the code at each step is located in the READ_ME.txt files within each bin subdirectory.