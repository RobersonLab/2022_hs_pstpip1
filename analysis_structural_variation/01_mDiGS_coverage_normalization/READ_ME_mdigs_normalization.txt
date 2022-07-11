The code in this folder generates capture-batch normalized copy number estimates for sliding window genomic regions
over the captured regions for each sample.

1 - 01.1_capture_create_coverage_files.sh:  This script runs samtools depth to evaluate coverage over the targeted
    capture regions for each individual in the HS cohort. 

2 - 01.2_capture_coverages_compile.py:  This script compiles all of the individual coverage files generated 
    by 01.1 into one large cohort coverages file.

3 - 01.3_mDiGS_coverage_normalization.py:  This script separates samples into capture-cohort batches and normalizes the
    coverages for all samples within capture batches to the mean coverage of the capture batch.  Coverages are poole into
	sliding windows of a size determined by the user (50bp for DJMH). *** NOTE *** At the end, the normalized relative
	coverages at each position are multiplied by 2 so that the final values returned are actually estimated copy-number
	based on assuming the average copy-number of the cohort is 2.
	
	
More detail descriptions of normalization strategies and coding documentation are in each script.
