The code in this folder performs copy number segmentation analysis on mDiGS normalized coverages from 01_mDiGS_coverage_normalization.
Identified segments are then put through a quality control pipeline to call what we believe are high quality CNV calls.

1 - 02.1_segmentation_calling_and_normalization.RMd : This script uses the R 'copynumber' package to call segmentations based on
    the mDiGS normalized coverages to look for copy number variation. Additionally, the segmentation coverages are median normalized
	and several other data manipulations are performed in order to perpare the segmentation data for 02.2.  This analysis is 
	performed on each capture-batch of samples separately

2 - 02.2_CNV_filtering.RMd :  Median-normalized segmentation data are put through a quality control filtering pipeline in order 
    to identify high quality CNV calls from 02.1.
	
3 - 02.3_pstpip1_copy_population_correction.py : Because the PSTPIP1 deletion CNV polymorphism was so prevalent in our cohort (>50%),
    the mean coverage normalization process is dramatically skewed as the average coverage is closer to 1-copy than 2-copies. This 
	skewing is so extreme that median-normalization of the CNV calls does not correct for the skewing, and in Batches 1 and 3 the only
	copy number populations called in the cohort were 0,2, and 4.  This CNV is called in the gnomAD structural database as a deletion,
	and not as an amplification.  Furthermore, it is highly unlikely that our cohort would have large populations of homozygous deletions
	and homozygous amplifications without observing one individual who was heterozygous for at deletion or insertion.  Rather what makes
	more sense is that as a result of the high prevalence of the CNV, the copy-number data has been normalized to the 1-copy population,
	making 1-copy individuals be reported as 2-copy individuals... and likewise, actual 2-copy individuals (who have 2X coverage as the
	1-copy individuals) are reporting as 4-copy because the 1-copies are reporting as 2-copies.  Thus I believe that the copy number calls
	need to be corrected from 0,2,4 to 0,1,2 respectively. 
	
	This python script is essentially equivalent to a manual correction of these values. It isn't the most elegant solution. But I strongly
	believe that the 0,2,4 copy population calls is a mis-interpretation of the data and needs to be corrected.

4 - 02.4_CNV_compilation.RMd :  The CNV results from 02.2 for each capture-batch are pooled together to produce one table
    of all CNV calls for the entire HS capture cohort.
	
5 - 02.5_CNV_COMPILATION_manual_populations_correction.RMd :  Almost the same as 02.4, but this produces the compiled table using the copy
    number population corrected files for Batch1 and Batch3 samples that are produced by 02.3

6 - 02.6_CNV_population_manual_correction_comparison.RMd :  This script produces a table of just the copy number population distributions 
    for the three sample batches but it produces for both the non-corrected and the corrected sample batches so that you can see how the 
	correction changes the data.

	
	
	
More detail descriptions of normalization strategies and coding documentation are in each script.
