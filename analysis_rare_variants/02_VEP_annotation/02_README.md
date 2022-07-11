
Annotation of genetic variants is performed using the web interface for Ensembl-VEP. One should be able to install VEP for use on the 
command line if desired.  Currently, the code is set up to:

1 - 02.1_variant_filtering.py : Filter variant sites files for the HS capture keep only variants that were observed in HS affected
    individuals included in the burdens analysis.  For our specific targeted capture, in addition to capturing HS affected individuals
	we also captured some samples that were from affected or unaffected relatives of the HS affected proband, or non-affected control samples
	for another HS project. In this script I filter variants that were observed only in samples that would not be included in the HS
	variant burdens analysis.
	
2 - 02.2_ensembl_format_conversion.py : The variant information that I extracted from the gnomad sites files or from our HS capture VCF
    needed to be converted to a format compatible for VEP annotation.  This script generates that VEP compatible file. 

    This generates two files to submit for VEP annotation: one for the gnomad variants, and one for the HS capture variants.  Both of
    these files are names '*_submission.txt' 	
	
3 - WEB BASED VEP ANNOTATION:  

	*** NOTE ***  MAKE  SURE  YOU  ARE  USING  THE GRCh37.p13 VERSION OF VEP !!!

	Annotate variants using VEP on Ensembl website. Use output files from 02.VEP.2_ensembl_format_conversion.py to
    run VEP annotation.  The command line equivalent of the options set for the web based annotation is: 
    
	./vep --appris --biotype --buffer_size 500 --check_existing --distance 0 --hgvs --mane --polyphen b --regulatory --sift b --species homo_sapiens --symbol --transcript_version --tsl --cache --input_file [input_data] --output_file [output_file] --port 3337

	
	(On the results page of the website, you can look under "Job Details" to see the command line equivalent. So if you re-run the analysis
	and want to compare back to my annotation settings, you can compare your job details to my original job above.)
	
	I don't think most of the options really matter.  The few important options are:
		--distance 0 : I don't really know that this will matter, but if the upstream/downstream distance is set > 0, then more variants will be 
		               annotated as being associated with particular genes.  These should all be non-coding, so it shouldn't affect the analysis
					   of rare-missense variant burden. But it might affect the summary statistics for how many variants are associated with each gene.
					   If one were interested in doing an analysis on (potential) regulatory affecting variant burdens, then this ought to be set to
					   something other than 0.
		
		--transcript_version : in order to do a burdens analysis on a particular isoform, the annotations need to be performed so that a variant is 
		               interpreted for all transcript versions.  Then in the next step, variant annotations can be filtered so that only annotations
					   for the transcripts of interest are kept.
					   
		--hgvs :       If this isn't set, then do set it so that you can report variants in tables using hgvs nomenclature.
		
		OTHER:  I think a lot of the other information is lagniappe.  I turned off having VEP report allele frequency information. For filtering variants
		        into noncoding/missense/synonymous you will need the annotation to report cDNA position, amino acid changes, etc, but all of that was default.
							
4 - Return VEP annotation file to the 'VEP_variant_files' directory as: "gnomad_vep_annotations.txt" and "hs_capture_vep_annotations.txt"			

5 - 02.3_vep_annotation_append.py :  This code takes in the filtered variant list from step #1 and the VEP annotation file from step #3. It appends VEP 
    annotation information to each variant ( gene that the variant is associated with, protein coding effect ) and it filters the variants into Non-coding, 
	Synonymous, and Missense variant lists/files. For the purposes of this analysis, missense variants is a slight misnomer.  In reality the 'missense' variants
	are variants that affect protein sequence of the proteins of interest... so, true missense variants, small indels, frameshifts, deletions spanning intron/exon
	junctions. 

VEP annotation for hs_capture variants:

./vep --appris --biotype --buffer_size 500 --check_existing --distance 0 --hgvs --mane --polyphen b --regulatory --sift b --species homo_sapiens --symbol --transcript_version --tsl --cache --input_file [input_data] --output_file [output_file] --port 3337

VEP annotations for gnomad variants:

./vep --appris --biotype --buffer_size 500 --check_existing --distance 0 --hgvs --mane --polyphen b --regulatory --sift b --species homo_sapiens --symbol --transcript_version --tsl --cache --input_file [input_data] --output_file [output_file] --port 3337
