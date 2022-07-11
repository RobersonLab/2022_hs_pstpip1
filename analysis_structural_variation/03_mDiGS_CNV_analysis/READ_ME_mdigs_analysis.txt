The code in this folder performs summary analyses of the HS capture CNV analysis and evaluates gnomAD structural
variants for the purposes of generating a summary variants table and supplementary CNV coverage plots for the paper.

1 - 03.1_structural_sites_extraction.sh :  Uses bcftools to extract structural variants from the gnomAD SV database
    that fall within our targeted capture region. 

2 - 03.2_gnomAD_polymorphic_sv_filter.py :  Identify polymorphic CNV 

3 - 03.3_exon_overlapping_HS_CNV.py : Identify any HS-cohort CNV that overlap exonic regions in genes of our targeted capture panel.

4 - 03.4_CNV_visualizations.RMd : Plot out examples of coverage plots for the called CNV for supplementary figures.
	
5 - 03.5_gnomAD_sv_table_variants.py : This script is very similar to 03.2.  It prints out variants from the gnomAD SV database,
    but it prints out variants for the papers final HS-capture/gnomAD CNV comparison table.  There were two CNV identified in our
	HS cohort that were not polymorphisms.  One is not in the gnomAD variant databse as far as I can tell. The other looks like
	it is in the gnomAD SV database, but with a slightly different CNV positional window.  This script prints out both the
	gnomAD polymorphisms and the non-polymorphic capture CNV.
    

	
More detail descriptions of normalization strategies and coding documentation are in each script.
 