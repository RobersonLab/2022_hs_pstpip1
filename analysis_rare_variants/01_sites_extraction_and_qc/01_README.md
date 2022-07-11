## Note 1
The initial sites extraction script throws the following error messsage for each gnomAD genome sites files, HOWEVER, it extracts all
    of the data that it is supposed to:

    [E::hts_hopen] Failed to open file ../../data/gnomad_sites_files/genome_sites/gnomad.genomes.r2.1.1.sites.11.vcf.bgz.tbi
    [E::hts_open_format] Failed to open file ../../data/gnomad_sites_files/genome_sites/gnomad.genomes.r2.1.1.sites.11.vcf.bgz.tbi
    Failed to open ../../data/gnomad_sites_files/genome_sites/gnomad.genomes.r2.1.1.sites.11.vcf.bgz.tbi: Exec format error

## Note 2
For the sites extraction script, the order in which the gnomAD population subgroups are extracted was intentional. The genome sites
    do not have SAS (south asian) information. I extracted SAS as the last subgroup for the exomes so that I could append SAS information of 
	0 alleles to the genomes data later, and then SAS would be last for both exomes and genomes and make it simpler for merging the two sets
	of data.
	
	So this is just to say that if you were to write a different gnomAD sites extraction script for a different targeted capture project,
	you should order SAS as the last group if you want everything to work smoothly with the rest of this pipeline.


01_sites_extraction_and_QC

The code in this folder collects and QC filters genetic variants from the targeted capture sequencing file and from the gnomad genomes and exomes
sites files. 

1 - 01.1_variant_sites_extraction.sh :  Uses bcftools to query and extract variant information from gnomad and targeted capture VCF.

2 - 01.2_gnomad_coverages_filter.py : This script filters the gnomAD coverages file to keep only sites within the targeted capture
    region.  The coverages file will be used to perform quality control filtering on gnomAD sites.  But is a very large file containing
	gnomAD coverage information for the entire genome.  The subsequent QC steps will be easier to run with a file containing only
	the coverages for the targeted capture region.

3 - 01.3_gnomad_genomic SAS_append.py : The gnomad genomes sites file does not contain information for South Asian (SAS) ancestry 
    due to the small sample size of SAS genomes in the gnomad collection.  In order to combine gnomAD genome and exome data in 
	subsequent steps it will be easier to have a place holder field for SAS in the genome sites file.  This script adds AC = 0 and 
	AN = 0 for SAS to the extracted gnomAD genome sites files. 
	
4 - 01.4_sites_quality_control.py : This script filters gnomad and targeted capture variants that do not pass quality control. 
    The quality control strategy requires variants to meet particular thresholds for depth, uncalled alleles, etc. The QC 
	strategy is described in more detail within the script.

5 - 01.5_sites_merge.py : This script compiles the gnomad exome and gnomad genome data.  AC and AN are pooled to generate 
    one alt allele frequency value for gnomad variants. 
