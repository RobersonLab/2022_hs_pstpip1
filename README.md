# 2022 hidradenitis suppurativa *PSTPIP1* genetics project

Preprint: [medRxiv](https://www.medrxiv.org/content/10.1101/2022.07.12.22277558v1)

Paper: [PubMed](https://pubmed.ncbi.nlm.nih.gov/37013170/)

Supporting data: [Western Blots](https://figshare.com/projects/2022_hidradenitis_suppurativa_PSTPIP1_missense_enrichment/139684)

## Project description
This project involved identifying genes with an enrichment of rare, missense variants in hidradenitis suppurativa. The repository contains the code used to analyze the genetic data that was generated for the work. The samples were sourced from retrospective biological archives, which were not broadly consented to share the raw FASTQ files.

The expected burden was calculated by simulating similar cohorts using gnomAD population-specific allele frequencies. While the raw FASTQ files aren't shareable, the simulations were all based on the number of detected rare variants given the PCA ancestry assignment of the samples. That information alone would be sufficient to reproduce the findings.

## Code description
There were two main arms of the project: 1) analysis of rare missense variants and 2) structural variation analysis. These two arms are in different top level folders with python or R code that would need to be run in order. The markdown folder contains R code primarily used to generate figures for the paper. **Beware** workflow management wasn't used for this project. However, I did independently clean up the code so that the entire analysis correctly ran from a separate user account in different directories.

### Annotation information
The annotations are based on the GRCh37 human genome build.

## Software
* Python
* R

### Python packages
* argparse
* gzip
* numpy
* random
* re
* setuptools
* sys

### R packages
* class
* copynumber
* cowplot
* dplyr
* ggplot2
* ggpubr
* ggrepel
* ggrepel
* gplots
* grid
* gridExtra
* here
* magrittr
* patchwork
* pheatmap
* plyr
* RColorBrewer
* readr
* reshape2
* sjPlot
* stringr
* tidyr
* tidyverse
* umap
* UpSetR

### Setting up custom python classes
Lots of the Python work requires custom classes and functions. They are included in the 'custom_functions_y_classes' directory. To prepare these for use, do the following.

```bash
cd custom_functions_y_classes
conda create -n psptip1 python=2.7.18
conda activate pstpip1
pip install .
```

