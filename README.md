-   [Installation](#installation)
    -   [System requirements](#system-requirements)
    -   [R packages](#r-packages)
    -   [Datasets](#datasets)
        -   [COSMIC](#cosmic)
        -   [CGC](#cgc)
        -   [Tumor Suppressors](#tumor-suppressors)
        -   [CDS](#cds)
-   [Comprehensive search for iSTOP targetable sites](#comprehensive-search-for-istop-targetable-sites)
-   [Off-target estimates](#off-target-estimates)
-   [Analysis of COSMIC nonsense mutations](#analysis-of-cosmic-nonsense-mutations)
-   [Reproducing Figures](#reproducing-figures)
-   [Session Information](#session-information)

This project accompanies the [iSTOP](github.com/CicciaLab/iSTOP) R package, and is intended to aid reproduction of figures for a forthcoming publication detailing the scope of possible nonsense mutations capable using [CRISPR mediated base editing](https://www.nature.com/nature/journal/v533/n7603/abs/nature17946.html).

Installation
============

To reproduce the analysis in [TBD](link/to/publication), begin by installing [R](https://cran.r-project.org) (~100 MB) and [RStudio](https://www.rstudio.com/products/rstudio/download/) (~500 MB).

Then, go to [github.com/CicciaLab/iSTOP-paper](github.com/CicciaLab/iSTOP-paper) and clone the project to your computer using the green `Clone or download>Download ZIP` button near the top of the page. Then unzip and open the project by clicking on the `iSTOP-paper.Rproj` file. This will open RStudio with the working directory set to the project's folder.

System requirements
-------------------

This analysis requires ~11 GB of disk space, at least 4 GB of memory, and some patience. While the analysis should be cross-platform, it has only been tested on macOS Sierra Version 10.12.5 with a 4 GHz processor and 32 GB of memory. Processing time on this system using 4 cores for just the Human genome is ~30 minutes to locate all iSTOP codons and targets, and an additional ~40 minutes for each RFLP annotation width (i.e. `add_RFLP(width = 150)` takes ~40 minutes).

R packages
----------

To install all necessary R packages, run the following R commands in the RStudio console.

``` r
# Source the Biocoductor installation tool - installs and loads the 
# BiocInstaller package which provides the biocLite function.
source("https://bioconductor.org/biocLite.R")

# Install the following required packages ~ 350 MB 
BiocInstaller::biocLite(c(
  # Packages from CRAN (cran.r-project.org)
  'tidyverse',
  'assertthat',
  'pbapply',
  'devtools',
  'fuzzyjoin',
  # Packages from Bioconductor (bioconductor.org)
  'BSgenome',
  'Biostrings',
  'GenomicRanges',
  'IRanges'
))

# Install genomes and annotation packages ~ 3 GB
BiocInstaller::biocLite(c(
  # Arabidopsis (not available from UCSC)
  'TxDb.Athaliana.BioMart.plantsmart28', # ~ 24  MB
  'org.At.tair.db',                      # ~ 239 MB
  # Genomes
  'BSgenome.Hsapiens.UCSC.hg38',         # ~ 802 MB
  'BSgenome.Celegans.UCSC.ce11',         # ~ 25  MB
  'BSgenome.Dmelanogaster.UCSC.dm6',     # ~ 36  MB
  'BSgenome.Drerio.UCSC.danRer10',       # ~ 343 MB
  'BSgenome.Mmusculus.UCSC.mm10',        # ~ 683 MB
  'BSgenome.Rnorvegicus.UCSC.rn6',       # ~ 719 MB
  'BSgenome.Scerevisiae.UCSC.sacCer3',   # ~ 3   MB
  'BSgenome.Athaliana.TAIR.TAIR9'        # ~ 34  MB
))

# If all goes well, install the iSTOP package hosted on GitHub
BiocInstaller::biocLite('CicciaLab/iSTOP')
```

Datasets
--------

### COSMIC

Download the "COSMIC Mutation Data" (~300 MB compressed) from the [Catalogue of Somatic Mutations in Cancer](https://cancer.sanger.ac.uk/cosmic/download). This dataset requires [registration with a valid email address](https://cancer.sanger.ac.uk/cosmic/register), then in your terminal (not the R console!) you can download the file with the following commands (substitute `your_email_address` with the email used to register).

    sftp "your_email_address"@sftp-cancer.sanger.ac.uk
    # You will be prompted for the password you provided when you registered
    get /files/grch38/cosmic/v80/CosmicMutantExport.tsv.gz

Move this file to the `data/COSMIC` directory of the iSTOP-paper project. It should already be named `CosmicMutantExport.tsv.gz`.

### CGC

Download the [Cancer Gene Census (CGC) dataset](http://cancer.sanger.ac.uk/census) by clicking on the `CSV` Export button. You will need to login with the same credentials used to download the COSMIC dataset. Move this file to the `data/COSMIC` directory of this project. Rename the file `CGC.csv`.

### Tumor Suppressors

Download supplementary table 4 from [PMID: 27911828](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5167163/) and place it in the `data/PMID27911828` folder. It should have the file name `pnas.1616440113.sd04.xlsx`

### CDS

Download CDS coordinates for each genome (~120 MB). Back in the RStudio console, run the following commands.

``` r
library(tidyverse)
library(iSTOP)

CDS_Athaliana_BioMart_plantsmart28() %>% 
  write_csv('data/CDS/Athaliana-plantsmart28.csv')
CDS_Celegans_UCSC_ce11() %>% 
  write_csv('data/CDS/Celegans-ce11.csv')
CDS_Dmelanogaster_UCSC_dm6() %>% 
  write_csv('data/CDS/Dmelanogaster-dm6.csv')
CDS_Drerio_UCSC_danRer10() %>% 
  write_csv('data/CDS/Drerio-danRer10.csv')
CDS_Hsapiens_UCSC_hg38() %>% 
  write_csv('data/CDS/Hsapiens-hg38.csv')
CDS_Mmusculus_UCSC_mm10() %>% 
  write_csv('data/CDS/Mmusculus-mm10.csv')
CDS_Rnorvegicus_UCSC_rn6() %>% 
  write_csv('data/CDS/Rnorvegicus-rn6.csv')
CDS_Scerevisiae_UCSC_sacCer3() %>% 
  write_csv('data/CDS/Scerevisiae-sacCer3.csv')
```

Comprehensive search for iSTOP targetable sites
===============================================

Load all required packages and functions with the following commands.

``` r
library(tidyverse)
library(stringr)
library(iSTOP)

# Source all R functions defined in this project
list.files('R/functions', '[.]R$', full.names = T) %>% walk(source)
```

Given CDS coordinates and genomes, search for all iSTOP sites with the following commands. Raw results will be saved to the `data/iSTOP` directory, and a compacted version with RFLP annotations will be saved to the `data/iSTOP-compact` directory. To dramatically reduce computation time, comment out the `add_RFLP` lines. Documentation for each function can be accessed using `?`. For example, to read detailed documentation for the `locate_codons` function, enter `?locate_codons` in the R console. In brief, previously downloaded CDS coordinates are passed to `locate_codons` which validates each transcript and determines genomic coordinates of user specified codons (CAA, CAG, CGA and TGG by default). These results are then passed to `locate_PAM`, which extracts genomic sequence context and searches for an appropriately spaced PAM (NGG, NGA, NGAG, NGCG, NNGRRT and NNNRRT by default). This dataset is saved, and then a compact version is constructed such that each targetable coordinate within a gene is represented by a single row. Finally, RFLP annotations are added with varying widths of unique cutting.

``` r
# Adjust the number of cores for parallel computation to suit your computer
# Assume that each core will require ~2.5 GB of Memory
# Set to 1 or 2 if you are unsure. Only 1 core is supported on Windows
cores = 1

# Only the human datasets are necessary to reproduce Figure 3
read_csv('data/CDS/Hsapiens-hg38.csv', col_types = 'cciccii') %>%
  locate_codons(BSgenome.Hsapiens.UCSC.hg38::Hsapiens, cores = cores) %>%
  locate_PAM(BSgenome.Hsapiens.UCSC.hg38::Hsapiens) %>%
  write_csv('data/iSTOP/Hsapiens-hg38.csv') %>%       # ~ 1 GB
  compact_iSTOP() %>%  # compresses information for multiple transcripts
  # Loss of cut site columns named RFLP_C_<width>
  add_RFLP(width = 150, cores = cores) %>%
  add_RFLP(width = 100, cores = cores) %>%
  add_RFLP(width = 50,  cores = cores) %>%
  # Gain of cut site columns named RFLP_T_<width>
  add_RFLP(width = 150, cores = cores, recognizes = 't') %>%
  add_RFLP(width = 100, cores = cores, recognizes = 't') %>%
  add_RFLP(width = 50,  cores = cores, recognizes = 't') %>%
  write_csv('data/iSTOP-compact/Hsapiens-hg38.csv')

# This script will run the above command for the non-human species
# Only RFLP width 50 is computed for these species
source('R/scripts/iSTOP-non-human-species.R')
```

Once complete, compute summaries by codon and untargetable datasets with the following commands.

``` r
# Writes two files each to 
# `data/iSTOP-by-codon` and `data/iSTOP-untargetable`
list.files('data/iSTOP', '[.]csv$', full.names = T) %>% 
  pbapply::pblapply(summarize_by_codon, cl = cores)

# Summarize targetability for all species on codon, ORF and gene levels
source('R/scripts/summarize-by-codon-ORF-gene.R')
```

Off-target estimates
====================

**Warning:** No memory optimization has been done for off-target searching, so this analysis requires &gt;6 GB of RAM and several hours of processing time. This section can be skipped if desired.

The following script will populate the `data/Off-target-counts` directory with estimates of the number of targets in the genome for each guide. The current settings restrict the search space by fixing positions 9 through 12 in the guide (i.e. no ambiguities are allowed in the guide's "seed" sequence). The PAM is also fixed allowing for ambiguity where specified (e.g. "NGG" allows ambiguity in the first position of the PAM). Up to two mismatches are tolerated in positions 1 through 8 of the guide. The resulting tables include two columns:

1.  **guide** - The guide sequence
2.  **n\_fuzzy\_matches** - The number of fuzzy matches in the genome. Each guide is expected to have 1 match in the genome. More than 1 match indicates there might be other target sites in the genome (when allowing two mismatches in the leading 8 bases of the guide)

``` r
cores = 1
source('R/scripts/Off-target-count-estimation.R')
```

Analysis of COSMIC nonsense mutations
=====================================

The raw COSMIC dataset can be cleaned and summarized by sourcing the `R/Clean-COSMIC.R` script. This will add three datasets to the `data/COSMIC` directory.

1.  `COSMIC-iSTOP.csv` - (~360 MB) All frameshift and substitution mutations for GRCh38, with aggregated cancer types, and annotated as to whether or not the mutation corresponds to an iSTOP targetable coordinate.
2.  `COSMIC-nonsense.csv` - (~7 MB) Only nonsense mutations from `COSMIC-iSTOP.csv`
3.  `COSMIC-summary-by-cancer.csv` - Summary by cancer type that details the frequency of nonsense, and targetability with iSTOP
4.  `COSMIC-summary-by-gene.csv` - Summary by gene that includes test results for frequent stoppers (likely tumor suppressors). Note that test results for "All cancers" are simply the smallest observed p-value across all cancer subtypes for a given gene.

``` r
source('R/scripts/Clean-COSMIC.R')
```

Cancer subtypes are defined by the following:

``` r
# Note: this is computed during source('R/scripts/Clean-COSMIC.R')
case_when(
  # Case definition                                              ~ Case name
  str_detect(primary_histology, 'melanoma')                      ~ 'Malignant melanoma',
  primary_site == 'large_intestine'                              ~ 'Colorectal',
  primary_site == 'endometrium'                                  ~ 'Endometrial',
  primary_site == 'lung'                                         ~ 'Lung',
  primary_site == 'liver'                                        ~ 'Liver',
  primary_site == 'skin'                                         ~ 'Non-melanoma skin',
  primary_site == 'breast'                                       ~ 'Breast',
  primary_site == 'stomach'                                      ~ 'Stomach',
  primary_site %in% c('upper_aerodigestive_tract', 'oesophagus') ~ 'Upper aerodigestive',
  primary_site == 'haematopoietic_and_lymphoid_tissue'           ~ 'Blood',
  primary_site == 'prostate'                                     ~ 'Prostate',
  primary_site == 'pancreas'                                     ~ 'Pancreatic',
  primary_site == 'urinary_tract'                                ~ 'Bladder',
  primary_site == 'kidney'                                       ~ 'Kidney',
  primary_histology == 'glioma'                                  ~ 'Glioma',
  primary_site == 'ovary'                                        ~ 'Ovarian',
  primary_site == 'cervix'                                       ~ 'Cervical',
  primary_site == 'thyroid'                                      ~ 'Thyroid',
  primary_site == 'bone'                                         ~ 'Bone',
  TRUE ~ 'Other'
)
```

Binomial p-values and multiple test corrected q-values for the frequent iSTOPer analysis were computed with the following:

``` r
p = binom.test(
  n_iSTOP_sites_in_gene_in_cancer, # Successes in gene in cancer type
  n_iSTOP_sites_in_gene,           # Trials in gene
  P_event[cancer_type],            # Probability of success in cancer type
  alternative = 'greater'          # One-sided test
)$p.value

q = p.adjust(p, method = 'fdr')

# To read the documentation on these functions
help(binom.test)
help(p.adjust)
```

Reproducing Figures
===================

The following script will write figures to the `figures` directory of the project. Final figures for the paper were edited in [Inkscape](https://inkscape.org/en/) to reduce the size of the files and improve readability of figure legends and axis labels.

``` r
source('R/figures/build-all.R')
```

Session Information
===================

This analysis was successfully performed with the following system, and package versions:

    ##  setting  value                       
    ##  version  R version 3.4.0 (2017-04-21)
    ##  system   x86_64, darwin15.6.0        
    ##  ui       X11                         
    ##  language (EN)                        
    ##  collate  en_US.UTF-8                 
    ##  tz       America/New_York            
    ##  date     2017-07-18                  
    ## 
    ##  package              * version  date       source                          
    ##  assertthat             0.2.0    2017-04-11 CRAN (R 3.4.0)                  
    ##  backports              1.1.0    2017-05-22 CRAN (R 3.4.0)                  
    ##  base                 * 3.4.0    2017-04-21 local                           
    ##  bindr                  0.1      2016-11-13 CRAN (R 3.4.0)                  
    ##  bindrcpp               0.2      2017-06-17 CRAN (R 3.4.0)                  
    ##  Biobase                2.36.2   2017-05-04 Bioconductor                    
    ##  BiocGenerics           0.22.0   2017-04-25 cran (@0.22.0)                  
    ##  BiocParallel           1.10.1   2017-05-03 Bioconductor                    
    ##  Biostrings             2.44.1   2017-06-01 Bioconductor                    
    ##  bitops                 1.0-6    2013-08-17 CRAN (R 3.4.0)                  
    ##  broom                  0.4.2    2017-02-13 CRAN (R 3.4.0)                  
    ##  BSgenome               1.44.0   2017-04-25 Bioconductor                    
    ##  cellranger             1.1.0    2016-07-27 CRAN (R 3.4.0)                  
    ##  colorspace             1.3-2    2016-12-14 CRAN (R 3.4.0)                  
    ##  compiler               3.4.0    2017-04-21 local                           
    ##  datasets             * 3.4.0    2017-04-21 local                           
    ##  DelayedArray           0.2.7    2017-06-03 Bioconductor                    
    ##  devtools               1.13.2   2017-06-02 CRAN (R 3.4.0)                  
    ##  digest                 0.6.12   2017-01-27 CRAN (R 3.4.0)                  
    ##  dplyr                * 0.7.1    2017-06-22 CRAN (R 3.4.0)                  
    ##  evaluate               0.10.1   2017-06-24 CRAN (R 3.4.0)                  
    ##  forcats                0.2.0    2017-01-23 CRAN (R 3.4.0)                  
    ##  foreign                0.8-69   2017-06-21 CRAN (R 3.4.0)                  
    ##  fuzzyjoin              0.1.3    2017-06-19 CRAN (R 3.4.1)                  
    ##  GenomeInfoDb           1.12.2   2017-06-09 Bioconductor                    
    ##  GenomeInfoDbData       0.99.0   2017-04-25 Bioconductor                    
    ##  GenomicAlignments      1.12.1   2017-05-12 Bioconductor                    
    ##  GenomicRanges          1.28.3   2017-05-25 Bioconductor                    
    ##  ggplot2              * 2.2.1    2016-12-30 CRAN (R 3.4.0)                  
    ##  glue                   1.1.1    2017-06-21 CRAN (R 3.4.0)                  
    ##  graphics             * 3.4.0    2017-04-21 local                           
    ##  grDevices            * 3.4.0    2017-04-21 local                           
    ##  grid                   3.4.0    2017-04-21 local                           
    ##  gtable                 0.2.0    2016-02-26 CRAN (R 3.4.0)                  
    ##  haven                  1.0.0    2016-09-23 CRAN (R 3.4.0)                  
    ##  hms                    0.3      2016-11-22 CRAN (R 3.4.0)                  
    ##  htmltools              0.3.6    2017-04-28 CRAN (R 3.4.0)                  
    ##  httr                   1.2.1    2016-07-03 CRAN (R 3.4.0)                  
    ##  IRanges                2.10.2   2017-05-25 Bioconductor                    
    ##  iSTOP                * 0.1.0    2017-07-06 Github (CicciaLab/iSTOP@dfa6e40)
    ##  jsonlite               1.5      2017-06-01 CRAN (R 3.4.0)                  
    ##  knitr                  1.16     2017-05-18 CRAN (R 3.4.0)                  
    ##  lattice                0.20-35  2017-03-25 CRAN (R 3.4.0)                  
    ##  lazyeval               0.2.0    2016-06-12 CRAN (R 3.4.0)                  
    ##  lubridate              1.6.0    2016-09-13 CRAN (R 3.4.0)                  
    ##  magrittr               1.5      2014-11-22 CRAN (R 3.4.0)                  
    ##  Matrix                 1.2-10   2017-04-28 CRAN (R 3.4.0)                  
    ##  matrixStats            0.52.2   2017-04-14 CRAN (R 3.4.0)                  
    ##  memoise                1.1.0    2017-04-21 CRAN (R 3.4.0)                  
    ##  methods              * 3.4.0    2017-04-21 local                           
    ##  mnormt                 1.5-5    2016-10-15 CRAN (R 3.4.0)                  
    ##  modelr                 0.1.0    2016-08-31 CRAN (R 3.4.0)                  
    ##  munsell                0.4.3    2016-02-13 CRAN (R 3.4.0)                  
    ##  nlme                   3.1-131  2017-02-06 CRAN (R 3.4.0)                  
    ##  parallel               3.4.0    2017-04-21 local                           
    ##  pbapply                1.3-3    2017-07-04 CRAN (R 3.4.1)                  
    ##  pkgconfig              2.0.1    2017-03-21 CRAN (R 3.4.0)                  
    ##  plyr                   1.8.4    2016-06-08 CRAN (R 3.4.0)                  
    ##  psych                  1.7.5    2017-05-03 CRAN (R 3.4.0)                  
    ##  purrr                * 0.2.2.2  2017-05-11 cran (@0.2.2.2)                 
    ##  R6                     2.2.2    2017-06-17 CRAN (R 3.4.0)                  
    ##  Rcpp                   0.12.11  2017-05-22 cran (@0.12.11)                 
    ##  RCurl                  1.95-4.8 2016-03-01 CRAN (R 3.4.0)                  
    ##  readr                * 1.1.1    2017-05-16 cran (@1.1.1)                   
    ##  readxl                 1.0.0    2017-04-18 CRAN (R 3.4.0)                  
    ##  reshape2               1.4.2    2016-10-22 CRAN (R 3.4.0)                  
    ##  rlang                  0.1.1    2017-05-18 cran (@0.1.1)                   
    ##  rmarkdown              1.6      2017-06-15 CRAN (R 3.4.0)                  
    ##  rprojroot              1.2      2017-01-16 CRAN (R 3.4.0)                  
    ##  Rsamtools              1.28.0   2017-04-25 Bioconductor                    
    ##  rtracklayer            1.36.3   2017-05-25 Bioconductor                    
    ##  rvest                  0.3.2    2016-06-17 CRAN (R 3.4.0)                  
    ##  S4Vectors              0.14.3   2017-06-03 Bioconductor                    
    ##  scales                 0.4.1    2016-11-09 CRAN (R 3.4.0)                  
    ##  stats                * 3.4.0    2017-04-21 local                           
    ##  stats4                 3.4.0    2017-04-21 local                           
    ##  stringi                1.1.5    2017-04-07 CRAN (R 3.4.0)                  
    ##  stringr              * 1.2.0    2017-02-18 CRAN (R 3.4.0)                  
    ##  SummarizedExperiment   1.6.3    2017-05-29 Bioconductor                    
    ##  tibble               * 1.3.3    2017-05-28 cran (@1.3.3)                   
    ##  tidyr                * 0.6.3    2017-05-15 cran (@0.6.3)                   
    ##  tidyverse            * 1.1.1    2017-01-27 CRAN (R 3.4.0)                  
    ##  tools                  3.4.0    2017-04-21 local                           
    ##  utils                * 3.4.0    2017-04-21 local                           
    ##  withr                  1.0.2    2016-06-20 CRAN (R 3.4.0)                  
    ##  XML                    3.98-1.9 2017-06-19 CRAN (R 3.4.1)                  
    ##  xml2                   1.1.1    2017-01-24 CRAN (R 3.4.0)                  
    ##  XVector                0.16.0   2017-04-25 cran (@0.16.0)                  
    ##  yaml                   2.1.14   2016-11-12 CRAN (R 3.4.0)                  
    ##  zlibbioc               1.22.0   2017-04-25 cran (@1.22.0)
