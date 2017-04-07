# COSMIC

This directory should contain data downloaded from [COSMIC](https://cancer.sanger.ac.uk/cosmic).

## Downloaded datasets

### Mutation data

Download the "COSMIC Mutation Data" (~300 MB compressed) from the [Catalogue of Somatic Mutations in Cancer](https://cancer.sanger.ac.uk/cosmic/download). This dataset requires [registration with a valid email address](https://cancer.sanger.ac.uk/cosmic/register), then in your terminal (not the R console!) you can download the file with the following commands (substitute `your_email_address` with the email used to register).

```
sftp "your_email_address"@sftp-cancer.sanger.ac.uk
# You will be prompted for the password you provided when you registered your email address
get /files/grch38/cosmic/v80/CosmicMutantExport.tsv.gz
```

Move this file to the `data/COSMIC` directory of this project. It should already be named `CosmicMutantExport.tsv.gz`.

Columns used in this study are:

- `ID_sample`
- `ID_tumour`
- `Primary site`
- `Mutation ID`
- `Primary histology`
- `Mutation zygosity`
- `Mutation Description`
- `Mutation CDS`
- `Mutation AA`
- `GRCh`
- `Mutation genome position`
- `Mutation strand`
- `FATHMM prediction`
- `FATHMM score`
- `Mutation somatic status`

### CGC

Download the [Cancer Gene Census (CGC) dataset](http://cancer.sanger.ac.uk/census) by clicking on the `CSV` Export button. You will need to login with the same credentials used to download the COSMIC dataset. Move this file to the `data/COSMIC` directory of this project. Rename the file `CGC.csv`.

The CGC dataset has the following columns (more details can be found at the COSMIC website):

- `Gene Symbol`
- `Name`
- `Entrez GeneId`
- `Genome Location`
- `Chr Band`
- `Somatic`
- `Germline`
- `Tumour Types(Somatic)`
- `Tumour Types(Germline)`
- `Cancer Syndrome`
- `Tissue Type`
- `Molecular Genetics`
- `Role in Cancer`
- `Mutation Types`
- `Translocation Partner`
- `Other Germline Mut`
- `Other Syndrome`
- `Synonyms`

## Computed Datasets

The raw COSMIC dataset can be cleaned and summarized by sourcing the `R/Clean-COSMIC.R` script. This will add four datasets to the `data/COSMIC` directory.

1. `COSMIC-iSTOP.csv` - (~360 MB) All frameshift and substitution mutations for GRCh38, with aggregated cancer types, and annotated as to whether or not the mutation corresponds to an iSTOP targetable coordinate.
2. `COSMIC-nonsense.csv` - (~7 MB) Only nonsense mutations from `COSMIC-iSTOP.csv`
3. `COSMIC-summary-by-cancer.csv` - Summary by cancer type that details the frequency of nonsense, and targetability with iSTOP
4. `COSMIC-summary-by-gene.csv` - Summary by gene that includes test results for frequent stoppers (likely tumor suppressors). Note that test results for "All cancers" are simply the smallest observed p-value across all cancer subtypes for a given gene.
