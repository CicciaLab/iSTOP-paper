# iSTOP targets

The tables in this directory consist of iSTOP targetable sites. Each row represents a gene, transcript, targetable coordinate. Note that `TGG` codons are targetable at two coordinates and therefore will have a separate entry for each coordinate.

| Column       | Type      | Description                          |
| :----------- | :-------- | :----------------------------------- |
| `tx`         | Character | Transcript ID                        |
| `gene`       | Character | Gene Name                            |
| `exon`       | Integer   | Exon number                          |
| `pep_length` | Integer   | Length of peptide                    |
| `cds_length` | Integer   | Length of coding sequence            |
| `chr`        | Character | Chromosome name                      |
| `strand`     | Character | Strand (+/-)                         |
| `sg_strand`  | Character | Guide strand (+/-)                   |
| `aa_target`  | Character | Amino Acid target (Q/R/W)            |
| `codon`      | Character | Codon (CAA/CAG/CGA/TGG)              |
| `aa_coord`   | Integer   | Amino Acid coordinate                |
| `cds_coord`  | Integer   | CDS coordinate                       |
| `genome_coord` | Integer | Genome coordinate                    |
| `NMD_pred`   | Logical   | NMD prediction (T/F)                 |
| `n_tx_in_gene` | Integer | Number of transcript isoforms in gene |
| `n_tx`       | Integer   | Number of transcript isoforms targeted at this coordinate |
| `percent_tx` | Numeric   | Percent of isoforms targeted at this coordinate |
| `searched`   | Character | Sequence context of the targeted base +/- 150 bp |
| `sgNGG`      | Character | sg sequence for NGG PAM (NA of not available) |
| `sgNGA`      | Character | sg sequence for NGA PAM (NA of not available) |
| `sgNGCG`     | Character | sg sequence for NGCG PAM (NA of not available) |
| `sgNGAG`     | Character | sg sequence for NGAG PAM (NA of not available) |
| `sgNNGRRT`   | Character | sg sequence for NNGRRT PAM (NA of not available) |
| `sgNNNRRT`   | Character | sg sequence for NNNRRT PAM (NA of not available) |
| `match_any`  | Logical   | Is there an sgRNA for any of the PAMs? |


Read these tables into an R session with:

```r
readr::read_csv('path/to/iSTOP.csv', col_types = '')
```
