# CDS coordinates

The tables in this directory consist of CDS coordinates in CSV format with the following structure.

| Column   | Type      | Description                          |
| :------- | :-------- | :----------------------------------- |
| `tx`     | Character | Transcript ID                        |
| `gene`   | Character | Gene Name                            |
| `exon`   | Integer   | Exon number                          |
| `chr`    | Character | Chromosome name                      |
| `strand` | Character | Strand (+/-)                         |
| `start`  | Integer   | Start coordinate of exon (inclusive) |
| `end`    | Integer   | End coordinate of exon (inclusive)   |

All CDS coordinates are downloaded from UCSC with the exception of *A. thaliana* which were downloaded from Plantsmart28 and BioMart via Bioconductor R packages.
