library(tidyverse)
library(stringr)
ultra_compact <- function(compact) {
  cmp <- read_csv(str_c('data/iSTOP-compact/', compact), col_types = cols(percent_tx = col_double()))
  # eliminated columns
  cmp$aa_target <- NULL # information available in codon
  cmp$strand   <- NULL  # information available in sg_strand
  cmp$n_tx     <- NULL # information available in percent_tx
  cmp$pep_lengths <- NULL # information available in cds_lengths
  cmp$aa_coords <- NULL # information available in cds_coords
  cmp$searched <- NULL
  cmp$RFLP_100 <- NULL
  cmp$RFLP_150 <- NULL
  write_csv(cmp, str_c('data/iSTOP-ultra-compact/', compact, '.gz'))
}

list.files('data/iSTOP-compact') %>% walk(ultra_compact)

#x <- read_csv('data/iSTOP-compact/Hsapiens-hg38.csv')
