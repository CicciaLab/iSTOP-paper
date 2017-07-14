library(tidyverse)
library(stringr)

ultra_compact <- function(compact) {
  cmp <- read_csv(str_c('data/iSTOP-compact/', compact), col_types = cols(percent_tx = col_double(), percent_NMD = col_double(), rel_pos_largest_isoform = col_double()))
  # eliminated columns
  cmp$txs         <- NULL # information available by looking up genome_coord
  cmp$cds_coords  <- NULL # information available by looking up genome_coord
  cmp$aa_target   <- NULL # information available in codon
  cmp$strand      <- NULL  # information available in sg_strand
  cmp$n_tx        <- NULL # information available in percent_tx
  cmp$pep_lengths <- NULL # information available in cds_lengths
  cmp$aa_coords   <- NULL # information available in cds_coords
  cmp$searched    <- NULL
  cmp$RFLP_C_100  <- NULL
  cmp$RFLP_C_150  <- NULL
  cmp$RFLP_T_100  <- NULL
  cmp$RFLP_T_150  <- NULL
  cmp$sgNGG_spacing <- NULL
  cmp$sgNGA_spacing <- NULL
  cmp$sgNGAG_spacing <- NULL
  cmp$sgNGCG_spacing <- NULL
  cmp$sgNNGRRT_spacing <- NULL
  cmp$sgNNNRRT_spacing <- NULL
  cmp$percent_tx <- round(cmp$percent_tx, 1)
  cmp$percent_NMD <- round(cmp$percent_NMD, 1)
  cmp$rel_pos_largest_isoform <- round(cmp$rel_pos_largest_isoform, 2)
  cmp$no_upstream_G <- ifelse(cmp$no_upstream_G, 1, 0)
  write_csv(cmp, str_c('data/iSTOP-ultra-compact/', compact, '.gz'), na = '')
}

list.files('data/iSTOP-compact') %>% walk(ultra_compact)
