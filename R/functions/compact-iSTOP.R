compact_iSTOP <- function(iSTOP) {

  compact <-
    iSTOP %>%
    filter(match_any) %>%
    group_by(
      gene, chr, strand, sg_strand, aa_target, codon, genome_coord, n_tx_in_gene, n_tx,
      percent_tx, searched, sgNGG, sgNGA, sgNGCG, sgNGAG, sgNNGRRT, sgNNNRRT
    ) %>%
    summarise(
      txs         = stringr::str_c(tx, collapse = ' | '),
      exons       = stringr::str_c(exon, collapse = ' | '),
      pep_lengths = stringr::str_c(pep_length, collapse = ' | '),
      cds_lengths = stringr::str_c(cds_length, collapse = ' | '),
      aa_coords   = stringr::str_c(aa_coord, collapse = ' | '),
      cds_coords  = stringr::str_c(cds_coord, collapse = ' | '),
      NMD_pred    = stringr::str_c(as.character(NMD_pred), collapse = ' | ')
    ) %>%
    select(txs, everything())
}
