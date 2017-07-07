compact_iSTOP <- function(iSTOP) {
  iSTOP %>%
    filter(match_any) %>%
    group_by(
      # Base edit site information
      gene, chr, strand, sg_strand, aa_target, codon, genome_coord,
      # Details
      n_tx_in_gene, n_tx, percent_tx, no_upstream_G, # match_any, not necessary after filter
      # Genomic Sequence Context
      searched,
      # Guides and guide information
      sgNGG, sgNGG_spacing,
      sgNGA, sgNGA_spacing,
      sgNGCG, sgNGCG_spacing,
      sgNGAG, sgNGAG_spacing,
      sgNNGRRT, sgNNGRRT_spacing,
      sgNNNRRT, sgNNNRRT_spacing
    ) %>%
    summarise(
      txs         = stringr::str_c(tx, collapse = ' | '),
      #exons       = stringr::str_c(exon, collapse = ' | '),
      pep_lengths = stringr::str_c(pep_length, collapse = ' | '),
      #cds_lengths = stringr::str_c(cds_length, collapse = ' | '),
      aa_coords   = stringr::str_c(aa_coord, collapse = ' | '),
      cds_coords  = stringr::str_c(cds_coord, collapse = ' | '),
      percent_NMD = (sum(NMD_pred) / unique(n_tx_in_gene)) * 100,
      #NMD_pred    = stringr::str_c(as.character(NMD_pred), collapse = ' | '),
      rel_pos_largest_isoform = mean(rel_position[which.max(cds_length)])
    ) %>%
    select(
      # Concatenated transcript details
      txs, pep_lengths, aa_coords, cds_coords,
      # Gene/Codon level details
      gene:percent_tx,
      # Prioritization details
      percent_NMD, rel_pos_largest_isoform, no_upstream_G,
      # Guides
      everything()
    )
}
