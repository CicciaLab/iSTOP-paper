summarize_by_codon <- function(path) {
  all <- read_csv(path, col_types = cols(percent_tx = col_number()), progress = F)

  untargetable <-
    all %>%
    group_by(gene, tx, pep_length, cds_length) %>%
    filter(is.na(codon)) %>%
    write_csv(file.path('data', 'iSTOP-untargetable', basename(path)))

  gene_tx_codon_summary <-
    all %>%
    # If aa_coord is NA than the transcript is untargetable. Set it's aa_coord to Inf
    mutate(aa_coord = if_else(is.na(aa_coord), Inf,  as.numeric(aa_coord))) %>%
    # TGG codons have two records (1 for each editable position) summarise these into match no match using aa_coord
    group_by(gene, tx, exon, chr, strand, pep_length, cds_length, codon, aa_coord) %>%
    summarise(
      genome_coord   = mean(genome_coord), # Just keep the first matching coordinate
      match_any      = any(match_any),
      match_sgNGG    = any(!is.na(sgNGG)),
      match_sgNGA    = any(!is.na(sgNGA)),
      match_sgNGCG   = any(!is.na(sgNGCG)),
      match_sgNGAG   = any(!is.na(sgNGAG)),
      match_sgNNGRRT = any(!is.na(sgNNGRRT)),
      match_sgNNNRRT = any(!is.na(sgNNNRRT))
    ) %>%
    write_csv(file.path('data', 'iSTOP-by-codon', basename(path)))
  return(invisible())
}
