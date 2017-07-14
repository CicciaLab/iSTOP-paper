library(tidyverse)
library(stringr)

website_version <- function(compact) {
  cmp <- read_csv(str_c('data/iSTOP-compact/', compact), col_types = cols(percent_tx = col_double(), percent_NMD = col_double(), rel_pos_largest_isoform = col_double()))
  # eliminated columns
  #cmp$txs         <- NULL # information available by looking up genome_coord
  #cmp$cds_coords  <- NULL # information available by looking up genome_coord
  #cmp$aa_target   <- NULL # information available in codon
  #cmp$strand      <- NULL  # information available in sg_strand
  cmp$n_tx        <- NULL # information available in percent_tx
  #cmp$pep_lengths <- NULL # information available in cds_lengths
  cmp$aa_coords   <- NULL # information available in cds_coords
  #cmp$searched    <- NULL
  cmp$RFLP_C_100  <- NULL
  cmp$RFLP_C_150  <- NULL
  cmp$RFLP_T_100  <- NULL
  cmp$RFLP_T_150  <- NULL
  #cmp$sgNGG_spacing <- NULL
  #cmp$sgNGA_spacing <- NULL
  #cmp$sgNGAG_spacing <- NULL
  #cmp$sgNGCG_spacing <- NULL
  #cmp$sgNNGRRT_spacing <- NULL
  #cmp$sgNNNRRT_spacing <- NULL
  cmp$percent_tx  <- round(cmp$percent_tx, 2)
  cmp$percent_NMD <- round(cmp$percent_NMD, 2)
  cmp$rel_pos_largest_isoform <- round(cmp$rel_pos_largest_isoform, 3)
  #cmp$no_upstream_G <- ifelse(cmp$no_upstream_G, 1, 0)
  cmp %>%
    select(
      # Coordinate annotations
      txs, pep_lengths, cds_coords, gene, chr, strand, sg_strand, genome_coord, aa_target, codon, n_tx_in_gene,
      # Prioritizations
      percent_tx, percent_NMD, rel_pos_largest_isoform, no_upstream_G, RFLP_Loss = RFLP_C_50, RFLP_Gain = RFLP_T_50,
      # Guides
      starts_with('sg'),
      genomic_sequence_context = searched
    ) %>%
    write_csv(str_c('data/iSTOP-website-version/', compact, '.gz'), na = '')
}

list.files('data/iSTOP-compact') %>% walk(website_version)


# ---- Special human dataset ----

# Cancer Data
nonsense <-
  read_csv('data/COSMIC/COSMIC-nonsense.csv') %>%
  filter(iSTOP) %>%
  select(gene, chr, strand, genome_coord, cancer_type) %>%
  group_by(gene, chr, strand, genome_coord) %>%
  summarise(cancer_type = stringr::str_c(unique(cancer_type), collapse = ' | ')) %>%
  ungroup %>%
  distinct #%>%
  #mutate(COSMIC_nonsense = TRUE)

read_csv('data/iSTOP-website-version/Hsapiens-hg38.csv.gz', col_types = cols()) %>%
  select(-matches('_matches$'), -cancer_type) %>%
  left_join(nonsense, by = c("gene", "chr", "strand", "genome_coord")) %>%
  mutate(guide = toupper(sgNGA)) %>%
  left_join(
    read_csv('data/Off-target-counts/Human-UCSC-hg38-NGA.csv', col_types = 'ci') %>%
      rename(sgNGA_matches = n_fuzzy_match) %>%
      mutate(guide = str_sub(guide, start = 1, end = 20)),
    by = 'guide'
  ) %>%
  mutate(guide = toupper(sgNGG)) %>%
  left_join(
    read_csv('data/Off-target-counts/Human-UCSC-hg38-NGG.csv', col_types = 'ci') %>%
      rename(sgNGG_matches = n_fuzzy_match) %>%
      mutate(guide = str_sub(guide, start = 1, end = 20)),
    by = 'guide'
  ) %>%
  mutate(guide = toupper(sgNGCG)) %>%
  left_join(
    read_csv('data/Off-target-counts/Human-UCSC-hg38-NGCG.csv', col_types = 'ci') %>%
      rename(sgNGCG_matches = n_fuzzy_match) %>%
      mutate(guide = str_sub(guide, start = 1, end = 20)),
    by = 'guide'
  ) %>%
  mutate(guide = toupper(sgNGAG)) %>%
  left_join(
    read_csv('data/Off-target-counts/Human-UCSC-hg38-NGAG.csv', col_types = 'ci') %>%
      rename(sgNGAG_matches = n_fuzzy_match) %>%
      mutate(guide = str_sub(guide, start = 1, end = 20)),
    by = 'guide'
  ) %>%
  mutate(guide = toupper(sgNNGRRT)) %>%
  left_join(
    read_csv('data/Off-target-counts/Human-UCSC-hg38-NNGRRT.csv', col_types = 'ci') %>%
      rename(sgNNGRRT_matches = n_fuzzy_match) %>%
      mutate(guide = str_sub(guide, start = 1, end = 20)),
    by = 'guide'
  ) %>%
  mutate(guide = toupper(sgNNNRRT)) %>%
  left_join(
    read_csv('data/Off-target-counts/Human-UCSC-hg38-NNNRRT.csv', col_types = 'ci') %>%
      rename(sgNNNRRT_matches = n_fuzzy_match) %>%
      mutate(guide = str_sub(guide, start = 1, end = 20)),
    by = 'guide'
  ) %>%
  select(
    txs:RFLP_Gain, cancer_type,
    starts_with('sgNGG'),
    sgNGA,
    starts_with('sgNGA_'),
    starts_with('sgNGCG'),
    starts_with('sgNGAG'),
    starts_with('sgNNGRRT'),
    starts_with('sgNNNRRT'),
    genomic_sequence_context
  ) %>%
  write_csv('data/iSTOP-website-version/Hsapiens-hg38.csv.gz', na = '')

# ---- Add Off target counts ----

add_off_targets <- function(path_main, path_NGA, path_NGG, path_NGAG, 
                            path_NGCG, path_NNGRRT, path_NNNRRT) {
  read_csv(path_main, col_types = cols(percent_tx = col_double(), percent_NMD = col_double(), rel_pos_largest_isoform = col_double())) %>%
    select(-matches('_matches$')) %>%
    mutate(guide = toupper(sgNGA)) %>%
    left_join(
      read_csv(path_NGA, col_types = 'ci') %>%
        rename(sgNGA_matches = n_fuzzy_match) %>%
        mutate(guide = str_sub(guide, start = 1, end = 20)),
      by = 'guide'
    ) %>%
    mutate(guide = toupper(sgNGG)) %>%
    left_join(
      read_csv(path_NGG, col_types = 'ci') %>%
        rename(sgNGG_matches = n_fuzzy_match) %>%
        mutate(guide = str_sub(guide, start = 1, end = 20)),
      by = 'guide'
    ) %>%
    mutate(guide = toupper(sgNGCG)) %>%
    left_join(
      read_csv(path_NGCG, col_types = 'ci') %>%
        rename(sgNGCG_matches = n_fuzzy_match) %>%
        mutate(guide = str_sub(guide, start = 1, end = 20)),
      by = 'guide'
    ) %>%
    mutate(guide = toupper(sgNGAG)) %>%
    left_join(
      read_csv(path_NGAG, col_types = 'ci') %>%
        rename(sgNGAG_matches = n_fuzzy_match) %>%
        mutate(guide = str_sub(guide, start = 1, end = 20)),
      by = 'guide'
    ) %>%
    mutate(guide = toupper(sgNNGRRT)) %>%
    left_join(
      read_csv(path_NNGRRT, col_types = 'ci') %>%
        rename(sgNNGRRT_matches = n_fuzzy_match) %>%
        mutate(guide = str_sub(guide, start = 1, end = 20)),
      by = 'guide'
    ) %>%
    mutate(guide = toupper(sgNNNRRT)) %>%
    left_join(
      read_csv(path_NNNRRT, col_types = 'ci') %>%
        rename(sgNNNRRT_matches = n_fuzzy_match) %>%
        mutate(guide = str_sub(guide, start = 1, end = 20)),
      by = 'guide'
    ) %>%
    select(
      txs:RFLP_Gain,
      starts_with('sgNGG'),
      sgNGA,
      starts_with('sgNGA_'),
      starts_with('sgNGCG'),
      starts_with('sgNGAG'),
      starts_with('sgNNGRRT'),
      starts_with('sgNNNRRT'),
      genomic_sequence_context
    ) %>%
    write_csv(path_main, na = '')
}

# ---- Arabadopsis ----
add_off_targets(
  path_main   = 'data/iSTOP-website-version/Athaliana-plantsmart28.csv.gz',
  path_NGA    = 'data/Off-target-counts/Arabidopsis-TAIR-TAIR9-NGA.csv',
  path_NGG    = 'data/Off-target-counts/Arabidopsis-TAIR-TAIR9-NGG.csv',
  path_NGAG   = 'data/Off-target-counts/Arabidopsis-TAIR-TAIR9-NGAG.csv',
  path_NGCG   = 'data/Off-target-counts/Arabidopsis-TAIR-TAIR9-NGCG.csv',
  path_NNGRRT = 'data/Off-target-counts/Arabidopsis-TAIR-TAIR9-NNGRRT.csv',
  path_NNNRRT = 'data/Off-target-counts/Arabidopsis-TAIR-TAIR9-NNNRRT.csv'
)

# ---- Celegans ----
add_off_targets(
  path_main   = 'data/iSTOP-website-version/Celegans-ce11.csv.gz',
  path_NGA    = 'data/Off-target-counts/Worm-UCSC-ce11-NGA.csv',
  path_NGG    = 'data/Off-target-counts/Worm-UCSC-ce11-NGG.csv',
  path_NGAG   = 'data/Off-target-counts/Worm-UCSC-ce11-NGAG.csv',
  path_NGCG   = 'data/Off-target-counts/Worm-UCSC-ce11-NGCG.csv',
  path_NNGRRT = 'data/Off-target-counts/Worm-UCSC-ce11-NNGRRT.csv',
  path_NNNRRT = 'data/Off-target-counts/Worm-UCSC-ce11-NNNRRT.csv'
)

# ---- Dmelanogaster ----
add_off_targets(
  path_main   = 'data/iSTOP-website-version/Dmelanogaster-dm6.csv.gz',
  path_NGA    = 'data/Off-target-counts/Fly-UCSC-dm6-NGA.csv',
  path_NGG    = 'data/Off-target-counts/Fly-UCSC-dm6-NGG.csv',
  path_NGAG   = 'data/Off-target-counts/Fly-UCSC-dm6-NGAG.csv',
  path_NGCG   = 'data/Off-target-counts/Fly-UCSC-dm6-NGCG.csv',
  path_NNGRRT = 'data/Off-target-counts/Fly-UCSC-dm6-NNGRRT.csv',
  path_NNNRRT = 'data/Off-target-counts/Fly-UCSC-dm6-NNNRRT.csv'
)

# ---- Drerio ----
add_off_targets(
  path_main   = 'data/iSTOP-website-version/Drerio-danRer10.csv.gz',
  path_NGA    = 'data/Off-target-counts/Zebrafish-UCSC-danRer10-NGA.csv',
  path_NGG    = 'data/Off-target-counts/Zebrafish-UCSC-danRer10-NGG.csv',
  path_NGAG   = 'data/Off-target-counts/Zebrafish-UCSC-danRer10-NGAG.csv',
  path_NGCG   = 'data/Off-target-counts/Zebrafish-UCSC-danRer10-NGCG.csv',
  path_NNGRRT = 'data/Off-target-counts/Zebrafish-UCSC-danRer10-NNGRRT.csv',
  path_NNNRRT = 'data/Off-target-counts/Zebrafish-UCSC-danRer10-NNNRRT.csv'
)

# ---- Mmusculus ----
add_off_targets(
  path_main   = 'data/iSTOP-website-version/Mmusculus-mm10.csv.gz',
  path_NGA    = 'data/Off-target-counts/Mouse-UCSC-mm10-NGA.csv',
  path_NGG    = 'data/Off-target-counts/Mouse-UCSC-mm10-NGG.csv',
  path_NGAG   = 'data/Off-target-counts/Mouse-UCSC-mm10-NGAG.csv',
  path_NGCG   = 'data/Off-target-counts/Mouse-UCSC-mm10-NGCG.csv',
  path_NNGRRT = 'data/Off-target-counts/Mouse-UCSC-mm10-NNGRRT.csv',
  path_NNNRRT = 'data/Off-target-counts/Mouse-UCSC-mm10-NNNRRT.csv'
)

# ---- Rnorvegicus ----
add_off_targets(
  path_main   = 'data/iSTOP-website-version/Rnorvegicus-rn6.csv.gz',
  path_NGA    = 'data/Off-target-counts/Rat-UCSC-rn6-NGA.csv',
  path_NGG    = 'data/Off-target-counts/Rat-UCSC-rn6-NGG.csv',
  path_NGAG   = 'data/Off-target-counts/Rat-UCSC-rn6-NGAG.csv',
  path_NGCG   = 'data/Off-target-counts/Rat-UCSC-rn6-NGCG.csv',
  path_NNGRRT = 'data/Off-target-counts/Rat-UCSC-rn6-NNGRRT.csv',
  path_NNNRRT = 'data/Off-target-counts/Rat-UCSC-rn6-NNNRRT.csv'
)

# ---- Scerevisiae ----
add_off_targets(
  path_main   = 'data/iSTOP-website-version/Scerevisiae-sacCer3.csv.gz',
  path_NGA    = 'data/Off-target-counts/Yeast-UCSC-sacCer3-NGA.csv',
  path_NGG    = 'data/Off-target-counts/Yeast-UCSC-sacCer3-NGG.csv',
  path_NGAG   = 'data/Off-target-counts/Yeast-UCSC-sacCer3-NGAG.csv',
  path_NGCG   = 'data/Off-target-counts/Yeast-UCSC-sacCer3-NGCG.csv',
  path_NNGRRT = 'data/Off-target-counts/Yeast-UCSC-sacCer3-NNGRRT.csv',
  path_NNNRRT = 'data/Off-target-counts/Yeast-UCSC-sacCer3-NNNRRT.csv'
)

