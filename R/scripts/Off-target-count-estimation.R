guides <- readr::read_csv('data/iSTOP-compact/Hsapiens-hg38.csv', col_types = cols())
genome <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
chromosomes <- BSgenome::seqnames(genome)[1:24]

PAMs <- c(
  'NGG'    = '.GG$',
  'NGA'    = '.GA$',
  'NGCG'   = '.GCG$',
  'NGAG'   = '.GAG$',
  'NNGRRT' = '..G[AG][AG]T$',
  'NNNRRT' = '...[AG][AG]T$'
)


walk(1:length(PAMs), function(i) {
  IUPAC <- names(PAMs)[i]
  regex <- PAMs[i]
  column <- stringr::str_c('sg', IUPAC)
  outfile <- stringr::str_c(
    'data/Off-target-counts/',
    BSgenome::commonName(genome), '-',
    BSgenome::provider(genome), '-',
    BSgenome::providerVersion(genome), '-', IUPAC, '.csv')

  guides_1_PAM <- stringr::str_c(unique(toupper(na.omit(guides[[column]]))), IUPAC)


  iSTOP::search_off_target(
    guides_1_PAM,
    genome,
    chromosomes,
    fixed_start  = 9,
    fixed_end    = 20,
    max_mismatch = 2,
    secondary_filter = regex,
    cores = cores
  ) %>%
    group_by(guide) %>%
    summarise(n_fuzzy_match = n()) %>%
    ungroup %>%
    #filter(n_fuzzy_match > 1) %>%
    arrange(desc(n_fuzzy_match)) %>%
    write_csv(outfile)
})
