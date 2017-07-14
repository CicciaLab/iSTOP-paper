#library(tidyverse)
#library(iSTOP)
#cores = 6

PAMs <- c(
  'NGG'    = '.GG$',
  'NGA'    = '.GA$',
  'NGCG'   = '.GCG$',
  'NGAG'   = '.GAG$',
  'NNGRRT' = '..G[AG][AG]T$',
  'NNNRRT' = '...[AG][AG]T$'
)

# ---- Count off targets in genome ----
count_off_target <- function(dir = 'data/Off-target-counts', guides, genome, chromosomes, PAMs) {

  walk(1:length(PAMs), function(i) {
    IUPAC <- names(PAMs)[i]
    regex <- PAMs[i]
    column <- stringr::str_c('sg', IUPAC)
    outfile <- stringr::str_c(
      dir, '/',
      BSgenome::commonName(genome), '-',
      BSgenome::provider(genome), '-',
      BSgenome::providerVersion(genome), '-', IUPAC, '.csv')

    # Add PAM to guides from appropriate guide column
    guides_1_PAM <- stringr::str_c(unique(toupper(na.omit(guides[[column]]))), IUPAC)

    # Filter out guides with ambiguity in fixed region
    guides_1_PAM_good <- guides_1_PAM[!stringr::str_detect(stringr::str_sub(guides_1_PAM, start = 9, end = 20), 'N')]

    iSTOP::search_off_target(
      guides_1_PAM_good,
      genome,
      chromosomes,
      fixed_start  = 9,
      fixed_end    = 20,
      max_mismatch = 2,
      secondary_filter = regex,
      cores = cores # !!! Change cores if necessary on your machine !!!
    ) %>%
      group_by(guide) %>%
      summarise(n_fuzzy_match = n()) %>%
      ungroup %>%
      arrange(desc(n_fuzzy_match)) %>%
      write_csv(outfile)
  })
}

# ---- Hsapiens ----
guides <- readr::read_csv('data/iSTOP-compact/Hsapiens-hg38.csv', col_types = cols(percent_tx = col_double(), percent_NMD = col_double(), rel_pos_largest_isoform = col_double()))
genome <- BSgenome.Hsapiens.UCSC.hg38::Hsapiens
chromosomes <- BSgenome::seqnames(genome)[1:24]

count_off_target(dir = 'data/Off-target-counts', guides, genome, chromosomes, PAMs)

# ---- Athaliana ----
guides <- readr::read_csv('data/iSTOP-compact/Athaliana-plantsmart28.csv', col_types = cols(percent_tx = col_double(), percent_NMD = col_double(), rel_pos_largest_isoform = col_double()))
genome <- BSgenome.Athaliana.TAIR.TAIR9::Athaliana
chromosomes <- BSgenome::seqnames(genome)[1:5]

count_off_target(dir = 'data/Off-target-counts', guides, genome, chromosomes, PAMs)

# ---- Celegans ----
guides <- readr::read_csv('data/iSTOP-compact/Celegans-ce11.csv', col_types = cols(percent_tx = col_double(), percent_NMD = col_double(), rel_pos_largest_isoform = col_double()))
genome <- BSgenome.Celegans.UCSC.ce11::Celegans
chromosomes <- BSgenome::seqnames(genome)[1:6]

count_off_target(dir = 'data/Off-target-counts', guides, genome, chromosomes, PAMs)

# ---- Dmelanogaster ----
guides <- readr::read_csv('data/iSTOP-compact/Dmelanogaster-dm6.csv', col_types = cols(percent_tx = col_double(), percent_NMD = col_double(), rel_pos_largest_isoform = col_double()))
genome <- BSgenome.Dmelanogaster.UCSC.dm6::Dmelanogaster
chromosomes <- BSgenome::seqnames(genome)[1:7]

count_off_target(dir = 'data/Off-target-counts', guides, genome, chromosomes, PAMs)

# ---- Mmusculus ----
guides <- readr::read_csv('data/iSTOP-compact/Mmusculus-mm10.csv', col_types = cols(percent_tx = col_double(), percent_NMD = col_double(), rel_pos_largest_isoform = col_double()))
genome <- BSgenome.Mmusculus.UCSC.mm10::Mmusculus
chromosomes <- BSgenome::seqnames(genome)[1:21]

count_off_target(dir = 'data/Off-target-counts', guides, genome, chromosomes, PAMs)

# ---- Rnorvegicus ----
guides <- readr::read_csv('data/iSTOP-compact/Rnorvegicus-rn6.csv', col_types = cols(percent_tx = col_double(), percent_NMD = col_double(), rel_pos_largest_isoform = col_double()))
genome <- BSgenome.Rnorvegicus.UCSC.rn6::Rnorvegicus
chromosomes <- BSgenome::seqnames(genome)[1:22]

count_off_target(dir = 'data/Off-target-counts', guides, genome, chromosomes, PAMs)

# ---- Scerevisiae ----
guides <- readr::read_csv('data/iSTOP-compact/Scerevisiae-sacCer3.csv', col_types = cols(percent_tx = col_double(), percent_NMD = col_double(), rel_pos_largest_isoform = col_double()))
genome <- BSgenome.Scerevisiae.UCSC.sacCer3::Scerevisiae
chromosomes <- BSgenome::seqnames(genome)[1:16]

count_off_target(dir = 'data/Off-target-counts', guides, genome, chromosomes, PAMs)


# ---- Drerio ----
guides <- readr::read_csv('data/iSTOP-compact/Drerio-danRer10.csv', col_types = cols(percent_tx = col_double(), percent_NMD = col_double(), rel_pos_largest_isoform = col_double()))
genome <- BSgenome.Drerio.UCSC.danRer10::Drerio
chromosomes <- BSgenome::seqnames(genome)[1:25]

count_off_target(dir = 'data/Off-target-counts', guides, genome, chromosomes, PAMs)



