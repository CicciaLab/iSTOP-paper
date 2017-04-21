library(tidyverse)

iSTOP  <- read_csv('data/iSTOP/Hsapiens-hg38.csv')
n_tx   <- length(unique(iSTOP$tx))
n_gene <- length(unique(iSTOP$gene))
codons <- read_csv('data/iSTOP-by-codon/Hsapiens-hg38.csv', col_types = cols(aa_coord = col_double()))
any(is.na(codons$genome_coord))
unique_codons <-
  codons %>%
  select(codon, chr, strand, genome_coord, starts_with('match')) %>%
  filter(!is.na(genome_coord)) %>%
  distinct()
unique_codons
unique_codons %>%
  filter(match_any)

tx_summary <-
  iSTOP %>%
  filter(!is.na(genome_coord)) %>%
  group_by(tx) %>%
  summarise(earliest_aa = min(aa_coord), NMD = any(NMD_pred))

earliest_aa_100 <- tx_summary %>%
  filter(earliest_aa <= 100)

nrow(earliest_aa_100) / n_tx

NMD <- tx_summary %>% filter(NMD)
nrow(NMD) / n_tx

rflp <- read_csv('data/iSTOP-compact/Hsapiens-hg38.csv')
rflp_gene <- rflp %>% filter(!is.na(RFLP_50)) %>% select(gene) %>% distinct()
nrow(rflp_gene) / n_gene
0.625 * 0.04 * 836715

all_sgRNA <-
  list.files('data/iSTOP-compact', '[.]csv', full.names = T) %>%
  map_df(~read_csv(., col_types = cols(percent_tx = col_double())) %>% select(gene, chr, strand, genome_coord) %>% distinct %>% filter(!is.na(genome_coord)))

all_sgRNA
