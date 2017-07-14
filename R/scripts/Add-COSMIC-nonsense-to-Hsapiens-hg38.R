library(tidyverse)

compact <- read_csv('~/Desktop/iSTOP-paper-master/data/iSTOP-compact/Hsapiens-hg38.csv')
#compact <- read_csv('data/iSTOP-compact/Hsapiens-hg38.csv')

nonsense <-
  read_csv('~/Desktop/iSTOP-paper-master/data/COSMIC/COSMIC-nonsense.csv') %>%
  filter(iSTOP) %>%
  select(gene, chr, strand, genome_coord, cancer_type) %>%
  group_by(gene, chr, strand, genome_coord) %>%
  summarise(cancer_type = stringr::str_c(unique(cancer_type), collapse = ' | ')) %>%
  ungroup %>%
  distinct %>%
  mutate(COSMIC_nonsense = TRUE)

result <- left_join(compact, nonsense)
write_csv(result, '~/Desktop/iSTOP-paper-master/data/COSMIC/compact-COSMIC-Hsapiens-hg38.csv')
