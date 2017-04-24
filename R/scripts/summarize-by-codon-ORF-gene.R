# ---- Codons ----
files <- list.files('data/iSTOP-by-codon', '[.]csv', full.names = T)

compute_sites_per_100aa <- function(path) {
  path %>%
    read_csv(col_types = cols(aa_coord = col_double())) %>%
    select(tx, pep_length, match_any) %>%
    mutate(match_any = if_else(is.na(match_any), FALSE, match_any)) %>%
    group_by(tx, pep_length) %>%
    summarise(n_matches = sum(match_any)) %>%
    ungroup %>%
    summarise(sites_per_100aa = mean(n_matches / (pep_length / 100))) %>%
    mutate(species = str_replace(basename(path), '-.*', ''))
}

compute_codon_summary <- function(path) {
  path %>%
    read_csv(col_types = cols(aa_coord = col_double())) %>%
    select(codon, chr, strand, genome_coord, match_any) %>%
    filter(!is.na(match_any)) %>% # only counting codons
    distinct %>%
    summarise(
      Total = n(),
      Targetable = length(which(match_any)),
      `Not Targetable` = length(which(!match_any)),
      `% Targetable` = (Targetable / n()) * 100
    ) %>%
    mutate(species = str_replace(basename(path), '-.*', ''))
}

codon_summary <- map_df(files, compute_codon_summary)
sites_summary <- map_df(files, compute_sites_per_100aa)

data <- left_join(codon_summary, sites_summary, by = 'species')

write_csv(data, 'data/Figure-data/Codon-Summary.csv')

# ---- ORFs ----

files <- list.files('data/iSTOP', '[.]csv', full.names = T)

compute_targetable_transcripts <- function(path) {
  #  path <- 'data/iSTOP/Hsapiens-hg38.csv'
  path %>%
    read_csv(col_types = cols(percent_tx = col_double())) %>%
    select(tx, match_any, genome_coord) %>%
    mutate(
      no_codons = is.na(genome_coord),
      match_any = if_else(no_codons, FALSE, match_any)
    ) %>%
    group_by(tx, no_codons) %>%
    summarise(match_any = any(match_any)) %>%
    ungroup %>%
    summarise(
      Total = n(),
      Targetable       = length(which(match_any)),
      `Not Targetable` = length(which(!match_any)),
      `No Codons`      = length(which(no_codons)),
      `% Targetable` = (Targetable / n()) * 100
    ) %>%
    mutate(species = str_replace(basename(path), '-.*', ''))
}

data <- map_df(files, compute_targetable_transcripts)

write_csv(data, 'data/Figure-data/ORF-Summary.csv')

# ---- Genes ----
files <- list.files('data/iSTOP', '[.]csv', full.names = T)

compute_targetable_genes <- function(path) {
  #  path <- 'data/iSTOP/Hsapiens-hg38.csv'
  path %>%
    read_csv(col_types = cols(percent_tx = col_double())) %>%
    select(gene, tx, match_any, genome_coord) %>%
    mutate(
      no_codons = is.na(genome_coord),
      match_any = if_else(no_codons, FALSE, match_any)
    ) %>%
    group_by(gene, tx) %>%
    summarise(targetable = any(match_any)) %>%
    group_by(gene) %>%
    summarise(n_tx = n(), n_targetable = sum(targetable)) %>%
    ungroup %>%
    summarise(
      Total = n(),
      one_or_more_targetable = length(which(n_targetable >= 1)),
      all_targetable = length(which(n_targetable == n_tx))
    ) %>%
    mutate(species = str_replace(basename(path), '-.*', ''))
}

data <- map_df(files, compute_targetable_genes)

write_csv(data, 'data/Figure-data/Gene-Summary.csv')
