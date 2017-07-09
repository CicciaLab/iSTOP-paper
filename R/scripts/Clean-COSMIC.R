# All CDS coordinates
CDS <- read_csv('data/CDS/Hsapiens-hg38.csv') %>% select(gene, chr, strand, start, end) %>% distinct

# These are all the contiguous intervals of coding sequence for every gene
CDS_reduced <-
  with(CDS, GenomicRanges::GRanges(chr, IRanges::IRanges(start, end), strand)) %>%
  IRanges::reduce() %>%
  IRanges::as.data.frame() %>%
  mutate(chr = as.character(seqnames)) %>%
  select(chr, strand, start, end) %>%
  fuzzyjoin::genome_left_join(CDS, c('chr', 'start', 'end')) %>%
  filter(strand.x == strand.y) %>%
  select(gene, chr = chr.x, strand = strand.x, start = start.x, end = end.x) %>%
  distinct

# All iSTOP coordinates (sites NOT codons)
iSTOP <- read_csv('data/iSTOP/Hsapiens-hg38.csv') %>%
  select(gene, chr, strand, coord = genome_coord, PAM = match_any) %>%
  filter(!is.na(PAM)) %>%
  mutate(iSTOP = T, mutation_class = 'Nonsense') %>%
  distinct

CDS_reduced_valid_genes <- filter(CDS_reduced, gene %in% unique(iSTOP$gene))

# Download "COSMIC Mutation Data" https://cancer.sanger.ac.uk/cosmic/download
COSMIC <-
  read_tsv(
    'data/COSMIC/CosmicMutantExport.tsv.gz',
    col_types = cols_only(
      'ID_sample' = col_integer(),
      'ID_tumour'  = col_integer(),
      'Primary site' = col_character(),
      'Mutation ID' = col_character(),
      'Primary histology' = col_character(),
      'Mutation zygosity' = col_character(),
      'Mutation Description' = col_character(),
      'Mutation CDS' = col_character(),
      'Mutation AA' = col_character(),
      'GRCh' = col_integer(),
      'Mutation genome position' = col_character(),
      'Mutation strand' = col_character(),
      'FATHMM prediction' = col_character(),
      'FATHMM score' = col_double(),
      'Mutation somatic status' = col_character()
    )
  ) %>%
  # Use hg38 coordinates and only consider frameshift (Indel) or substitutions
  filter(
    GRCh == 38,
    `Mutation Description` %>% str_detect('Frameshift|Substitution')
  ) %>%
  # Extract chromosome coordinates from 'Mutation genome position' string
  mutate(
    chr   = str_c('chr', str_extract(`Mutation genome position`, '.*(?=:)')),
    start = str_extract(`Mutation genome position`, '(?<=:).*(?=-)') %>% as.integer(),
    end   = str_extract(`Mutation genome position`, '(?<=-).*') %>% as.integer()
  )

# Annotate with Gene information from UCSC using genome coordinates
COSMIC_annotated <-
  COSMIC  %>%
  select(
    primary_site      = `Primary site`,
    primary_histology = `Primary histology`,
    #zygosity          = `Mutation zygosity`,
    #somatic           = `Mutation somatic status`,
    mutation_id       = `Mutation ID`,
    mutation_desc     = `Mutation Description`,
    mutation_cds      = `Mutation CDS`,
    mutation_aa       = `Mutation AA`,
    chr,
    strand            = `Mutation strand`,
    start,
    end
  ) %>%
  # Add gene name based on UCSC coordinates of CDS
  fuzzyjoin::genome_left_join(CDS_reduced_valid_genes, by = c('chr', 'start', 'end')) %>%
  filter(strand.x == strand.y) %>% # Strand information should correspond to the gene's strand
  distinct %>%
  mutate(
    cancer_type = case_when(
      str_detect(.$primary_histology, 'melanoma')                      ~ 'Malignant melanoma',
      .$primary_site == 'large_intestine'                              ~ 'Colorectal',
      .$primary_site == 'endometrium'                                  ~ 'Endometrial',
      .$primary_site == 'lung'                                         ~ 'Lung',
      .$primary_site == 'liver'                                        ~ 'Liver',
      .$primary_site == 'skin'                                         ~ 'Non-melanoma skin',
      .$primary_site == 'breast'                                       ~ 'Breast',
      .$primary_site == 'stomach'                                      ~ 'Stomach',
      .$primary_site %in% c('upper_aerodigestive_tract', 'oesophagus') ~ 'Upper aerodigestive',
      .$primary_site == 'haematopoietic_and_lymphoid_tissue'           ~ 'Blood',
      .$primary_site == 'prostate'                                     ~ 'Prostate',
      .$primary_site == 'pancreas'                                     ~ 'Pancreatic',
      .$primary_site == 'urinary_tract'                                ~ 'Bladder',
      .$primary_site == 'kidney'                                       ~ 'Kidney',
      .$primary_histology == 'glioma'                                  ~ 'Glioma',
      .$primary_site == 'ovary'                                        ~ 'Ovarian',
      .$primary_site == 'cervix'                                       ~ 'Cervical',
      .$primary_site == 'thyroid'                                      ~ 'Thyroid',
      .$primary_site == 'bone'                                         ~ 'Bone',
      TRUE ~ 'Other'),
    mutation_class = case_when(
      str_detect(.$mutation_desc, 'Frameshift') ~ 'Frameshift',
      str_detect(.$mutation_desc, 'Missense')   ~ 'Missense',
      str_detect(.$mutation_desc, 'Nonsense')   ~ 'Nonsense',
      str_detect(.$mutation_desc, 'silent')     ~ 'Silent',
      TRUE ~ 'Other'
    )
  ) %>%
  select(gene, cancer_type, mutation_id, mutation_class, mutation_cds, mutation_aa, chr = chr.x, strand = strand.x, start = start.x, end = end.x) %>%
  full_join(iSTOP, by = c('mutation_class', 'gene', 'chr', 'strand', 'start' = 'coord')) %>%
  mutate(
    iSTOP = ifelse(is.na(iSTOP), F, iSTOP), # if no matching coordinate for iSTOP then consider it not an iSTOP site
    PAM   = ifelse(iSTOP, PAM, NA)          # leave PAM information missing if not iSTOP
  ) %>%
  distinct

write_csv(COSMIC_annotated, 'data/COSMIC/COSMIC-iSTOP.csv')

# ---- By cancer type mutation class summary ----
by_cancer <-
  COSMIC_annotated %>%
  select(cancer_type, mutation_class, chr, strand, start, end, PAM, iSTOP) %>%
  filter(!is.na(cancer_type)) %>% # only COSMIC when looking across all cancers
  bind_rows(mutate(., cancer_type = 'All cancers')) %>% # Add an all cancer group
  distinct() %>%
  group_by(cancer_type, mutation_class) %>%
  summarise(
    n_nonsense_in_cancer = n(),
    n_iSTOP_sites_in_cancer = sum(iSTOP),
    n_iSTOP_sites_in_cancer_with_PAM = sum(PAM, na.rm = T)
  )

write_csv(by_cancer, 'data/COSMIC/COSMIC-summary-by-cancer.csv')


# ---- By gene iSTOP codons summary ----
by_cancer_nonsense <- filter(by_cancer, mutation_class == 'Nonsense') %>% select(-mutation_class)
n_iSTOP_sites      <- select(iSTOP, chr, strand, coord) %>% distinct %>% nrow
P_event            <- by_cancer_nonsense$n_iSTOP_sites_in_cancer / n_iSTOP_sites
names(P_event)     <- by_cancer_nonsense$cancer_type

COSMIC_nonsense <-
  COSMIC_annotated %>%
  filter(mutation_class == 'Nonsense', !is.na(cancer_type)) %>%
  select(gene, cancer_type, mutation_aa, chr, strand, genome_coord = start, PAM, iSTOP) %>%
  distinct

write_csv(COSMIC_nonsense, 'data/COSMIC/COSMIC-nonsense.csv')

n_sites <-
  iSTOP %>%
  group_by(gene) %>%
  summarise(n_iSTOP_sites_in_gene = n())

by_gene <-
  COSMIC_nonsense %>%
  select(gene, cancer_type, chr, strand, genome_coord, PAM, iSTOP) %>%
  bind_rows(mutate(., cancer_type = 'All cancers')) %>%
  distinct %>%
  group_by(gene, cancer_type) %>%
  summarise(
    n_nonsense_in_gene_in_cancer             = n(),
    n_iSTOP_sites_in_gene_in_cancer          = sum(iSTOP),
    n_iSTOP_sites_in_gene_in_cancer_with_PAM = sum(PAM, na.rm = T)
  ) %>%
  left_join(n_sites) %>%
  left_join(by_cancer_nonsense) %>%
  group_by(gene, cancer_type) %>%
  mutate(
    p = binom.test(
      n_iSTOP_sites_in_gene_in_cancer, # Successes in gene in cancer type
      n_iSTOP_sites_in_gene,           # Trials in gene
      P_event[unique(cancer_type)],    # Probability of success in cancer type
      alternative = 'greater'          # One-sided test
    )$p.value
  ) %>%
  # For all cancers dataset, highest P for any single cancer subtype
  group_by(gene) %>%
  mutate(
    min_p = min(p[cancer_type != 'All cancers']),
    p = ifelse(cancer_type == 'All cancers', min_p, p)
  ) %>%
  # Multiple test correction across all tests
  ungroup %>%
  mutate(q = p.adjust(p, method = 'fdr'), `-log10(q)` = -log10(q))

write_csv(by_gene, 'data/COSMIC/COSMIC-summary-by-gene.csv')

# Clean up workspace
rm(list = ls())

# ---- By cancer type: Tumor suppressor  and oncogenes ----
#CGC <- read_csv("data/COSMIC/Census_allMon Mar 27 23_57_23 2017.csv")
#TSG <- filter(CGC, `Role in Cancer` == 'TSG')
#by_cancer_TSG <-
#  COSMIC_annotated %>%
#  filter(gene %in% TSG$`Gene Symbol`) %>%
#  select(cancer_type, mutation_class, chr, strand, start, end, PAM, iSTOP) %>%
#  filter(!is.na(cancer_type)) %>% # only COSMIC when looking across all cancers
#  bind_rows(mutate(., cancer_type = 'All cancers')) %>% # Add an all cancer group
#  distinct() %>%
#  group_by(cancer_type, mutation_class) %>%
#  summarise(
#    n_nonsense_in_cancer = n(),
#    n_iSTOP_sites_in_cancer = sum(iSTOP),
#    n_iSTOP_sites_in_cancer_with_PAM = sum(PAM, na.rm = T)
#  )
