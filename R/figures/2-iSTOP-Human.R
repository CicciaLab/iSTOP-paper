# ---- Number of targetable sites by codon/PAM ----
codons_tx <- read_csv('data/iSTOP-by-codon/Hsapiens-hg38.csv', col_types = cols(aa_coord = col_double()))

# Targetable codons are defined by codon, chromosome, strand, and genome coordinate
# Any given site may be present in multiple transcripts/genes
codons <- codons_tx %>%
  select(codon, chr, strand, genome_coord, starts_with('match')) %>%
  distinct %>%
  filter(!is.na(match_any)) %>% # Not concerned about untargetable genes/transcripts
  rename(match_Targetable = match_any)

codon_summary_data <- codons %>%
  bind_rows(codons %>% mutate(codon = 'CAG\nTGG\nCAA\nCGA')) %>%
  select(-chr, -genome_coord, -strand) %>%
  mutate('_Total' = TRUE) %>%
  group_by(codon) %>%
  summarise_all(sum) %>%
  gather(match_type, count, -codon) %>%
  mutate(
    codon = ordered(codon, levels = c('CGA', 'CAA', 'TGG', 'CAG', 'CAG\nTGG\nCAA\nCGA')),
    match_type = str_extract(match_type, '(?<=_).*'),
    match_type = str_replace(match_type, 'sg', ''),
    match_type = ordered(match_type, levels = c('NGCG', 'NNGRRT', 'NGAG', 'NNNRRT', 'NGNG', 'NGG', 'NGA', 'Targetable', 'Total'))
  )

codon_summary_human <-
  codon_summary_data %>%
  ggplot(aes(x = codon, y = count / 1000, color = match_type, fill = match_type)) +
  geom_col(aes(alpha = match_type), position = 'dodge') +
  scale_y_continuous(breaks = seq(0, 900, by = 100)) +
  scale_alpha_discrete(range = c(0.2, 0.8)) +
  coord_flip(ylim = c(0, 900), expand = F) +
  guides(alpha = guide_legend(title = '', reverse = T), fill = guide_legend(title = '', reverse = T), color = guide_legend(title = '', reverse = T)) +
  theme_bw() +
  geom_vline(xintercept = seq(1.5, 4.5, by = 1), size = 0.5, color = 'grey') +
  labs(y = 'Count (x1000)') +
  theme(
    axis.title.y = element_blank(),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 'cm'),
    panel.grid.major.y = element_blank(),
    aspect.ratio = 1
  )
codon_summary_human

ggsave('figures/Codon-summary-human.pdf', codon_summary, height = 5, width = 5)

# ---- Protein length vs. Targetable codons ----
match_counts_by_pep_length_data <- codons_tx %>%
  select(tx, pep_length, codon, starts_with('match_')) %>%
  group_by(tx, pep_length) %>%
  summarise(
    n_codons   = length(which(!is.na(codon))), # number of non NA codon rows is number of codons
    n_match    = length(which(match_any)),
    n_sgNGG    = length(which(match_sgNGG)),
    n_sgNGA    = length(which(match_sgNGA)),
    n_sgNGCG   = length(which(match_sgNGCG)),
    n_sgNGAG   = length(which(match_sgNGAG)),
    n_sgNNGRRT = length(which(match_sgNNGRRT)),
    n_sgNNNRRT = length(which(match_sgNNNRRT))
  )

match_counts_by_pep_length <-
  match_counts_by_pep_length_data %>%
  ggplot(aes(x = pep_length, y = n_match)) +
  geom_point(alpha = 0.2) +
  geom_point(alpha = 0.5, color = 'red', data = filter(match_counts_by_pep_length_data, n_match == 0L)) +
  geom_density2d(color = '#ffe371') +
  scale_x_log10(breaks = c(10, 100, 1000, 10000, 100000), labels = c(10, 100, '1,000', '10,000', '100,000')) +
  scale_y_log10(breaks = c(0.5, 1, 10, 100, 1000, 2000), labels = c(0, 1, 10, 100, '1,000', '')) +
  annotation_logticks() +
  stat_smooth(method = 'lm', color = '#99b1df', se = F) +
  coord_cartesian(xlim = c(10, 100000), ylim = c(0.5, 2000), expand = F) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), aspect.ratio = 1) +
  annotate(
    'text',
    color = 'red',
    hjust = 'left',
    x = 1050,
    y = 2,
    label = str_c(
      nrow(filter(match_counts_by_pep_length_data, n_match == 0L)),
      ' (',
      round((nrow(filter(match_counts_by_pep_length_data, n_match == 0L)) / nrow(match_counts_by_pep_length_data)) * 100, digits = 2),
      '%) untargetable')
  ) +
  annotate('segment', x = 1000, y = 1.6, xend = 150, yend = 0.6, size = 1, color = 'red', arrow = arrow(angle = 15, length = unit(0.15, "inches"), type = 'closed')) +
  labs(x = 'Protein length (AA)', y = 'Targetable codons')
match_counts_by_pep_length

ggsave('figures/Match-counts-by-pep-length.pdf', match_counts_by_pep_length, height = 5, width = 5)

# ---- ECDF earliest targetable codon ----

# Currently only used for order (not color). Uncomment scale_color_manual
colors <- c(
  Total  = 'red',
  All    = '#242424',
  NGA    = '#faa11f',
  NGG    = '#009f77',
  NNNRRT = '#f4e449',
  NGAG   = '#0078b5',
  NNGRRT = '#f16424',
  NGCG   = '#e27fad'
)

gene_tx_codon <-
  codons_tx %>%
  mutate(relative_position = aa_coord / pep_length)

earliest_by_tx <-
  gene_tx_codon %>%
  select(tx, genome_coord, pep_length, aa_coord, relative_position, starts_with('match_')) %>%
  group_by(tx, pep_length) %>%
  summarise(
    untargetable                      = all(is.infinite(relative_position)),
    earliest_relative_position_All    = min(relative_position[match_any]),
    earliest_relative_position_NGG    = min(relative_position[match_sgNGG]),
    earliest_relative_position_NGA    = min(relative_position[match_sgNGA]),
    earliest_relative_position_NGCG   = min(relative_position[match_sgNGCG]),
    earliest_relative_position_NGAG   = min(relative_position[match_sgNGAG]),
    earliest_relative_position_NNGRRT = min(relative_position[match_sgNNGRRT]),
    earliest_relative_position_NNNRRT = min(relative_position[match_sgNNNRRT])
  ) %>%
  ungroup %>%
  # if not targetable then set relative position to 2 (outside range of [0, 1])
  mutate_at(vars(starts_with('earliest')), funs(ifelse(is.infinite(.), 2, .)))


earliest_by_tx_gathered <-
  earliest_by_tx %>%
  gather(PAM, earliest_relative_position, starts_with('earliest_relative_position')) %>%
  mutate(PAM = str_extract(PAM, '(?<=(position_)).*')) %>%
#  bind_rows(data_frame(PAM = 'Total', earliest_relative_position = 2)) %>%  # used only to get matching color scale with codon_summary_human
  mutate(PAM = PAM %>% ordered(levels = rev(names(colors))))

ecdf_earliest_targetable_position <-
  ggplot(earliest_by_tx_gathered) +
  stat_ecdf(
    aes(earliest_relative_position, color = PAM),
    size  = 1.75,
    alpha = 0.75
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = F) +
  #scale_color_manual(values = colors) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), labels = scales::percent_format()) +
  labs(x = 'Relative position in peptide',
       y = 'Earliest targetable position < relative position \n(% of protein coding transcripts)') +
  theme_bw() +
  theme(aspect.ratio = 1, panel.grid.minor = element_blank())
ecdf_earliest_targetable_position

ggsave('figures/ECDF-earliest-targetable-position.pdf', ecdf_earliest_targetable_position, height = 5, width = 5)

# ---- Targetable homologs ----
targetable_homologs_data <-
  read_csv('data/Figure-data/Targetable-Homologs.csv') %>%
  summarise(
    'Untargetable in Human' = n(),
    'Mouse homolog'     = length(which(!is.na(mouse_ucsc))),
    'Rat homolog'       = length(which(!is.na(rat_gene))),
    'Fish homolog'      = length(which(!is.na(fish_gene))),
    'Worm homolog'      = length(which(!is.na(worm_gene))),
    'Fly homolog'       = length(which(!is.na(fly_gene))),
    'Yeast homolog'     = length(which(!is.na(yeast_gene))),
    'Combined homologs' = length(which(
      (!is.na(rat_gene)) |
        (!is.na(mouse_ucsc)) |
        (!is.na(fish_gene)) |
        (!is.na(worm_gene)) |
        (!is.na(fly_gene)) |
        (!is.na(yeast_gene))
    ))
  ) %>%
  gather(type, count) %>%
  mutate(type = ordered(type, levels = type[order(count)])) %>%
  mutate(group = ifelse(type == 'Untargetable in Human', 'Untargetable', 'Targetable'))

targetable_homologs <-
  targetable_homologs_data %>%
  ggplot(aes(y = count, x = type, fill = group)) +
  geom_col(color = 'black', position = 'identity', width = 0.6) +
  geom_text(aes(label = count, color = group), fontface = 'bold', hjust = 1.2) +
  scale_y_continuous(breaks = seq(0, 70, by = 10), expand = c(0, 0), labels = scales::comma) +
  scale_fill_manual(values = c('#ce9bf9', 'white')) +
  scale_color_manual(values = c('white', 'black')) +
  coord_flip(ylim = c(0, 70)) +
  labs(x = NULL, y = 'Genes') +
  guides(fill = guide_legend(title = NULL, reverse = T), color = 'none') +
  theme_bw() +
  theme(
    axis.text.y = element_text(face = 'plain', color = 'black'),
    plot.margin = unit(c(1,0.5,1,0.5),"cm"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = 'black', linetype = 'dotted'),
    legend.position = c(0.95, 0.05),
    legend.justification = c(1, 0),
    legend.direction = 'vertical',
    legend.background = element_rect(colour = 'black', fill = 'white')
  )
targetable_homologs

ggsave('figures/Targetable-homologs.pdf', targetable_homologs, width = 5, height = 5)

# ---- Restriction enzymes
RFLP <- 'data/iSTOP-compact/Hsapiens-hg38.csv' %>%
  read_csv(col_types = cols_only(
    chr          = col_character(),
    strand       = col_character(),
    genome_coord = col_integer(),
    RFLP_150     = col_character(),
    RFLP_100     = col_character(),
    RFLP_50      = col_character())
  ) %>%
  distinct

RFLP_summary_data <-
  sites %>%
  summarise(
    total = n(),
    RFLP_150 = length(which(!is.na(RFLP_150))),
    RFLP_100 = length(which(!is.na(RFLP_100))),
    RFLP_50  = length(which(!is.na(RFLP_50)))
  ) %>%
  gather(distance, count, -total) %>%
  mutate(
    distance = gsub('RFLP_', '', distance) %>% as.integer,
    percent = count / total
  )

RFLP_summary <-
  RFLP_summary_data %>%
  ggplot(aes(y = percent, x = distance)) +
  geom_col(aes(fill = 'orange', color = 'orange'), alpha = 0.5) +
  #geom_line(aes(color = 'orange'), size = 2) +
  #geom_point(aes(color = 'orange'), size = 5) +
  #geom_point(aes(color = 'white'), size = 2) +
  scale_color_manual(breaks = c('orange', 'white'), values = c('#f2a880', 'white')) +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  scale_x_continuous(breaks = c(50, 100, 150)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(y = 'iSTOP verifiable by RFLP\n(loss of cutting after edit)', x = 'Range of unique cutting\n(+/- bases from edited site)') +
  theme_bw() +
  theme(panel.grid.major.y = element_line(color = 'black', linetype = 'dotted')) +
  guides(color = 'none', fill = 'none')
RFLP_summary

ggsave('figures/RFLP-summary.pdf', RFLP_summary, width = 5, height = 5)

# ---- Combining into a single figure ----
#gt_Fig3A <- Fig3A %>% ggplot_build %>% ggplot_gtable
#gt_Fig3B <- Fig3B %>% ggplot_build %>% ggplot_gtable
#gt_Fig3C <- Fig3C %>% ggplot_build %>% ggplot_gtable
#gt_Fig3D <- Fig3D %>% ggplot_build %>% ggplot_gtable
#gt_Fig3E <- Fig3E %>% ggplot_build %>% ggplot_gtable
#gt_Fig3F <- Fig3F %>% ggplot_build %>% ggplot_gtable
#
#gt_Fig3AB <- cbind(gt_Fig3A, gt_Fig3B)
#gt_Fig3CD <- cbind(gt_Fig3C, gt_Fig3D)
#gt_Fig3EF <- cbind(gt_Fig3E, gt_Fig3F)
#gt_Fig3   <- rbind(gt_Fig3AB, gt_Fig3CD, gt_Fig3EF)
#
#pdf('figures/Figure-3.pdf', width = 10, height = 12)
#plot(gt_Fig3) # Too big! Go home!
#dev.off()

