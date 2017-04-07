# ---- Figure 3A: ----
Fig3A <- ggplot() + theme_bw()

# ---- Figure 3B: Number of targetable sites by codon/PAM ----
codons_tx <- read_csv('data/iSTOP-by-codon/Hsapiens-hg38.csv', col_types = cols(aa_coord = col_double()))

# Targetable sites are defined by chromosome, genome coordinate, and codon
# Any given site may be present in multiple transcripts/genes
codons <- codons_tx %>%
  select(chr, strand, genome_coord, codon, starts_with('match')) %>%
  distinct %>%
  filter(!is.na(match_any)) %>% # Not concerned about untargetable genes/transcripts
  rename(match_Targetable = match_any)

codon_summary <- codons %>%
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

Fig3B <-
  codon_summary %>%
  ggplot(aes(x = codon, y = count / 1000, color = match_type, fill = match_type)) +
  geom_col(aes(alpha = match_type), position = 'dodge') +
  #geom_point(position = position_dodge(width = 1)) +
  scale_y_continuous(breaks = seq(0, 900, by = 100)) +
  scale_alpha_discrete(range = c(0.2, 0.8)) +
  coord_flip(ylim = c(0, 900), expand = F) +
  #scale_y_log10() +
  guides(alpha = guide_legend(title = '', reverse = T), fill = guide_legend(title = '', reverse = T), color = guide_legend(title = '', reverse = T)) +
  theme_bw() +
  #facet_grid(codon ~ .) +
  geom_vline(xintercept = seq(1.5, 4.5, by = 1), size = 0.5, color = 'grey') +
  labs(y = 'Count (x1000)') +
  theme(
    axis.title.y = element_blank(),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), 'cm'),
    panel.grid.major.y = element_blank(),
    aspect.ratio = 1
  )
Fig3B

ggsave('figures/Figure-3B.pdf', Fig3B, height = 5, width = 5)

# ---- Figure 3C: Protein length vs. Targetable codons ----
match_counts_by_tx <- codons_tx %>%
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

Fig3C <- match_counts_by_tx %>%
  ggplot(aes(x = pep_length, y = n_match)) +
  geom_point(alpha = 0.2) +
  geom_point(alpha = 0.5, color = 'red', data = filter(match_counts_by_tx, n_match == 0L)) +
  geom_density2d(color = 'orange') +
  scale_x_log10(breaks = c(10, 100, 1000, 10000, 100000), labels = c(10, 100, '1,000', '10,000', '100,000')) +
  scale_y_log10(breaks = c(0.5, 1, 10, 100, 1000, 2000), labels = c(0, 1, 10, 100, '1,000', '')) +
  annotation_logticks() +
  stat_smooth(method = 'lm', color = 'green', se = F) +
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
      nrow(filter(match_counts_by_tx, n_match == 0L)),
      ' (',
      round((nrow(filter(match_counts_by_tx, n_match == 0L)) / nrow(match_counts_by_tx)) * 100, digits = 2),
      '%) untargetable')
  ) +
  annotate('segment', x = 1000, y = 1.6, xend = 150, yend = 0.6, size = 1, color = 'red', arrow = arrow(angle = 15, length = unit(0.15, "inches"), type = 'closed')) +
  labs(x = 'Protein length (AA)', y = 'Targetable codons')
Fig3C

ggsave('figures/Figure-3C.pdf', Fig3C, height = 5, width = 5)

# ---- Figure 3D: ECDF earliest targetable codon ----

# Currently only used for order (not color). Uncomment scale_color_manual
colors <- c(
  Total  = 'red',
  All    = '#242424', # Black
  NGA    = '#faa11f', # Orange
  NGG    = '#009f77', # Bluish Green
  #NGNG   = '#00b8ea', # Sky Blue
  NNNRRT = '#f4e449', # Yellow
  NGAG   = '#0078b5', # Blue
  NNGRRT = '#f16424', # Vermillion
  NGCG   = '#e27fad' # Reddish Purple
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
#  bind_rows(data_frame(PAM = 'Total', earliest_relative_position = 2)) %>%  # add only to get matching color scale with Fig 3B
  mutate(PAM = PAM %>% ordered(levels = rev(names(colors))))

Fig3D <-
  ggplot(earliest_by_tx_gathered) +
  #geom_hline(yintercept = c(0, 1)) +
  #geom_hline(yintercept = frac_targetable, linetype = 'dotted') +
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
Fig3D

ggsave('figures/Figure-3D.pdf', Fig3D, height = 5, width = 5)

# ---- Figure 3E: Targetable homologs ----
Fig3E_data <-
  read_csv('data/Figure-data/Figure-3E.csv') %>%
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

Fig3E <-
  Fig3E_data %>%
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
Fig3E

ggsave('figures/Figure-3E.pdf', Fig3E, width = 5, height = 5)

# ---- Figure 3F: Restriction enzymes
sites <- 'data/iSTOP-compact/Hsapiens-hg38.csv' %>%
  read_csv(col_types = cols_only(
    chr          = col_character(),
    strand       = col_character(),
    genome_coord = col_integer(),
    RFLP_150     = col_character(),
    RFLP_100     = col_character(),
    RFLP_50      = col_character())
  ) %>%
  distinct

sites_summary <-
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

Fig3F <-
  sites_summary %>%
  ggplot(aes(y = percent, x = distance)) +
  geom_line(aes(color = 'orange'), size = 2) +
  geom_point(aes(color = 'orange'), size = 5) +
  geom_point(aes(color = 'white'), size = 2) +
  scale_color_manual(breaks = c('orange', 'white'), values = c('#f2a880', 'white')) +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  scale_x_continuous(breaks = c(50, 100, 150)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(y = 'iSTOP verifiable by RFLP\n(loss of cutting after edit)', x = 'Range of unique cutting\n(+/- bases from edited site)') +
  theme_bw() +
  theme(panel.grid.major.y = element_line(color = 'black', linetype = 'dotted'))
Fig3F

ggsave('figures/Figure-3F.pdf', Fig3F, width = 5, height = 5)

# ---- Figure 3: Bringing it all together! Hold on to your butts! ----
gt_Fig3A <- Fig3A %>% ggplot_build %>% ggplot_gtable
gt_Fig3B <- Fig3B %>% ggplot_build %>% ggplot_gtable
gt_Fig3C <- Fig3C %>% ggplot_build %>% ggplot_gtable
gt_Fig3D <- Fig3D %>% ggplot_build %>% ggplot_gtable
gt_Fig3E <- Fig3E %>% ggplot_build %>% ggplot_gtable
gt_Fig3F <- Fig3F %>% ggplot_build %>% ggplot_gtable

gt_Fig3AB <- cbind(gt_Fig3A, gt_Fig3B)
gt_Fig3CD <- cbind(gt_Fig3C, gt_Fig3D)
gt_Fig3EF <- cbind(gt_Fig3E, gt_Fig3F)
gt_Fig3   <- rbind(gt_Fig3AB, gt_Fig3CD, gt_Fig3EF)

pdf('figures/Figure-3.pdf', width = 10, height = 12)
plot(gt_Fig3) # Too big! Go home!
dev.off()

