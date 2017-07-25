library(tidyverse)
figure = 'figures/ECDF-n-sites-per-gene.pdf'

data <- read_csv('data/iSTOP/Hsapiens-hg38.csv', col_types = cols(aa_coord = col_double()))

summary <- data %>%
  select(gene, chr, strand, genome_coord, match_any, starts_with('sgN'), -ends_with('_spacing')) %>%
  rename_at(vars(starts_with('sgN')), funs(gsub('sg', '', .))) %>%
  distinct() %>%
  mutate_at(vars(NGG:NNNRRT), funs(!is.na(.))) %>%
  mutate(All = ifelse(is.na(match_any), FALSE, match_any)) %>%
  gather(PAM, guide, NGG:All) %>%
  group_by(gene, PAM) %>%
  summarise(n = sum(guide))

frac3  <- 0.0349
median <- 23

summary %>%
  mutate(
    PAM = ordered(PAM, levels = c('All', 'NGA', 'NGG', 'NNNRRT', 'NGAG', 'NNGRRT', 'NGCG')),
    n = ifelse(n == 0, 0.1, n) # for compatibility with log scale
  ) %>%
  ggplot(aes(x = n)) +
  stat_ecdf(aes(color = PAM)) +
  scale_x_log10(breaks = c(1, 3, 10, median, 100, 1000), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, frac3, 0.25, 0.5, 0.75, 1), labels = scales::percent) +
  annotation_logticks(sides = 'b') +
  geom_segment(aes(x = median, xend = median, y = -1, yend = 0.5), linetype = 'dotted', data = data_frame()) +
  geom_segment(aes(x = 1, xend = median, y = 0.5, yend = 0.5), linetype = 'dotted', data = data_frame()) +
  geom_segment(aes(x = 3, xend = 3, y = -1, yend = frac3), linetype = 'dotted', data = data_frame()) +
  geom_segment(aes(x = 1, xend = 3, y = frac3, yend = frac3), linetype = 'dotted', data = data_frame()) +
  geom_hline(yintercept = c(0, 1), linetype = 'dotted') +
  scale_color_manual(values = c('#cc9ff9', '#49a7f8', '#67c8cd', '#55b96f', '#91b43e', '#c49532', '#e87c70')) +
  coord_cartesian(xlim = c(1, 1000), ylim = c(0, 1)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), aspect.ratio = 1) +
  labs(x = 'Number of iSTOP targets', y = 'Cumulative % of genes') +
  ggsave(figure, width = 5, height = 5)

