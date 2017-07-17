RFLP <-
  read_csv('data/iSTOP-compact/Hsapiens-hg38.csv') %>%
  select(chr, strand, genome_coord, RFLP_C_150, RFLP_T_150) %>%
  distinct()

enzymes <- read_csv(system.file('db/restriction-enzymes.csv', package = 'iSTOP'))$enzyme %>% unique %>% set_names(.)

counts_per_enzyme_C <-
  map_int(enzymes, ~length(which(str_detect(RFLP$RFLP_C_150, str_c('^', ., '$|', '^', ., ' | ', ., '$| ', . ,' ')) & !is.na(RFLP$RFLP_C_150))))

counts_per_enzyme_T <-
  map_int(enzymes, ~length(which(str_detect(RFLP$RFLP_T_150, str_c('^', ., '$|', '^', ., ' | ', ., '$| ', . ,' ')) & !is.na(RFLP$RFLP_T_150))))

enzyme_performance <-
  data_frame(
    enzyme = enzymes,
    n_loss_cut = counts_per_enzyme_C,
    n_gain_cut = counts_per_enzyme_T
  ) %>%
  mutate(
    n_sites = nrow(RFLP),
    percent_loss_cut = (n_loss_cut / n_sites) * 100,
    percent_gain_cut = (n_gain_cut / n_sites) * 100
  ) %>%
  arrange(-n_loss_cut)

write_csv(select(enzyme_performance, enzyme, n_loss_cut, percent_loss_cut, n_gain_cut, percent_gain_cut), 'data/Figure-data/RFLP-enzyme-performance.csv')

RFLP_loss_enzyme_performance_top10 <- enzyme_performance %>%
  #top_n(50, percent_loss_cut) %>%
  top_n(10, percent_loss_cut) %>%
  mutate(enzyme = ordered(enzyme, levels = enzyme[order(percent_loss_cut)])) %>%
  ggplot(aes(x = enzyme, y = percent_loss_cut)) +
  geom_col(color = 'black', fill = 'orange', alpha = 0.5) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip(ylim = c(0, 3)) +
  labs(y = "% of iSTOP sites with unique cutting +/- 150 bp") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_line(color = 'black', linetype = 'dotted'), panel.grid.major.y = element_blank()) +
  theme(axis.title.y = element_blank())

RFLP_loss_enzyme_performance_top10

ggsave('figures/RFLP-loss-enzyme-performance-top10.png', RFLP_loss_enzyme_performance_top10, width = 5, height = 7)
#ggsave('figures/RFLP-loss-enzyme-performance-top10.pdf', RFLP_loss_enzyme_performance_top10, width = 5, height = 7)

RFLP_gain_enzyme_performance_top10 <- enzyme_performance %>%
  #top_n(50, percent_loss_cut) %>%
  top_n(10, percent_gain_cut) %>%
  mutate(enzyme = ordered(enzyme, levels = enzyme[order(percent_gain_cut)])) %>%
  ggplot(aes(x = enzyme, y = percent_gain_cut)) +
  geom_col(color = 'black', fill = 'orange', alpha = 0.5) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip(ylim = c(0, 3)) +
  labs(y = "% of iSTOP sites with unique cutting +/- 150 bp") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_line(color = 'black', linetype = 'dotted'), panel.grid.major.y = element_blank()) +
  theme(axis.title.y = element_blank())

RFLP_gain_enzyme_performance_top10

ggsave('figures/RFLP-gain-enzyme-performance-top10.png', RFLP_gain_enzyme_performance_top10, width = 5, height = 7)
#ggsave('figures/RFLP-gain-enzyme-performance-top10.pdf', RFLP_gain_enzyme_performance_top10, width = 5, height = 7)

RFLP_enzyme_performance_top50 <- enzyme_performance %>%
  top_n(50, percent_loss_cut) %>%
  #top_n(10, percent_loss_cut) %>%
  mutate(enzyme = ordered(enzyme, levels = enzyme[order(percent_loss_cut)])) %>%
  ggplot(aes(x = enzyme, y = percent_loss_cut)) +
  geom_col(color = 'black', fill = 'orange', alpha = 0.5) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip(ylim = c(0, 5)) +
  labs(y = "% of iSTOP sites with unique cutting +/- 150 bp") +
  theme_bw() +
  theme(axis.title.y = element_blank())

RFLP_enzyme_performance_top50

ggsave('figures/RFLP-enzyme-performance-top50.png', RFLP_enzyme_performance_top50, width = 5, height = 7)

