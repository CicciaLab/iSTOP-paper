RFLP <-
  read_csv('data/iSTOP-compact/Hsapiens-hg38.csv') %>%
  select(chr, strand, genome_coord, RFLP_150) %>%
  distinct()

enzymes <- read_csv(system.file('db/restriction-enzymes.csv', package = 'iSTOP'))$enzyme %>% unique %>% set_names(.)

counts_per_enzyme <-
  map(enzymes, ~length(which(str_detect(RFLP$RFLP_150, str_c('^', ., '$|', '^', ., ' | ', ., '$| ', . ,' ')) & !is.na(RFLP$RFLP_150))))

enzyme_performance <-
  data_frame(
    enzyme = names(counts_per_enzyme),
    n      = unlist(counts_per_enzyme)
  ) %>%
  mutate(
    n_sites = nrow(RFLP),
    percent = (n / n_sites) * 100
  )

write_csv(select(enzyme_performance, enzyme, n, percent) %>% arrange(-n), 'data/Figure-data/enzyme-performance.csv')

FigS3A <- enzyme_performance %>%
  #top_n(50, percent) %>%
  top_n(10, percent) %>%
  mutate(enzyme = ordered(enzyme, levels = enzyme[order(percent)])) %>%
  ggplot(aes(x = enzyme, y = percent)) +
  geom_col(color = 'black', fill = 'orange', alpha = 0.5) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  coord_flip(ylim = c(0, 0.05)) +
  labs(y = "% of iSTOP sites with unique cutting +/- 150 bp") +
  theme_bw() +
  theme(axis.title.y = element_blank())
FigS3A
ggsave('figures/Figure-S3A.svg', FigS3A, width = 5, height = 7)

iSTOP <- read_csv('data/iSTOP/Hsapiens-hg38.csv')

NMD <-
  iSTOP %>%
  select(gene, tx, NMD_pred, match_any) %>%
  group_by(gene, tx) %>%
  summarise(NMD_possible = any(NMD_pred & match_any)) %>%
  mutate(NMD_possible = ifelse(is.na(NMD_possible), F, NMD_possible))

(nrow(NMD) - nrow(NMD[!NMD$NMD_possible,])) / nrow(NMD)


