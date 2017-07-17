library(tidyverse)
figure = 'figures/Off-target-distribution.pdf'

data <- read_csv('data/iSTOP-website-version/Hsapiens-hg38.csv.gz')

NGG <-
  data %>%
  select(guide = sgNGG, matches = sgNGG_matches) %>%
  filter(complete.cases(.)) %>%
  mutate(guide = toupper(guide))

NGG %>%
  ggplot(aes(x = matches)) +
  stat_ecdf() +
  scale_x_log10(breaks = c(1, 10, 100, 1000, 10000, 100000), expand = c(0, 0), labels = scales::comma) +
  scale_y_continuous(labels = scales::percent) +
  annotation_logticks(sides = 'b') +
  #geom_segment(aes(x = median, xend = median, y = -1, yend = 0.5), linetype = 'dotted', data = data_frame()) +
  #geom_segment(aes(x = 1, xend = median, y = 0.5, yend = 0.5), linetype = 'dotted', data = data_frame()) +
  #geom_segment(aes(x = 3, xend = 3, y = -1, yend = frac3), linetype = 'dotted', data = data_frame()) +
  #geom_segment(aes(x = 1, xend = 3, y = frac3, yend = frac3), linetype = 'dotted', data = data_frame()) +
  geom_hline(yintercept = c(0, 1), linetype = 'dotted') +
  coord_cartesian(xlim = c(0.9, max(NGG$matches)), ylim = c(0, 1)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), aspect.ratio = 1) +
  labs(
    x = "Number of mapped sites in genome\n(Allowing 2 ambiguities oustide of guide's seed sequence)",
    y = 'Cumulative % of human NGG sgSTOPs'
  ) +
  ggsave(figure, width = 5, height = 5)
