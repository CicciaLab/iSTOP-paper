library(tidyverse)
figure = 'figures/ECDF-n-sites-per-gene.pdf'

data <- read_csv('data/iSTOP-website-version/Hsapiens-hg38.csv.gz')

summary <- data %>%
  select(gene, chr, strand, genome_coord, sgNGG, sgNGA, sgNGAG, sgNGCG, sgNNGRRT, sgNNNRRT) %>%
  mutate(sgAll = T) %>%
  gather(PAM, guide, sgNGG:sgAll, na.rm = T) %>%
  mutate(PAM = gsub('sg', '', PAM)) %>%
  group_by(gene, PAM) %>%
  summarise(n = n())


frac3 <- with(filter(summary, PAM == 'All'), sum(n <= 3) / length(n))
median <- median(filter(summary, PAM == 'All')$n)

summary %>%
  mutate(
    PAM = ordered(PAM, levels = c('All', 'NGA', 'NGG', 'NNNRRT', 'NGAG', 'NNGRRT', 'NGCG'))
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
  coord_cartesian(xlim = c(1, 1000), ylim = c(0, 1)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), aspect.ratio = 1) +
  labs(x = 'Number of iSTOP targets', y = 'Cumulative % of genes') +
  ggsave(figure, width = 5, height = 5)

