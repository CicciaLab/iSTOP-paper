# ---- Figure 1A: Number of codons that can be targeted with CAS9 base editing ----
Fig1A_data <- read_csv('data/Figure-data/Figure-1A.csv')

Fig1A <-
  Fig1A_data %>%
  mutate(group = forcats::fct_rev(group)) %>%
  ggplot(aes(x = strand, y = n / 64, fill = group)) +
  geom_col(color = 'black', alpha = 0.5, width = 0.4) +
  geom_text(aes(label = n), position = 'stack', vjust = 1.4, fontface = 'bold') +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  scale_fill_manual(values = c('white', 'grey', '#e55b4a')) +
  labs(y = '% of codons') +
  guides(fill = guide_legend(title = NULL)) +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank()
  )

Fig1A
ggsave('figures/Figure-1A.pdf', Fig1A, width = 4, height = 3)

# ---- Figure 1B: Number of amino acids that can be mutated to/from with CAS9 base editing ----
Fig1D_data <- read_csv('data/Figure-data/Figure-1D.csv')
AA_order <- c(
  "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His",
  "Ile", "Met", "Leu", "Lys", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val", "Stop")
AA_order <- rev(AA_order)
AA_color <- c('red', rep('black', times = 20))

Fig1D <-
  Fig1D_data %>%
  mutate(
    AA = ordered(AA, levels = AA_order),
    n_AA = ifelse(group == 'Mutated', n_AA * -1, n_AA)
  ) %>%
  ggplot(aes(x = AA, y = n_AA, fill = group)) +
  geom_col(position = 'identity', color = 'black', width = 0.6) +
  geom_hline(yintercept = 0) +
  scale_y_continuous(breaks = -6:3, expand = c(0, 0), labels = c(6:0, 1:3)) +
  labs(y = 'Number of amino acids') +
  scale_fill_manual(values = c('#d02c1e', '#77cd99'), labels = c('To', 'From')) +
  guides(fill = guide_legend(title = 'Base edit', label.position = 'right', direction = 'vertical')) +
  theme_bw() +
  theme(
    axis.text.x = element_text(color = AA_color, angle = 45, face = 'bold', hjust = 1),
    axis.title.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.025, 0.05),
    legend.justification = c(0, 0),
    legend.background = element_rect(color = 'black')
  ) +
  coord_fixed(ylim = c(-7, 4))
Fig1D
ggsave('figures/Figure-1D.pdf', Fig1D, height = 3, width = 5)
