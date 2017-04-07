library(tidyverse)
library(grid)
library(gridExtra)
library(gtable)

cancer_types <- read_csv('data/COSMIC/COSMIC-summary-by-cancer.csv') %>% filter(mutation_class == 'Nonsense') %>% mutate(frac_nonsense_cancer = n_iSTOP_sites_in_cancer / n_nonsense_in_cancer)
plot_data <-
  read_csv('data/COSMIC/COSMIC-summary-by-gene.csv') %>%
  mutate(
    cancer_type = ordered(cancer_type, levels = unique(c(cancer_types$cancer_type[order(cancer_types$frac_nonsense_cancer)], 'Other', 'All cancers'), fromLast = T)),
    cancer_type = forcats::fct_rev(cancer_type),
    frac_sites = n_iSTOP_sites_in_gene_in_cancer / n_iSTOP_sites_in_gene
  ) %>%
  group_by(cancer_type) %>%
  mutate(
    rank = rank(q),
    frac_nonsense_iSTOP = n_iSTOP_sites_in_gene_in_cancer / n_nonsense_in_gene_in_cancer,
    frac_nonsense_iSTOP_targetable = n_iSTOP_sites_in_gene_in_cancer_with_PAM / n_nonsense_in_gene_in_cancer
  )

keep <-
  plot_data %>%
  group_by(gene) %>%
  filter(cancer_type != 'All cancers', `-log10(q)` > 15 | rank <= 3, frac_sites > 0.05) %>%
  ungroup

#gene_order <- plot_data %>% filter(gene %in% unique(keep$gene)) %>% group_by(gene) %>% summarise(n = sum(frac_sites)) %>% arrange(desc(n)) %>% .$gene
gene_order <- plot_data %>% filter(gene %in% unique(keep$gene), cancer_type == 'All cancers') %>% arrange(desc(frac_nonsense_iSTOP)) %>% .$gene

gg1_data <-
  plot_data %>%
  filter(gene %in% unique(keep$gene)) %>%
  mutate(
    #  cancer_type = ordered(cancer_type, levels = cancer_type_order),
    #cancer_type = forcats::fct_rev(cancer_type),
    gene        = ordered(gene, levels = gene_order),
    #gene        = forcats::fct_rev(gene),
    frac_sites  = ifelse(frac_sites < 0.05, 0, frac_sites)
  )

Fig5C <-
  gg1_data %>%
  ggplot(aes(y = cancer_type, x = gene, color = cancer_type)) +
  geom_point(aes(size = frac_sites, alpha = frac_sites)) +
  scale_alpha_continuous(breaks = c(0.05, 0.25, 0.5, 0.75, 1), range = c(0.2, 1),  labels = scales::percent, limits = c(0, 1)) +
  scale_size_continuous(breaks  = c(0.05, 0.25, 0.5, 0.75, 1), range = c(0.1, 10), labels = scales::percent, limits = c(0, 1)) +
  guides(
    color = 'none',
    #size = 'none',
    #alpha = 'none'
    size  = guide_legend(title = 'Premature stop observed in cancer (% of possible iSTOP sites in gene)', reverse = T, label.position = 'bottom', title.position = 'bottom', label.hjust = 0.5),
    alpha = guide_legend(title = 'Premature stop observed in cancer (% of possible iSTOP sites in gene)', reverse = T, label.position = 'bottom', title.position = 'bottom', label.hjust = 0.5)
  ) +
  #  labs(title = 'Percent of possible deamination sites observed in cancer', subtitle = expression(bold('Source:')~'COSMIC v79 ('*italic('cancer.sanger.ac.uk/cosmic')*')')) +
  #coord_fixed() +
  theme_bw() +
  theme(
    axis.text = element_text(color = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank(),
    legend.position = 'bottom',
    legend.margin = margin(t = -0.7, unit = 'cm'),
    legend.justification = 'left',
    legend.background = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(0,0,0.25,0.25),"cm")
  )
Fig5C
ggsave('figures/Figure-5C.pdf', Fig5C, width = 10, height = 5)


gg2_data <-
  cancer_types %>%
  select(cancer_type, n_nonsense_in_cancer, n_iSTOP_sites_in_cancer, n_iSTOP_sites_in_cancer_with_PAM) %>%
  mutate(
    cancer_type = ordered(cancer_type, levels = unique(c(cancer_types$cancer_type[order(cancer_types$frac_nonsense_cancer)], 'Other', 'All cancers'), fromLast = T)),
    cancer_type = forcats::fct_rev(cancer_type),
    Nonsense  = 1,
    iSTOP     = n_iSTOP_sites_in_cancer / n_nonsense_in_cancer,
    iSTOP_PAM = n_iSTOP_sites_in_cancer_with_PAM / n_nonsense_in_cancer
  )

gg2_data_1 <-
  gg2_data %>%
  select(cancer_type, Nonsense, iSTOP, iSTOP_PAM) %>%
  gather(group, frac, -cancer_type) %>%
  mutate(group = ordered(group, levels = c('Nonsense', 'iSTOP', 'iSTOP_PAM'), labels = c('Nonsense', 'iSTOP', 'iSTOP + PAM')))

gg2_data_2 <-
  gg2_data %>%
  select(cancer_type, count = n_iSTOP_sites_in_cancer_with_PAM, frac = iSTOP_PAM)

Fig5D <-
  gg2_data_1 %>%
  ggplot(aes(x = cancer_type, color = cancer_type, fill = cancer_type)) +
  geom_col(aes(alpha = group, y = frac), position = 'identity') +
  geom_text(aes(y = frac, x = cancer_type, label = scales::comma(count)), inherit.aes = F, color = 'white', hjust = 1.2, data = gg2_data_2) +
  #geom_text(aes(y = frac + gg_right1_data$frac[which(gg_right1_data$group == 'n_nonsense_and_iSTOP_codon_targetable')], x = cancer_type_n, label = scales::comma(count)),inherit.aes = F, color = 'black', hjust = 1.2, data = gg_right1_data %>% filter(group == 'n_anot_targetable')) +
  scale_y_continuous(breaks = c(0.2, 0.4, 0.6, 0.8), expand = c(0, 0), labels = scales::percent) +
  guides(color = 'none', fill = 'none', alpha = guide_legend(title = NULL, reverse = T, label.position = 'bottom')) +
  labs(y = 'Nonsense in cancer') +
  theme_bw() +
  coord_flip(ylim = c(0, 1)) +
  theme(
    axis.title.y = element_blank(),
    #axis.text.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(0,0.25,0.25,0),"cm"),
    legend.position = 'bottom',
    legend.margin = margin(t = 0, unit = 'cm'),
    legend.justification = 'left',
    legend.background = element_blank()
  )
Fig5D

ggsave('figures/Figure-5D.pdf', Fig5D, width = 10, height = 5)

gg3_data <-
  gg1_data %>%
  ungroup %>%
  filter(cancer_type == 'All cancers') %>%
  select(gene, frac_nonsense_iSTOP, frac_nonsense_iSTOP_targetable) %>%
  mutate(
    frac_nonsense    = 1#,
    #iSTOP_cancer     = n_iSTOP_sites_in_gene_in_cancer          / n_iSTOP_sites_in_gene,
    #iSTOP_cancer_PAM = n_iSTOP_sites_in_gene_in_cancer_with_PAM / n_iSTOP_sites_in_gene
  ) %>%
  #select(gene, iSTOP_gene, iSTOP_cancer, iSTOP_cancer_PAM) %>%
  gather(group, frac, -gene) %>%
  mutate(
    group = ordered(group, levels = c('frac_nonsense', 'frac_nonsense_iSTOP', 'frac_nonsense_iSTOP_targetable'), labels = c('Nonsense', 'iSTOP', 'iSTOP + PAM'))
    # group = forcats::fct_rev(group)
  )

gg3_data_text <-
  gg1_data %>%
  ungroup %>%
  filter(cancer_type == 'All cancers') %>%
  select(gene, count = n_iSTOP_sites_in_gene_in_cancer_with_PAM, frac = frac_nonsense_iSTOP_targetable)

Fig5A <-
  gg3_data %>%
  ggplot(aes(x = gene, y = frac, color = 'All cancers', fill = 'All cancers', alpha = group)) +
  geom_col(position = 'identity') +
  geom_text(aes(y = frac, x = gene, label = scales::comma(count)), angle = 0, inherit.aes = F, color = 'white', hjust = 0.5, vjust = 1.2, data = gg3_data_text) +
  #geom_text(aes(y = count + gg_top1_data$count[which(gg_top1_data$group == 'n_nonsense_and_iSTOP_codon_targetable')], x = gene, label = scales::comma(count)),inherit.aes = F, color = 'black', hjust = 0.5, vjust = 1.2, data = gg_top1_data %>% filter(group == 'n_anot_targetable')) +
  #scale_alpha_manual(values = c(0.1, 0.25, 0.5)) +
  #scale_x_continuous(breaks = 1:21) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  #scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000), expand = c(0,0), labels = scales::comma) +
  guides(color = 'none', fill = 'none', alpha = guide_legend(title = NULL, label.position = 'left', label.vjust = 0.5, label.hjust = 0.5)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(y = 'Nonsense in gene') +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    #axis.text.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(0.25,0.25,0.25,0.25),"cm"),
    legend.position = 'left',
    legend.margin = margin(t = 0, unit = 'cm'),
    legend.justification = 'left',
    legend.background = element_blank(),
    legend.text = element_text(angle = 90),
    legend.key.height = unit(0.09, 'npc')
  )

Fig5A
ggsave('figures/Figure-5A.pdf', Fig5A, width = 10, height = 5)

# ---- Upper right panel: Frequent stoppers ----

#`-log10(q)` = ifelse(`-log10(q)` > 35, 35, `-log10(q)`)
plot_data_limit <-
  plot_data %>%
  mutate(`-log10(q)` = ifelse(`-log10(q)` > 35, 35, `-log10(q)`))

Fig5B <-
  plot_data_limit %>%
  ggplot(aes(x = cancer_type, y = `-log10(q)`, color = cancer_type)) +
  geom_point(color = 'grey', position = position_jitter(),   data = filter(plot_data_limit, `-log10(q)` <= 3)) +
  scale_y_continuous(breaks = c(0, 3, 5, 10, 15, 20, 25, 30, 35), expand = c(0, 0), labels = c(0, 3, 5, 10, 15, 20, 25, 30, '>35')) +
  geom_point(alpha = 0.5, position = position_jitter(),      data = filter(plot_data_limit, `-log10(q)` > 3, `-log10(q)` <= 12)) +
  geom_point(alpha = 0.5, position = position_jitterdodge(), data = filter(plot_data_limit, `-log10(q)` > 12)) +
  geom_vline(xintercept = seq(0.5, 22, by = 1)) +
  geom_hline(yintercept = 3, linetype = 'dashed') +
  coord_flip(ylim = c(0, 40)) +
  #coord_cartesian(ylim = c(0, 40)) +
  #  coord_polar() +
  labs(y = expression(-log[10](q))) +
  guides(color = 'none') +
  theme_bw() +
  theme(
    panel.grid.minor.x = element_blank(),
    #panel.grid.major.x = element_blank(),
    axis.title.y = element_blank()
    #axis.text.x = element_text(angle = 45, hjust = 1)
  )
Fig5B
ggsave('figures/Figure-5B.pdf', Fig5B, width = 5, height = 5)

g5A  <- ggplotGrob(Fig5A + guides(alpha = 'none'))
g5B  <- ggplotGrob(Fig5B)
g5AB <- cbind(g5A, g5B)
g5C  <- ggplotGrob(Fig5C)
g5D  <- ggplotGrob(Fig5D)
g5CD <- cbind(g5C, g5D)
g5   <- rbind(g5AB, g5CD)
plot(g5)
ggsave('figures/Figure-5.pdf', g5, width = 15, height = 10)
ggsave('figures/Figure-5-narrow.pdf', g5, width = 7.5, height = 10)
